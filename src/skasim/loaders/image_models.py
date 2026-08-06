"""Image-model validation, preview, and CASA injection helpers."""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
import warnings
import astropy.units as u
from astropy.wcs import FITSFixedWarning
from astropy.utils.iers import conf
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from loguru import logger
import re

from ..config import (
    CasaTaylorTermsModelEntry,
    ComponentSkyModelEntry,
    ContinuumIAlphaModelEntry,
    ImgConfig,
    ModelEntry,
    ObsConfig,
    SimConfig,
    SpectralCubeModelEntry,
    StaticStokesMapsModelEntry,
    has_spectral_cube_model,
    spectral_cube_model_entries,
)
from ..imaging import write_fits_preview
from ..manifest import RunContext
from ..runtime import require_casacore


# suppress fits formatting fixes
warnings.simplefilter("ignore", category=FITSFixedWarning)
# suppress polar motion fallback warnings
warnings.filterwarnings("ignore", message=".*polar motions.*")

@dataclass(frozen=True)
class FitsCubeInfo:
    """Small summary of an accepted 3D spectral-cube model."""

    path: Path
    shape: tuple[int, int, int]
    spatial_shape: tuple[int, int]
    unit: str
    n_channels: int
    channel_width_hz: float
    start_frequency_hz: float
    reference_frequency_hz: float


@dataclass(frozen=True)
class FitsImageInfo:
    """Small summary of one accepted FITS image model plane."""

    path: Path
    spatial_shape: tuple[int, int]
    unit: str | None
    celestial_header: dict[str, object]
    center: SkyCoord | None


@dataclass(frozen=True)
class CasaModelProduct:
    """CASA-ready model product generated from one model entry."""

    model_paths: list[Path]
    nterms: int
    reffreq: str
    intermediates: list[Path]
    cube_data: np.ndarray | None = None
    header: fits.Header | None = None
    freq_axis: int | None = None
    model_dir: Path | None = None


def component_model_entries(config: SimConfig) -> list[ComponentSkyModelEntry]:
    return [
        entry for entry in config.models if isinstance(entry, ComponentSkyModelEntry)
    ]


def image_model_entries(config: SimConfig) -> list[ModelEntry]:
    return [
        entry
        for entry in config.models
        if isinstance(
            entry,
            (
                ContinuumIAlphaModelEntry,
                CasaTaylorTermsModelEntry,
                StaticStokesMapsModelEntry,
                SpectralCubeModelEntry,
            ),
        )
    ]


def image_model_center(entries: list[ModelEntry]) -> SkyCoord | None:
    """Return the centre of the first image model with usable celestial WCS."""
    for entry in entries:
        path = primary_model_fits_path(entry)
        if path is None:
            continue
        try:
            info = read_fits_image_info(path)
        except Exception as exc:
            logger.debug(f"image_model_center: failed to read {path}: {exc}")
            continue
        if info.center is not None:
            return info.center
    return None


def primary_model_fits_path(entry: ModelEntry) -> Path | None:
    """Return the representative FITS image for previews and phase-centre inference."""
    if isinstance(entry, ContinuumIAlphaModelEntry):
        return Path(entry.stokes_i).expanduser().resolve()
    if isinstance(entry, StaticStokesMapsModelEntry):
        for value in (entry.stokes_i, entry.stokes_q, entry.stokes_u, entry.stokes_v):
            if value:
                return Path(value).expanduser().resolve()
    return None


def read_fits_image_info(path: Path) -> FitsImageInfo:
    """Read FITS image metadata used by validation and reporting."""
    with fits.open(path) as hdul:
        hdu = hdul[0]
        if hdu.data is None:
            raise ValueError(f"{path} has no image data")
        data = np.asarray(hdu.data).squeeze()
        if data.ndim < 2:
            raise ValueError(f"{path} is not a spatial FITS image")
        spatial_shape = tuple(int(v) for v in data.shape[-2:])
        header = hdu.header.copy()
        unit = header.get("BUNIT")

    try:
        wcs = WCS(header).celestial
        celestial_header = dict(wcs.to_header())
        center_y = (spatial_shape[0] - 1) / 2.0
        center_x = (spatial_shape[1] - 1) / 2.0
        center = wcs.pixel_to_world(center_x, center_y)
        if not isinstance(center, SkyCoord):
            center = None
    except Exception as exc:
        logger.debug(f"read_fits_image_info: WCS construction failed for {path}: {exc}")
        celestial_header = {}
        center = None

    return FitsImageInfo(
        path=path,
        spatial_shape=spatial_shape,
        unit=unit,
        celestial_header=celestial_header,
        center=center,
    )


# canonical aliases accepted for Jy/pixel (per-pixel flux density) inputs.
_ACCEPTED_JY_PER_PIXEL_UNITS = frozenset(
    {
        "jy/pixel",
        "jy pix-1",
        "jy/pix",
        "jy",
        "jy px-1",
        "jy/px",
        "jy pixels-1",
        "jy pixel-1",
        "jy pixel^-1",
        "jy pix^-1",
        "jy px^-1",
    }
)


def validate_continuum_i_alpha(entry: ContinuumIAlphaModelEntry) -> dict:
    """Validate the continuum image contract and return report metadata."""
    stokes_info = read_fits_image_info(Path(entry.stokes_i).expanduser().resolve())
    alpha_info = read_fits_image_info(Path(entry.alpha).expanduser().resolve())
    if stokes_info.spatial_shape != alpha_info.spatial_shape:
        raise ValueError(
            "continuum_i_alpha requires matching spatial dimensions: "
            f"{stokes_info.path} has {stokes_info.spatial_shape}, "
            f"{alpha_info.path} has {alpha_info.spatial_shape}"
        )
    if stokes_info.celestial_header != alpha_info.celestial_header:
        raise ValueError("continuum_i_alpha requires matching celestial WCS.")
    unit = (stokes_info.unit or "").strip().lower()
    if unit not in _ACCEPTED_JY_PER_PIXEL_UNITS:
        raise ValueError(
            f"{stokes_info.path} must declare Jy/pixel-compatible BUNIT; "
            f"found {stokes_info.unit!r}"
        )
    alpha_unit = (alpha_info.unit or "").strip().lower()
    if alpha_unit not in {"", "1", "dimensionless", "none"}:
        raise ValueError(
            f"{alpha_info.path} must be dimensionless; found BUNIT={alpha_info.unit!r}"
        )
    return {
        "stokes_i": str(stokes_info.path),
        "alpha": str(alpha_info.path),
        "spatial_shape": list(stokes_info.spatial_shape),
        "unit": stokes_info.unit,
        "reference_frequency_hz": entry.reference_frequency_hz,
    }


# ---------------------------------------------------------------------------
# spectral cube helpers
# ---------------------------------------------------------------------------


def _find_frequency_axis(header: fits.Header) -> int:
    """Return the 1-based FITS axis index whose CTYPE is FREQ.

    FITS axis k (CTYPEk/NAXISk) corresponds to numpy axis -k (k-th from the end).
    """
    naxis = int(header.get("NAXIS", 3))
    for axis in range(1, naxis + 1):
        ctype = str(header.get(f"CTYPE{axis}", "")).strip().upper()
        if ctype.startswith("FREQ"):
            return axis
    raise ValueError("spectral cube has no FREQ axis in CTYPE1..CTYPEn")


def _fits_axis_to_numpy(axis: int, ndim: int = 3) -> int:
    """Return the NumPy axis index for a 1-based FITS axis number.

    FITS axis ``k`` is stored as the ``ndim - k`` NumPy axis (the ``k``-th axis
    from the slowest-varying end).  For a 3D cube this means:
    ``NAXIS1`` -> axis 2, ``NAXIS2`` -> axis 1, ``NAXIS3`` -> axis 0.
    """
    return ndim - axis


def _squeeze_degenerate_axes(data: np.ndarray, header: fits.Header) -> tuple[np.ndarray, fits.Header]:
    """Drop length-1 axes from a FITS cube, keeping the spectral axis intact.

    Spectral cubes are sometimes stored as 4D (RA, DEC, FREQ, STOKES) with a
    single Stokes axis.  This helper squeezes those degenerate axes and updates
    the header so downstream code sees a standard 3D cube.
    """
    if data.ndim == 3:
        return data, header
    if data.ndim < 3 or data.ndim > 4:
        raise ValueError(f"spectral_cube must be 3D or 4D with a degenerate axis, got ndim={data.ndim}")

    freq_axis = _find_frequency_axis(header)
    single_axes = [axis for axis in range(1, data.ndim + 1) if header[f"NAXIS{axis}"] == 1]

    if not single_axes:
        raise ValueError(
            f"spectral_cube has {data.ndim} dimensions but no degenerate (length-1) axis to squeeze"
        )

    if len(single_axes) > 1:
        raise ValueError(
            f"spectral_cube has multiple degenerate axes {single_axes}; cannot disambiguate"
        )

    squeeze_axis = single_axes[0]
    if squeeze_axis == freq_axis:
        raise ValueError(
            f"spectral_cube frequency axis (FREQ in FITS axis {freq_axis}) has length 1; "
            "cannot squeeze the spectral axis"
        )

    # NumPy axis to drop.
    np_axis = _fits_axis_to_numpy(squeeze_axis, data.ndim)
    data = np.squeeze(data, axis=np_axis)

    # Build a clean 3D header preserving the remaining axes in FITS order.
    new_header = header.copy()
    new_header["NAXIS"] = 3

    old_axes = [a for a in range(1, data.ndim + 2) if a != squeeze_axis]
    for new_axis, old_axis in enumerate(old_axes, start=1):
        for key in ("NAXIS", "CTYPE", "CRPIX", "CRVAL", "CDELT", "CUNIT"):
            old_key = f"{key}{old_axis}"
            new_key = f"{key}{new_axis}"
            if old_key in new_header:
                new_header[new_key] = new_header[old_key]

    # Remove leftover 4th-axis keys and any dangling higher-axis keys.
    for key in list(new_header.keys()):
        match = re.match(r"(NAXIS|CTYPE|CRPIX|CRVAL|CDELT|CUNIT)\d+", key)
        if match:
            axis_num = int(match.group(0)[-1])
            if axis_num > 3:
                del new_header[key]

    return data, new_header


def read_fits_cube_info(path: Path) -> dict:
    """Read metadata from a 3D FITS spectral cube.

    Returns a dict with keys: shape, spatial_shape, unit, n_channels,
    channel_width_hz, start_frequency_hz, reference_frequency_hz.
    """
    with fits.open(path) as hdul:
        hdu = hdul[0]
        if hdu.data is None:
            raise ValueError(f"{path} has no image data")
        data = np.asarray(hdu.data)
        header = hdu.header.copy()
        data, header = _squeeze_degenerate_axes(data, header)
        if data.ndim != 3:
            raise ValueError(f"spectral_cube must be 3D, got ndim={data.ndim}")
        freq_axis = _find_frequency_axis(header)
        n_freq = int(header[f"NAXIS{freq_axis}"])
        freqs = _freq_axis_centres(header, n_freq, axis=freq_axis)
        channel_width_hz = float(np.diff(freqs).mean()) if n_freq > 1 else 0.0
        unit = str(header.get("BUNIT") or "Jy/px").strip().lower()

        # spatial shape in (n_ra, n_dec) order regardless of axis order.
        spatial_shape = [None, None]
        for axis in range(1, 4):
            if axis == freq_axis:
                continue
            ctype = str(header.get(f"CTYPE{axis}", "")).strip().upper()
            if ctype.startswith("RA") or ctype.startswith("GLON"):
                spatial_shape[0] = int(header[f"NAXIS{axis}"])
            elif ctype.startswith("DEC") or ctype.startswith("GLAT"):
                spatial_shape[1] = int(header[f"NAXIS{axis}"])
        if None in spatial_shape:
            raise ValueError("spectral cube spatial axes are not labelled as RA/DEC or GLON/GLAT")

    return {
        "path": path,
        "shape": tuple(int(v) for v in data.shape),
        "spatial_shape": tuple(spatial_shape),
        "unit": unit,
        "n_channels": n_freq,
        "channel_width_hz": channel_width_hz,
        "start_frequency_hz": float(freqs[0]),
        "reference_frequency_hz": float(freqs[n_freq // 2]),
        "freq_axis": freq_axis,
    }


def _freq_axis_centres(header: fits.Header, n_channels: int, axis: int = 3) -> np.ndarray:
    """Return the frequency axis centre positions in Hz for the given 1-based axis."""
    crpix = float(header.get(f"CRPIX{axis}", 1.0))
    crval = float(header[f"CRVAL{axis}"])
    cdelt = float(header[f"CDELT{axis}"])
    cunit = str(header.get(f"CUNIT{axis}") or "").strip().lower()
    if cunit == "mhz":
        crval *= 1e6
        cdelt *= 1e6
    return crval + (np.arange(n_channels) - (crpix - 1.0)) * cdelt


def _reorder_cube_axes(data: np.ndarray, header: fits.Header) -> np.ndarray:
    """Reorder a 3D FITS array so its numpy axes become (freq, dec, ra).

    FITS axis ``k`` (CTYPEk/NAXISk) maps to numpy axis ``ndim - k`` (the
    ``k``-th axis from the slowest-varying end).  This function looks at
    CTYPE1..CTYPE3 and moves the frequency axis to the front, followed by DEC
    and RA.
    """
    ndim = data.ndim
    axis_map: dict[str, int] = {}
    for axis in range(1, ndim + 1):
        ctype = str(header.get(f"CTYPE{axis}", "")).strip().upper()
        if ctype.startswith("FREQ"):
            label = "freq"
        elif ctype.startswith("RA") or ctype.startswith("GLON"):
            label = "ra"
        elif ctype.startswith("DEC") or ctype.startswith("GLAT"):
            label = "dec"
        else:
            raise ValueError(f"spectral cube has unsupported CTYPE{axis}={ctype!r}")
        axis_map[label] = _fits_axis_to_numpy(axis, ndim)

    for label in ("freq", "dec", "ra"):
        if label not in axis_map:
            raise ValueError(f"spectral cube is missing {label} axis")

    target_order = ["freq", "dec", "ra"]
    source_axes = [axis_map[label] for label in target_order]
    return np.moveaxis(data, source_axes, [0, 1, 2])


def _read_ms_phase_center(visibility_path: Path) -> tuple[float, float] | None:
    """Return (ra_deg, dec_deg) from the MS FIELD table, or None if unavailable."""
    try:
        from casacore.tables import table
    except Exception as exc:
        logger.debug(
            f"_read_ms_phase_center: casacore not available for {visibility_path}: {exc}"
        )
        return None
    try:
        with table(str(visibility_path / "FIELD"), ack=False) as tb:
            phase_dir = tb.getcol("PHASE_DIR")
            # PHASE_DIR shape can be (2, nrows) or (2, nrows, npoly)
            if phase_dir.ndim == 2:
                ra_rad = float(phase_dir[0, 0])
                dec_rad = float(phase_dir[1, 0])
            else:
                ra_rad = float(phase_dir[0, 0, 0])
                dec_rad = float(phase_dir[1, 0, 0])
            return float(np.degrees(ra_rad)), float(np.degrees(dec_rad))
    except Exception as exc:
        logger.debug(
            f"_read_ms_phase_center: failed to read PHASE_DIR from {visibility_path}: {exc}"
        )
        return None


def _read_ms_spectral_window(visibility_path: Path) -> tuple[np.ndarray, float] | None:
    """Return (chan_freqs_hz, chan_width_hz) from the MS SPECTRAL_WINDOW table."""
    try:
        from casacore.tables import table
    except Exception as exc:
        logger.debug(
            f"_read_ms_spectral_window: casacore not available for {visibility_path}: {exc}"
        )
        return None
    try:
        with table(str(visibility_path / "SPECTRAL_WINDOW"), ack=False) as tb:
            num_chan = np.asarray(tb.getcol("NUM_CHAN"))
            n_spw = len(num_chan)
            n_chan_first = int(num_chan[0])

            chan_freq = np.asarray(tb.getcol("CHAN_FREQ"))
            chan_width = np.asarray(tb.getcol("CHAN_WIDTH"))

            # CHAN_FREQ shape is usually (n_spw, n_chan) but can be transposed.
            # Choose orientation that matches n_chan of the first SPW.
            if chan_freq.shape[0] == n_spw and chan_freq.shape[1] == n_chan_first:
                freqs = chan_freq[0, :]
                widths = chan_width[0, :]
            elif chan_freq.shape[1] == n_spw and chan_freq.shape[0] == n_chan_first:
                freqs = chan_freq[:, 0]
                widths = chan_width[:, 0]
            else:
                # Fallback: flatten and take the first n_chan_first values
                freqs = chan_freq.flatten()[:n_chan_first]
                widths = chan_width.flatten()[:n_chan_first]

            return freqs, float(np.mean(widths))
    except Exception as exc:
        logger.debug(
            f"_read_ms_spectral_window: failed to read SPECTRAL_WINDOW from {visibility_path}: {exc}"
        )
        return None


def _reorder_casa_image_to_radecfreqstokes(image_path: Path) -> None:
    """Ensure the CASA image table has axis order (RA, DEC, FREQ, STOKES).

    .. deprecated::
        Kept for backward compatibility but no longer used by the spectral-cube
        pipeline, which now uses WSClean -predict.
    """
    try:
        from casatools import image
    except Exception as exc:
        logger.debug(f"_reorder_casa_image_to_radecfreqstokes: casatools not available: {exc}")
        return
    try:
        ia = image()
        ia.open(str(image_path))
        coordsys = ia.coordsys()
        names = [str(n).lower() for n in coordsys.names()]
        order = []
        for target in ("ra", "dec", "freq", "stokes"):
            for idx, name in enumerate(names):
                if name.startswith(target):
                    order.append(idx)
                    break
        if len(order) == 4:
            ia.reorder(order)
        ia.close()
    except Exception as exc:
        logger.debug(f"_reorder_casa_image_to_radecfreqstokes: failed for {image_path}: {exc}")


def _resample_spectral_axis_to_ms_channels(
    freqs_in_hz: np.ndarray,
    data_in: np.ndarray,
    freqs_out_hz: np.ndarray,
    df_out_hz: float,
) -> np.ndarray:
    """Resample a 3D spectral cube to the MS channel grid conserving integrated flux.

    data_in shape is (n_freq, n_y, n_x).  For each output channel, integrate the
    input pixels that overlap the output channel width, then divide by the output
    channel width.  This preserves total flux when the input channels are finer
    than (or misaligned with) the output channels.
    """
    n_in, n_y, n_x = data_in.shape
    n_out = len(freqs_out_hz)
    half_out = abs(df_out_hz) / 2.0

    data_in_r = data_in.reshape(n_in, n_y * n_x)
    data_out_r = np.zeros((n_out, n_y * n_x), dtype=np.float32)

    df_in = np.diff(freqs_in_hz)
    if df_in.size == 0:
        df_in = np.array([abs(df_out_hz)])
    edges_in = np.empty(n_in + 1, dtype=freqs_in_hz.dtype)
    edges_in[0] = freqs_in_hz[0] - df_in[0] / 2.0
    edges_in[1:-1] = freqs_in_hz[:-1] + df_in / 2.0
    edges_in[-1] = freqs_in_hz[-1] + df_in[-1] / 2.0

    for out_idx in range(n_out):
        out_lo = freqs_out_hz[out_idx] - half_out
        out_hi = freqs_out_hz[out_idx] + half_out

        in_lo_idx = max(0, np.searchsorted(edges_in, out_lo, side="right") - 1)
        in_hi_idx = min(n_in, np.searchsorted(edges_in, out_hi, side="right"))
        if in_lo_idx >= in_hi_idx:
            continue

        overlap_lo = np.maximum(edges_in[in_lo_idx:in_hi_idx], out_lo)
        overlap_hi = np.minimum(edges_in[in_lo_idx + 1 : in_hi_idx + 1], out_hi)
        weights = np.clip(overlap_hi - overlap_lo, 0.0, None)
        total_weight = float(weights.sum())
        if total_weight <= 0:
            continue

        band = data_in_r[in_lo_idx:in_hi_idx, :]
        weighted_sum = np.average(band, axis=0, weights=weights) * total_weight
        data_out_r[out_idx, :] = weighted_sum / abs(df_out_hz)

    return data_out_r.reshape(n_out, n_y, n_x)


def validate_spectral_cube(
    entry: SpectralCubeModelEntry,
    obs: ObsConfig,
    img: ImgConfig,
) -> dict:
    """Validate the spectral-cube contract against observation/imaging config."""
    path = Path(entry.cube).expanduser().resolve()
    info = read_fits_cube_info(path)

    _accepted_units = _ACCEPTED_JY_PER_PIXEL_UNITS
    if info["unit"] not in _accepted_units:
        raise ValueError(
            f"{path} must declare Jy/pixel-compatible BUNIT; found {info['unit']!r}"
        )

    if info["spatial_shape"] != (img.pixels, img.pixels):
        raise ValueError(
            f"spectral_cube spatial dimensions {info['spatial_shape']} do not match "
            f"imaging pixels {img.pixels}"
        )

    if obs.bandwidth_mhz is None or obs.n_channels is None or obs.channel_width_mhz is None:
        raise ValueError("observation spectral grid is incomplete")

    obs_bw_hz = obs.bandwidth_mhz * 1e6
    obs_center_hz = obs.frequency_mhz * 1e6
    obs_min_hz = obs_center_hz - obs_bw_hz / 2.0
    obs_max_hz = obs_center_hz + obs_bw_hz / 2.0

    n_channels = info["n_channels"]
    channel_width_hz = info["channel_width_hz"]
    cube_center_hz = info["start_frequency_hz"] + (n_channels - 1) * channel_width_hz / 2.0
    cube_min_hz = cube_center_hz - n_channels * channel_width_hz / 2.0
    cube_max_hz = cube_center_hz + n_channels * channel_width_hz / 2.0

    edge_tol_hz = 0.01 * obs_bw_hz
    if cube_min_hz < obs_min_hz - edge_tol_hz or cube_max_hz > obs_max_hz + edge_tol_hz:
        raise ValueError(
            f"spectral_cube frequency range [{cube_min_hz:.3e}, {cube_max_hz:.3e}] Hz "
            f"extends beyond the observation band "
            f"[{obs_min_hz:.3e}, {obs_max_hz:.3e}] Hz"
        )

    logger.info(
        f"validate_spectral_cube: cube {n_channels} channels x {channel_width_hz:.3e} Hz "
        f"covering [{cube_min_hz:.3e}, {cube_max_hz:.3e}] Hz inside observation band "
        f"[{obs_min_hz:.3e}, {obs_max_hz:.3e}] Hz"
    )

    return {
        "cube": str(path),
        "shape": info["shape"],
        "spatial_shape": list(info["spatial_shape"]),
        "unit": info["unit"],
        "n_channels": n_channels,
        "channel_width_hz": channel_width_hz,
        "reference_frequency_hz": cube_center_hz,
        "frequency_range_hz": [cube_min_hz, cube_max_hz],
    }


def prepare_spectral_cube_for_casa(
    ctx: RunContext,
    entry: SpectralCubeModelEntry,
    index: int,
    report: dict,
) -> CasaModelProduct:
    """Resample a 3D FITS spectral cube to the MS spectral grid.

    Returns the resampled cube array in numpy order (freq, dec, ra) together
    with a 3D FITS header that describes the per-channel model images.  The
    caller is responsible for splitting this into WSClean per-channel model
    images and running ``wsclean -predict``.
    """
    source_path = Path(entry.cube).expanduser().resolve()
    prefix = f"model_entry_{index + 1:02d}_spectral_cube"
    cube_fits = ctx.work_dir / f"{prefix}.fits"

    obs = ctx.config.observation
    if obs.n_channels is None:
        raise ValueError("observation n_channels is required; cannot resample cube")

    with fits.open(source_path) as hdul:
        hdu = hdul[0]
        if hdu.data is None:
            raise ValueError(f"{source_path} has no image data")
        data_in = np.asarray(hdu.data, dtype=np.float32)
        header_in = hdu.header.copy()
        data_in, header_in = _squeeze_degenerate_axes(data_in, header_in)
        if data_in.ndim != 3:
            raise ValueError(f"spectral_cube must be 3D, got ndim={data_in.ndim}")

    freq_axis = _find_frequency_axis(header_in)
    spatial_axes = [i for i in range(1, 4) if i != freq_axis]
    n_freq_in = int(header_in[f"NAXIS{freq_axis}"])
    n_dec_in = int(header_in[f"NAXIS{spatial_axes[1]}"])
    n_ra_in = int(header_in[f"NAXIS{spatial_axes[0]}"])

    freqs_in_hz = _freq_axis_centres(header_in, n_freq_in, axis=freq_axis)
    cdelt_freq = float(header_in[f"CDELT{freq_axis}"])
    cunit_freq = str(header_in.get(f"CUNIT{freq_axis}") or "").strip().lower()
    if cunit_freq == "mhz":
        cdelt_freq *= 1e6
    if cdelt_freq < 0:
        sort_idx = np.argsort(freqs_in_hz)
        slicer = [slice(None)] * 3
        slicer[_fits_axis_to_numpy(freq_axis, data_in.ndim)] = sort_idx
        freqs_in_hz = freqs_in_hz[sort_idx]
        data_in = data_in[tuple(slicer)]

    visibility_path = ctx.visibility_path
    ms_grid = None
    if visibility_path is not None:
        ms_grid = _read_ms_spectral_window(visibility_path)
    if ms_grid is not None:
        freqs_out_hz, obs_df_hz = ms_grid
        n_freq_out = len(freqs_out_hz)
    else:
        if obs.channel_width_mhz is None or obs.frequency_mhz is None:
            raise ValueError(
                "Cannot determine output spectral grid from MS or config"
            )
        n_freq_out = obs.n_channels
        obs_df_hz = obs.channel_width_mhz * 1e6
        obs_center_hz = obs.frequency_mhz * 1e6
        freqs_out_hz = obs_center_hz + (
            np.arange(n_freq_out) - (n_freq_out - 1) / 2.0
        ) * obs_df_hz

    if obs_df_hz == 0:
        raise ValueError("output channel width is 0; cannot build output grid")

    reverse_output = False
    if obs_df_hz < 0:
        sort_idx = np.argsort(freqs_out_hz)
        freqs_out_hz = freqs_out_hz[sort_idx]
        reverse_output = True

    cube_min_hz, cube_max_hz = float(freqs_in_hz[0]), float(freqs_in_hz[-1])
    ms_min_hz, ms_max_hz = float(freqs_out_hz[0]), float(freqs_out_hz[-1])
    logger.info(
        f"spectral cube input: {n_freq_in} channels, "
        f"[{cube_min_hz/1e6:.6f}, {cube_max_hz/1e6:.6f}] MHz, "
        f"data min={float(data_in.min()):.3e} max={float(data_in.max()):.3e} mean={float(data_in.mean()):.3e}"
    )
    logger.info(
        f"MS output grid: {n_freq_out} channels x {obs_df_hz/1e6:.6f} MHz, "
        f"[{ms_min_hz/1e6:.6f}, {ms_max_hz/1e6:.6f}] MHz"
    )
    if cube_max_hz < ms_min_hz or cube_min_hz > ms_max_hz:
        raise ValueError(
            f"spectral cube [{cube_min_hz/1e6:.3f}, {cube_max_hz/1e6:.3f}] MHz "
            f"does not overlap MS grid [{ms_min_hz/1e6:.3f}, {ms_max_hz/1e6:.3f}] MHz"
        )

    data_ordered = _reorder_cube_axes(data_in, header_in)

    logger.info(f"resampling spectral cube from {n_freq_in} to {n_freq_out} channels")
    data_out = _resample_spectral_axis_to_ms_channels(
        freqs_in_hz, data_ordered, freqs_out_hz, obs_df_hz
    )
    if reverse_output:
        data_out = data_out[::-1, :, :]
        obs_df_hz = abs(obs_df_hz)
        freqs_out_hz = freqs_out_hz[::-1]

    nonzero_out = int(np.count_nonzero(data_out))
    logger.info(
        f"resampled cube: min={float(data_out.min()):.3e} max={float(data_out.max()):.3e} "
        f"mean={float(data_out.mean()):.3e} nonzero_pixels={nonzero_out}/{data_out.size}"
    )

    phase_center = None
    if visibility_path is not None:
        phase_center = _read_ms_phase_center(visibility_path)
    if phase_center is None:
        try:
            wcs = WCS(header_in).celestial
            ra_fits_axis = spatial_axes[0]
            dec_fits_axis = spatial_axes[1]
            center_x = (n_ra_in - 1) / 2.0
            center_y = (n_dec_in - 1) / 2.0
            sky = wcs.pixel_to_world(center_x, center_y)
            if isinstance(sky, SkyCoord):
                phase_center = (float(sky.ra.deg), float(sky.dec.deg))
        except Exception as exc:
            logger.debug(
                f"prepare_spectral_cube_for_casa: WCS phase centre failed for {source_path}: {exc}"
            )
            phase_center = None
    if phase_center is None:
        try:
            ra_val = float(header_in.get(f"CRVAL{spatial_axes[0]}", 0.0))
            dec_val = float(header_in.get(f"CRVAL{spatial_axes[1]}", 0.0))
            phase_center = (ra_val, dec_val)
        except Exception as exc:
            raise ValueError(
                f"Cannot determine spatial phase centre for {source_path}"
            ) from exc

    ra_deg, dec_deg = phase_center
    cdelt1 = float(header_in.get(f"CDELT{spatial_axes[0]}", 1.0))
    cdelt2 = float(header_in.get(f"CDELT{spatial_axes[1]}", 1.0))
    cunit1 = str(header_in.get(f"CUNIT{spatial_axes[0]}") or "").strip().lower()
    cunit2 = str(header_in.get(f"CUNIT{spatial_axes[1]}") or "").strip().lower()
    if cunit1 == "rad":
        cdelt1 = np.degrees(cdelt1)
    elif cunit1 in {"asec", "arcsec"}:
        cdelt1 = cdelt1 / 3600.0
    if cunit2 == "rad":
        cdelt2 = np.degrees(cdelt2)
    elif cunit2 in {"asec", "arcsec"}:
        cdelt2 = cdelt2 / 3600.0

    header = fits.Header()
    header["NAXIS"] = 3
    header["NAXIS1"] = n_ra_in
    header["NAXIS2"] = n_dec_in
    header["NAXIS3"] = n_freq_out

    header["CTYPE1"] = "RA---SIN"
    header["CRPIX1"] = float(n_ra_in) / 2.0
    header["CRVAL1"] = float(ra_deg)
    header["CDELT1"] = float(cdelt1)
    header["CUNIT1"] = "deg"

    header["CTYPE2"] = "DEC--SIN"
    header["CRPIX2"] = float(n_dec_in) / 2.0
    header["CRVAL2"] = float(dec_deg)
    header["CDELT2"] = float(cdelt2)
    header["CUNIT2"] = "deg"

    header["CTYPE3"] = "FREQ"
    header["CRPIX3"] = 1.0
    header["CRVAL3"] = float(freqs_out_hz[0])
    header["CDELT3"] = float(obs_df_hz)
    header["CUNIT3"] = "Hz"

    bunit = str(header_in.get("BUNIT") or "Jy/px").strip().lower()
    if bunit in _ACCEPTED_JY_PER_PIXEL_UNITS:
        bunit = "Jy/pixel"
    header["BUNIT"] = bunit
    header["RESTFRQ"] = float(np.mean(freqs_out_hz))
    header["SPECSYS"] = "LSRK"
    header["SSYSOBS"] = "LSRK"
    header["HISTORY"] = "resampled to MS spectral grid by skasim.prepare_spectral_cube_for_casa"

    # wsclean -predict expects one 2D FITS per channel named <prefix>-NNNN-model.fits.
    # We keep the full resampled cube as a reference FITS but do not import to CASA.
    fits.writeto(cube_fits, data_out, header, overwrite=True)

    obs_center_hz = float(np.mean(freqs_out_hz))
    return CasaModelProduct(
        model_paths=[cube_fits],
        nterms=1,
        reffreq=f"{obs_center_hz}Hz",
        intermediates=[cube_fits],
        cube_data=data_out,
        header=header,
        freq_axis=3,
        model_dir=ctx.work_dir,
    )


def run_moment8_for_spectral_cube(
    ctx: RunContext,
    work_dir: Path,
    output_prefix: str,
    tag: str,
) -> None:
    """Generate a moment-8 (peak intensity) map and an average spectrum plot from the stacked WSClean clean cube.

    Uses pure NumPy over the FITS cube; no CASA required.
    """
    cube_path = work_dir / f"{output_prefix}-cube-image.fits"
    if not cube_path.exists():
        logger.warning(f"No cleaned cube found at {cube_path}; skipping moment-8")
        return

    with fits.open(cube_path) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header.copy()

    if data.ndim != 3:
        logger.warning(f"Cube {cube_path} has shape {data.shape}; expected 3D")
        return

    nchan = data.shape[0]
    crpix3 = header.get("CRPIX3", 1.0)
    crval3 = header.get("CRVAL3", 0.0)
    cdelt3 = header.get("CDELT3", 1.0)
    cunit3 = (header.get("CUNIT3") or "Hz").strip()
    freq_axis = crval3 + (np.arange(nchan) - (crpix3 - 1.0)) * cdelt3
    restfreq = header.get("RESTFRQ") or header.get("RESTFREQ") or header.get("RESTWAV")
    if restfreq:
        restfreq = float(restfreq)
        velocities = 299792.458 * (1.0 - freq_axis / restfreq)  # km/s
        x_label = "Velocity (km/s)"
        x_values = velocities
    else:
        velocities = None
        x_label = f"Frequency ({cunit3})"
        x_values = freq_axis

    moment8 = np.nanmax(data, axis=0)

    # average spectrum (mean over all spatial pixels)
    avg_spectrum = np.nanmean(data.reshape(nchan, -1), axis=1)
    png_spectrum = work_dir / f"{output_prefix}-avg_spectrum.png"
    try:
        import matplotlib
        matplotlib.use("Agg", force=True)
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 4))
        ax.plot(x_values, avg_spectrum * 1000.0, color="#0969da", linewidth=1.0)
        ax.set_xlabel(x_label)
        ax.set_ylabel("Mean intensity (mJy/pixel)")
        ax.set_title("Average spectrum (clean cube)")
        ax.grid(True, color="0.85", linestyle=":", linewidth=0.8)
        fig.tight_layout()
        fig.savefig(png_spectrum, dpi=150)
        plt.close(fig)
        ctx.manifest.add_output(
            "image_product",
            str(png_spectrum.relative_to(ctx.work_dir)),
            image_product_id=output_prefix,
            imager="wsclean",
            role="avg_spectrum_plot",
            metadata={"tag": tag},
        )
    except Exception as exc:
        logger.debug(f"Average spectrum plot failed: {exc}")

    base_header = header.copy()
    for key in list(base_header.keys()):
        if key in ("NAXIS", "NAXIS3") or key.startswith("NAXIS3"):
            base_header.remove(key, ignore_missing=True)
        if key in ("CRPIX3", "CRVAL3", "CDELT3", "CTYPE3", "CUNIT3"):
            base_header.remove(key, ignore_missing=True)
    base_header["NAXIS"] = 2
    base_header["NAXIS1"] = header["NAXIS1"]
    base_header["NAXIS2"] = header["NAXIS2"]

    out_fits = work_dir / f"{output_prefix}-moment8.fits"
    base_header["BUNIT"] = header.get("BUNIT") or "Jy/beam"
    base_header["MOMENT"] = 8
    base_header["HISTORY"] = "produced by skasim.run_moment8_for_spectral_cube"
    if out_fits.exists():
        out_fits.unlink()
    fits.writeto(out_fits, np.asarray(moment8, dtype=np.float32), base_header, overwrite=True)
    ctx.manifest.add_output(
        "image_product",
        str(out_fits.relative_to(ctx.work_dir)),
        image_product_id=output_prefix,
        imager="wsclean",
        role="moment8",
        metadata={"tag": tag},
    )

    png_path = out_fits.with_suffix(".png")
    try:
        from ..imaging import write_fits_preview
        write_fits_preview(out_fits, png_path, "Moment 8 (peak)")
        ctx.manifest.add_output(
            "image_product",
            str(png_path.relative_to(ctx.work_dir)),
            image_product_id=output_prefix,
            imager="wsclean",
            role="moment8_preview",
            metadata={"tag": tag},
        )
    except Exception as exc:
        logger.debug(f"Preview for moment8 failed: {exc}")


def write_spectral_cube_input_preview(
    ctx: RunContext,
    fov: u.Quantity,
) -> None:
    """Render a peak-intensity (moment-8) PNG preview of the input spectral cube.

    The preview is added to the manifest as a sky-model plot so it appears in
    the weblog's Sky Model section.
    """
    from ..config import SpectralCubeModelEntry

    cube_entries = [
        m for m in ctx.config.models if isinstance(m, SpectralCubeModelEntry)
    ]
    if not cube_entries:
        return

    for index, entry in enumerate(cube_entries, start=1):
        cube_path = Path(entry.cube).expanduser().resolve()
        if not cube_path.exists():
            logger.warning(f"spectral_cube input not found for preview: {cube_path}")
            continue
        try:
            with fits.open(cube_path) as hdul:
                hdu = hdul[0]
                data = np.asarray(hdu.data, dtype=np.float32)
                header = hdu.header.copy()
                data, header = _squeeze_degenerate_axes(data, header)
        except Exception as exc:
            logger.warning(f"Failed to read spectral cube for preview: {exc}")
            continue

        if data.ndim != 3:
            logger.warning(f"spectral_cube preview expects 3D data, got {data.shape}")
            continue

        # Determine the frequency axis dynamically; the raw FITS may have an
        # arbitrary axis ordering and ``_reorder_cube_axes`` is designed for
        # the pipeline's internal (freq, dec, ra) representation, not for raw
        # preview data.
        freq_axis = _find_frequency_axis(header)
        moment8 = np.nanmax(data, axis=_fits_axis_to_numpy(freq_axis))

        # Build a minimal 2D header for the preview FITS
        out_header = header.copy()
        for key in list(out_header.keys()):
            if key in ("NAXIS", "NAXIS3") or key.startswith("NAXIS3"):
                out_header.remove(key, ignore_missing=True)
            if key in ("CRPIX3", "CRVAL3", "CDELT3", "CTYPE3", "CUNIT3"):
                out_header.remove(key, ignore_missing=True)
        out_header["NAXIS"] = 2
        out_header["NAXIS1"] = header["NAXIS1"]
        out_header["NAXIS2"] = header["NAXIS2"]
        out_header["BUNIT"] = header.get("BUNIT") or "Jy/pixel"
        out_header["MOMENT"] = 8

        suffix = "" if len(cube_entries) == 1 else f"_{index:02d}"
        fits_name = f"{ctx.work_dir.name}_input_cube_moment8{suffix}.fits"
        png_name = f"{ctx.work_dir.name}_input_cube_moment8{suffix}.png"
        fits_path = ctx.work_dir / fits_name
        png_path = ctx.work_dir / png_name
        fits.writeto(
            fits_path,
            np.asarray(moment8, dtype=np.float32),
            out_header,
            overwrite=True,
        )

        recenter = None
        try:
            wcs = WCS(header).celestial
            pix = np.array([header["NAXIS1"] / 2.0, header["NAXIS2"] / 2.0])
            sky = wcs.pixel_to_world(*pix)
            recenter = (sky.ra.deg, sky.dec.deg, fov.to(u.deg).value)
        except Exception as exc:
            logger.debug(
                f"write_spectral_cube_input_preview: WCS recenter failed for {fits_path}: {exc}"
            )
            recenter = None

        write_fits_preview(
            fits_path,
            png_path,
            "Input spectral cube — Moment 8 (peak)",
            recenter=recenter,
            scale_factor=1000.0,
            bunit="mJy/pixel",
            colorbar_label="mJy/pixel",
        )
        ctx.manifest.add_output(
            "plot",
            png_name,
            role="input_cube_moment8",
            metadata={
                "model_entry_index": index - 1,
                "model_type": entry.type,
                "source_fits": str(cube_path),
                "preview_fits": str(fits_path),
            },
        )


def write_image_model_previews(
    ctx: RunContext,
    center: SkyCoord,
    fov: u.Quantity,
) -> None:
    """Write FITS model previews for the weblog sky-model section."""
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    for index, entry in enumerate(entries, start=1):
        image_path = primary_model_fits_path(entry)
        export_path = None
        if image_path is None and isinstance(entry, CasaTaylorTermsModelEntry):
            image_path = Path(entry.tt0).expanduser().resolve()
            export_path = ctx.work_dir / f"model_entry_{index:02d}_casa_taylor.tt0.fits"
        if image_path is None:
            continue
        # use the FITS image's own WCS center to avoid recentering NaN
        # when the model and sky-catalog coordinates differ
        try:
            info = read_fits_image_info(image_path)
            if info.center is None:
                recenter = None
            else:
                assert isinstance(info.center, SkyCoord)  # narrow for type checker
                recenter = (
                    info.center.ra.deg,
                    info.center.dec.deg,
                    fov.to(u.deg).value,
                )
        except Exception as exc:
            logger.debug(
                f"write_image_model_previews: WCS recenter failed for {image_path}: {exc}"
            )
            recenter = None
        suffix = "" if len(entries) == 1 else f"_{index:02d}"
        png_name = f"{ctx.work_dir.name}_fits_model{suffix}.png"
        png_path = ctx.work_dir / png_name
        preview_source = image_path
        if export_path is not None:
            run_casa_exportfits(ctx.work_dir, image_path, export_path)
            preview_source = export_path
        write_fits_preview(
            preview_source,
            png_path,
            "FITS Model",
            recenter=recenter,
            scale_factor=1000.0,
            bunit="mJy/pixel",
            colorbar_label="mJy/pixel",
        )
        ctx.manifest.add_output(
            "plot",
            png_name,
            role="fits_model",
            metadata={
                "model_entry_index": index - 1,
                "model_type": entry.type,
                "source_fits": str(image_path),
                "preview_fits": str(preview_source),
            },
        )

    # For spectral-cube inputs, also render a moment-8 (peak) preview of the raw cube.
    write_spectral_cube_input_preview(ctx, fov)


def inject_image_models(ctx: RunContext, visibility_path: Path) -> None:
    """Inject configured image models into an existing Measurement Set."""
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    backends = {"casa_ft", "wsclean_predict"}
    ctx.add_milestone(
        "image_injection_started",
        "started",
        details={"n_model_entries": len(entries), "backends": sorted(backends)},
    )

    for index, entry in enumerate(entries):
        if isinstance(entry, StaticStokesMapsModelEntry):
            raise NotImplementedError(
                "static_stokes_maps is schema-ready, but the CASA backend path is "
                "planned for the next implementation phase."
            )
        if isinstance(entry, ContinuumIAlphaModelEntry):
            report = validate_continuum_i_alpha(entry)
            product = prepare_continuum_i_alpha_for_casa(ctx, entry, index)

            if entry.injection_backend == "casa_ft":
                logger.warning(
                    "continuum_i_alpha injection_backend='casa_ft' is deprecated; "
                    "prefer 'wsclean_predict'."
                )
                run_casa_ft(
                    visibility_path=visibility_path,
                    model_paths=product.model_paths,
                    nterms=product.nterms,
                    reffreq=product.reffreq,
                    incremental=index > 0,
                )
                backend = "casa_ft"
            else:
                img_config = next(
                    (img for img in ctx.config.imaging if img.imager == "wsclean"),
                    ctx.config.imaging[0],
                )
                from .wsclean_predict import inject_continuum_i_alpha_with_wsclean_predict

                report_predict = inject_continuum_i_alpha_with_wsclean_predict(
                    ctx,
                    entry,
                    index,
                    visibility_path,
                    img_config,
                    product,
                )
                backend = "wsclean_predict"
                report = {**report, **report_predict}
        elif isinstance(entry, CasaTaylorTermsModelEntry):
            report = validate_casa_taylor_terms(entry)
            product = prepare_casa_taylor_terms(ctx, entry, index)
            logger.warning(
                "casa_taylor_terms uses the deprecated CASA ft backend; "
                "consider migrating to continuum_i_alpha with wsclean_predict."
            )
            run_casa_ft(
                visibility_path=visibility_path,
                model_paths=product.model_paths,
                nterms=product.nterms,
                reffreq=product.reffreq,
                incremental=index > 0,
            )
            backend = "casa_ft"
        elif isinstance(entry, SpectralCubeModelEntry):
            img_config = next(
                (img for img in ctx.config.imaging if img.imager == "wsclean"),
                ctx.config.imaging[0],
            )
            report = validate_spectral_cube(entry, ctx.config.observation, img_config)
            product = prepare_spectral_cube_for_casa(ctx, entry, index, report)
            from .wsclean_predict import inject_spectral_cube_with_wsclean_predict

            assert product.cube_data is not None
            assert product.header is not None
            assert product.freq_axis is not None
            report_predict = inject_spectral_cube_with_wsclean_predict(
                ctx,
                entry,
                index,
                visibility_path,
                img_config,
                product.cube_data,
                product.header,
                product.freq_axis,
            )
            backend = "wsclean_predict"
            report = {**report, **report_predict}
        else:
            continue
        ctx.manifest.add_output(
            "sky_model",
            product.model_paths[0].name,
            role="casa_model_image",
            metadata={
                "model_entry_index": index,
                "model_type": entry.type,
                "nterms": product.nterms,
                "reffreq": product.reffreq,
                "all_model_paths": [path.name for path in product.model_paths],
            },
        )
        ctx.add_milestone(
            "image_model_injected",
            "completed",
            details={
                "model_entry_index": index,
                "model_type": entry.type,
                "backend": backend,
                **report,
            },
        )

    merge_model_data_into_data(visibility_path)
    ctx.add_milestone(
        "image_injection_completed",
        "completed",
        details={"visibility_path": str(visibility_path), "model_data_merged": True},
    )


def validate_casa_taylor_terms(entry: CasaTaylorTermsModelEntry) -> dict:
    """Validate an existing CASA Taylor-term image model entry."""
    model_paths = [
        Path(path).expanduser().resolve()
        for path in (entry.tt0, entry.tt1)
        if path is not None
    ]
    if not model_paths:
        raise ValueError("casa_taylor_terms requires at least tt0.")
    for path in model_paths:
        if not path.is_dir():
            raise ValueError(f"{path} must be a CASA image table directory.")
        if not (path / "table.dat").exists():
            raise ValueError(f"{path} does not look like a CASA image table.")
    return {
        "model_paths": [str(path) for path in model_paths],
        "nterms": len(model_paths),
        "reference_frequency_hz": entry.reference_frequency_hz,
    }


def prepare_casa_taylor_terms(
    ctx: RunContext,
    entry: CasaTaylorTermsModelEntry,
    index: int,
) -> CasaModelProduct:
    """Copy CASA Taylor-term images into the run and align their spectral reference.

    The reference frequency is adjusted to the observation band centre.
    For nterms≥2, tt0 pixel data is scaled:  tt0' = tt0 · (ν_obs / ν_old)^α
    where α = mean(tt1) / mean(tt0).  tt1 pixel data is unchanged.
    For nterms=1, only CRVAL4 is updated.
    """
    new_ref_hz = ctx.config.observation.frequency_mhz * 1e6
    old_ref_hz = entry.reference_frequency_hz

    source_paths = [
        Path(path).expanduser().resolve()
        for path in (entry.tt0, entry.tt1)
        if path is not None
    ]
    prefix = f"model_entry_{index + 1:02d}_casa_taylor"
    model_paths = []
    for term_index, source_path in enumerate(source_paths):
        target_path = ctx.work_dir / f"{prefix}.tt{term_index}.image"
        if target_path.exists():
            shutil.rmtree(target_path)
        shutil.copytree(source_path, target_path)
        model_paths.append(target_path)

    nterms = len(model_paths)

    if nterms >= 2:
        # alpha map: element-wise tt1/tt0, with safeguard for tt0==0 (those pixels
        # stay zero regardless of spectral index).
        casacore_table = require_casacore()
        with casacore_table(str(model_paths[0]), readonly=True, ack=False) as tbl0:
            tt0_data = np.asarray(tbl0.getcol("map"))
        with casacore_table(str(model_paths[1]), readonly=True, ack=False) as tbl1:
            tt1_data = np.asarray(tbl1.getcol("map"))

        # compute element-wise alpha, safeguarding divide-by-zero
        with np.errstate(divide="ignore", invalid="ignore"):
            alpha_map = np.where(
                tt0_data != 0,
                tt1_data / tt0_data,
                0.0,
            )
        alpha_mean = float(np.mean(alpha_map))
        logger.info(
            f"prepare_casa_taylor_terms: alpha_map_mean={alpha_mean:.6f} "
            f"from element-wise tt1/tt0"
        )

        adjust_spectral_reference(
            model_paths[0],
            old_ref_hz,
            new_ref_hz,
            alpha_map=alpha_map,
        )
        # tt1 must be scaled by the same factor as tt0 so that the ratio
        # tt1/tt0 (the per-pixel spectral index) is preserved when the
        # reference frequency changes.
        adjust_spectral_reference(
            model_paths[1],
            old_ref_hz,
            new_ref_hz,
            alpha_map=alpha_map,
        )
    else:
        # nterms=1: spectrally flat, only set CRVAL4
        adjust_spectral_reference(
            model_paths[0],
            old_ref_hz,
            new_ref_hz,
            alpha_map=None,
        )

    ctx.add_milestone(
        "adjusted_spectral_reference",
        "completed",
        details={
            "model_type": "casa_taylor_terms",
            "old_reference_frequency_hz": old_ref_hz,
            "new_reference_frequency_hz": new_ref_hz,
            "nterms": nterms,
        },
    )

    return CasaModelProduct(
        model_paths=model_paths,
        nterms=nterms,
        reffreq=f"{new_ref_hz}Hz",
        intermediates=model_paths,
    )


def prepare_continuum_i_alpha_for_casa(
    ctx: RunContext,
    entry: ContinuumIAlphaModelEntry,
    index: int,
) -> CasaModelProduct:
    """Create CASA image products for a continuum I+alpha model.

    Adjusts the spectral reference to the observation band centre using
    the explicit spectral index from the model entry.
    """
    new_ref_hz = ctx.config.observation.frequency_mhz * 1e6
    old_ref_hz = entry.reference_frequency_hz

    stokes_path = Path(entry.stokes_i).expanduser().resolve()
    alpha_path = Path(entry.alpha).expanduser().resolve()
    prefix = f"model_entry_{index + 1:02d}_continuum"
    tt0_fits = ctx.work_dir / f"{prefix}.tt0.fits"
    tt1_fits = ctx.work_dir / f"{prefix}.tt1.fits"
    tt0_image = ctx.work_dir / f"{prefix}.tt0.image"
    tt1_image = ctx.work_dir / f"{prefix}.tt1.image"

    with fits.open(stokes_path) as stokes_hdul, fits.open(alpha_path) as alpha_hdul:
        stokes_data = np.asarray(stokes_hdul[0].data, dtype=float)
        alpha_data = np.asarray(alpha_hdul[0].data, dtype=float)
        header = stokes_hdul[0].header.copy()
        header["BUNIT"] = stokes_hdul[0].header.get("BUNIT", "Jy/pixel")

        # ensure both model FITS have a degenerate 4D shape (STOKES, FREQ)
        # so CASA importfits creates images with a frequency axis that ft can map.
        stokes_4d, header_4d = _as_4d_image(stokes_data, header, new_ref_hz)
        alpha_4d = _broadcast_to_4d(alpha_data)

        # Scale data to the new reference frequency
        factor = (new_ref_hz / old_ref_hz) ** alpha_4d
        stokes_4d_scaled = stokes_4d * factor
        # CASA Taylor terms: tt0 is I(nu_ref), tt1 is I(nu_ref) * alpha
        tt1_data = stokes_4d_scaled * alpha_4d

        fits.writeto(tt0_fits, stokes_4d_scaled, header=header_4d, overwrite=True)
        fits.writeto(tt1_fits, tt1_data, header=header_4d, overwrite=True)

    for imagename in (tt0_image, tt1_image):
        if imagename.exists():
            shutil.rmtree(imagename)
    casa_tasks = import_casa_tasks()
    if casa_tasks is not None:
        importfits, _ = casa_tasks
        importfits(fitsimage=str(tt0_fits), imagename=str(tt0_image), overwrite=True)
        importfits(fitsimage=str(tt1_fits), imagename=str(tt1_image), overwrite=True)
    else:
        run_casa_importfits(
            ctx.work_dir,
            [(tt0_fits, tt0_image), (tt1_fits, tt1_image)],
        )

    # Pixel data was already scaled; we just use adjust_spectral_reference
    # with alpha_map=None to log the new frequency and enforce CRVAL4.
    alpha_mean = float(np.mean(alpha_data))
    adjust_spectral_reference(
        tt0_image,
        old_ref_hz,
        new_ref_hz,
        alpha_map=None,
    )
    adjust_spectral_reference(
        tt1_image,
        old_ref_hz,
        new_ref_hz,
        alpha_map=None,
    )

    ctx.add_milestone(
        "adjusted_spectral_reference",
        "completed",
        details={
            "model_type": "continuum_i_alpha",
            "old_reference_frequency_hz": old_ref_hz,
            "new_reference_frequency_hz": new_ref_hz,
            "alpha_mean": alpha_mean,
            "nterms": 2,
        },
    )

    return CasaModelProduct(
        model_paths=[tt0_image, tt1_image],
        nterms=2,
        reffreq=f"{new_ref_hz}Hz",
        intermediates=[tt0_fits, tt1_fits],
    )


def _as_4d_image(data: np.ndarray, header: fits.Header, reference_freq_hz: float) -> tuple[np.ndarray, fits.Header]:
    """Return a 4D view (FREQ, STOKES, Y, X) of a 2D FITS image plus an updated header."""
    # collapse any leading degenerate axes and validate spatial dimensions
    flat = np.asarray(data).squeeze()
    if flat.ndim != 2:
        raise ValueError(f"continuum_i_alpha expects a 2D spatial image; got shape {data.shape}")
    ny, nx = flat.shape
    image_4d = flat.reshape(1, 1, ny, nx)

    new_header = header.copy()
    new_header["NAXIS"] = 4
    new_header["NAXIS1"] = nx
    new_header["NAXIS2"] = ny
    new_header["NAXIS3"] = 1
    new_header["NAXIS4"] = 1
    new_header["CTYPE1"] = header.get("CTYPE1", "RA---SIN")
    new_header["CTYPE2"] = header.get("CTYPE2", "DEC--SIN")
    new_header["CTYPE3"] = "STOKES"
    new_header["CTYPE4"] = "FREQ"
    new_header["CRPIX1"] = header.get("CRPIX1", nx / 2.0)
    new_header["CRPIX2"] = header.get("CRPIX2", ny / 2.0)
    new_header["CRPIX3"] = 1.0
    new_header["CRPIX4"] = 1.0
    new_header["CRVAL3"] = 1.0
    new_header["CRVAL4"] = reference_freq_hz
    new_header["CDELT3"] = 1.0
    new_header["CDELT4"] = 1.0
    new_header["CUNIT3"] = ""
    new_header["CUNIT4"] = "Hz"
    return image_4d, new_header


def _broadcast_to_4d(data: np.ndarray) -> np.ndarray:
    """Collapse leading degenerate axes and broadcast a 2D spatial map to 4D."""
    flat = np.asarray(data).squeeze()
    if flat.ndim != 2:
        raise ValueError(f"continuum_i_alpha alpha map must be 2D spatial; got shape {data.shape}")
    return flat.reshape(1, 1, *flat.shape)

def adjust_spectral_reference(
    image_path: Path,
    old_ref_hz: float,
    new_ref_hz: float,
    alpha_map: np.ndarray | None = None,
) -> float:
    """Adjust the spectral reference of a CASA image to the observation band centre.

    For nterms=1 (alpha_map is None), set CRVAL4 to new_ref_hz only — the model
    is spectrally flat and no pixel-data correction is needed.

    For nterms≥2 (alpha_map provided), correct the pixel data element-wise following
    CASA's Taylor-series convention::

        tt0'(x,y) = tt0(x,y) · (ν_new / ν_old) ^ α(x,y)

    where α(x,y) = tt1(x,y) / tt0(x,y)   for each pixel.  CRVAL4 is also set.

    Returns the adjusted reference frequency in Hz (always new_ref_hz).
    """
    if alpha_map is not None:
        casacore_table = require_casacore()

        # element-wise scaling
        factor = (new_ref_hz / old_ref_hz) ** alpha_map
        with casacore_table(str(image_path), readonly=False, ack=False) as tbl:
            data = np.asarray(tbl.getcol("map"))
            corrected = data * factor
            tbl.putcol("map", corrected)
        logger.info(
            f"Adjusted pixel data in {image_path}: ref_freq {old_ref_hz:.3e}Hz → "
            f"{new_ref_hz:.3e}Hz (factor range: [{np.min(factor):.4f}, {np.max(factor):.4f}])"
        )
    else:
        logger.info(
            f"Spectrally flat image {image_path}: updating CRVAL4 = {new_ref_hz:.3e}Hz "
            "(no pixel-data correction)"
        )

    _set_crval4_via_script(image_path.parent, [image_path], new_ref_hz)
    return new_ref_hz


def _image_has_spectral_axis(image_path: Path) -> bool:
    """Return True if the CASA image has a frequency/spectral axis."""
    try:
        from casacore.tables import table
    except Exception as exc:
        raise RuntimeError(
            "python-casacore is required to inspect CASA image tables."
        ) from exc

    with table(str(image_path), readonly=True, ack=False) as tbl:
        # check primary coordinate system dictionary for a 'spectralN' key
        if "coords" in tbl.keywordnames():
            coords = tbl.getkeyword("coords")
            if any(key.startswith("spectral") for key in coords.keys()):
                return True

        # fallback to map column dimensional names
        if "map" in tbl.colnames():
            map_kws = tbl.getcolkeywords("map")
            dim_names = [str(n).lower() for n in map_kws.get("dimnames", [])]
            return any("freq" in name for name in dim_names)

    return False



def _set_crval4_via_script(
    work_dir: Path,
    image_paths: list[Path],
    frequency_hz: float,
) -> None:
    """Set CRVAL4 on CASA images — subprocess fallback for environments without casatasks."""
    # 2D images have no spectral axis; CASA ft will treat them as spectrally flat
    # with the reference frequency passed as reffreq, so CRVAL4 cannot (and need not) be set.
    image_paths = [p for p in image_paths if _image_has_spectral_axis(p)]
    if not image_paths:
        logger.warning(
            "Spectrally flat 2D image(s): skipping CRVAL4 update; ft reffreq carries the reference frequency"
        )
        return

    try:
        from casatasks import imhead

        for image_path in image_paths:
            imhead(
                imagename=str(image_path),
                mode="put",
                hdkey="crval4",
                hdvalue=f"{frequency_hz}Hz",
            )
        return
    except Exception as exc:
        logger.debug(f"_set_crval4_via_script: in-process imhead failed: {exc}")
        # fall through to batch mode

    run_casa_set_spectral_coordinate(work_dir, image_paths, frequency_hz)


def import_casa_tasks():
    """Return in-process CASA tasks when they are importable in this Python env."""
    try:
        from casatasks import ft, importfits
    except Exception as exc:
        logger.debug(f"import_casa_tasks: casatasks not available: {exc}")
        return None
    return importfits, ft


def import_casa_tools():
    """Return in-process CASA simulator tool when importable."""
    try:
        from casatools import simulator
    except Exception as exc:
        logger.debug(f"import_casa_tools: casatools not available: {exc}")
        return None
    return simulator


def require_casa_tasks():
    """Import CASA tasks lazily and provide a clear runtime error if unavailable."""
    casa_tasks = import_casa_tasks()
    if casa_tasks is None:
        raise RuntimeError(
            "CASA casatasks.importfits and casatasks.ft are required for "
            "in-process image-model injection. Install casatasks in this "
            "environment or make the CASA executable available on PATH for "
            "batch-mode fallback."
        )
    return casa_tasks


def require_casa_tools():
    """Import CASA simulator tool lazily and provide a clear runtime error if unavailable."""
    simulator_cls = import_casa_tools()
    if simulator_cls is None:
        raise RuntimeError(
            "CASA casatools.simulator is required for in-process "
            "sm.predict image-model injection. Install casatools in this "
            "environment or make the CASA executable available on PATH for "
            "batch-mode fallback."
        )
    return simulator_cls


def require_casa_executable() -> Path:
    """Return a CASA executable for batch-mode fallback."""
    executable = shutil.which("casa")
    if executable is None:
        raise RuntimeError(
            "CASA image-model injection requires either importable casatasks "
            "or a casa executable on PATH."
        )
    return Path(executable)


def run_casa_importfits(
    work_dir: Path,
    images: list[tuple[Path, Path]],
) -> None:
    """Run CASA importfits in batch mode for prepared FITS images."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_importfits.py"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import importfits",
        "except Exception as _e:",
        "    raise RuntimeError('casatasks.importfits is not available: ' + str(_e))",
    ]
    for fitsimage, imagename in images:
        lines.append(
            "importfits(fitsimage={!r}, imagename={!r}, overwrite=True)".format(
                str(fitsimage),
                str(imagename),
            )
        )
    run_casa_script(executable, script_path, lines)
    # verify side-effects; if not present, casatasks ran but failed silently
    for _fitsimage, imagename in images:
        table_dat = imagename / "table.dat"
        if not table_dat.exists():
            raise RuntimeError(
                f"CASA importfits did not create {imagename} as expected."
            )


def run_casa_exportfits(
    work_dir: Path,
    imagename: Path,
    fitsimage: Path,
) -> None:
    """Run CASA exportfits in batch mode for a CASA image table."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_exportfits.py"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import exportfits",
        "except Exception:",
        "    pass",
        "exportfits(imagename={!r}, fitsimage={!r}, overwrite=True)".format(
            str(imagename),
            str(fitsimage),
        ),
    ]
    run_casa_script(executable, script_path, lines)


def run_casa_set_spectral_coordinate(
    work_dir: Path,
    image_paths: list[Path],
    frequency_hz: float,
) -> None:
    """Set the single-channel spectral coordinate of CASA images to the run reference."""
    executable = require_casa_executable()
    script_path = work_dir / "skasim_casa_set_spectral_coordinate.py"
    hdvalue = f"{frequency_hz}Hz"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import imhead",
        "except Exception:",
        "    pass",
    ]
    for image_path in image_paths:
        lines.append(
            "imhead(imagename={!r}, mode='put', hdkey='crval4', hdvalue={!r})".format(
                str(image_path),
                hdvalue,
            )
        )
    run_casa_script(executable, script_path, lines)


def run_casa_ft(
    visibility_path: Path,
    model_paths: list[Path],
    nterms: int,
    reffreq: str,
    incremental: bool,
) -> None:
    """Run CASA ft into MODEL_DATA for one prepared model entry."""
    logger.info(
        f"CASA ft model={[str(path) for path in model_paths]} "
        f"nterms={nterms} reffreq={reffreq} incremental={incremental}"
    )
    casa_tasks = import_casa_tasks()
    if casa_tasks is not None:
        _, ft = casa_tasks
        ft(
            vis=str(visibility_path),
            model=[str(path) for path in model_paths],
            nterms=nterms,
            reffreq=reffreq,
            incremental=incremental,
            usescratch=True,
        )
        return

    executable = require_casa_executable()
    script_path = visibility_path.parent / "skasim_casa_ft.py"
    model_literal = "[" + ", ".join(repr(str(path)) for path in model_paths) + "]"
    lines = [
        "# Auto-generated by skasim; safe to delete after the run.",
        "try:",
        "    from casatasks import ft",
        "except Exception as _e:",
        "    raise RuntimeError('casatasks.ft is not available: ' + str(_e))",
        "ft(",
        f"    vis={str(visibility_path)!r},",
        f"    model={model_literal},",
        f"    nterms={int(nterms)},",
        f"    reffreq={reffreq!r},",
        f"    incremental={bool(incremental)!r},",
        "    usescratch=True,",
        ")",
    ]
    run_casa_script(executable, script_path, lines)


def run_casa_predict(
    visibility_path: Path,
    model_paths: list[Path],
    incremental: bool,
) -> None:
    """Run CASA simulator.predict into MODEL_DATA for a spectral-cube model.

    .. deprecated::
        This CASA path has been replaced by ``wsclean -predict`` for spectral
        cubes.  It is kept as a visible tombstone so existing callers get a clear
        error instead of a missing-symbol failure.
    """
    raise RuntimeError(
        "CASA sm.predict path is deprecated; use wsclean -predict instead."
    )


def run_casa_script(executable: Path, script_path: Path, lines: list[str]) -> None:
    """Write and execute one CASA batch script, surfacing useful failure output."""
    script_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    log_path = script_path.with_suffix(".log")
    command = [
        str(executable),
        "--nologger",
        "--nogui",
        "--log2term",
        "-c",
        str(script_path),
    ]
    logger.info(f"CASA batch command: {' '.join(command)}")
    with log_path.open("w", encoding="utf-8") as log_file:
        result = subprocess.run(
            command,
            cwd=str(script_path.parent),
            text=True,
            stdin=subprocess.DEVNULL,
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
        if result.returncode != 0:
            tail = "\n".join(log_path.read_text(encoding="utf-8").splitlines()[-80:])
            raise RuntimeError(
                f"CASA batch command failed with exit code {result.returncode}: "
                f"{script_path}\n{tail}"
            )


def merge_model_data_into_data(visibility_path: Path) -> None:
    """Add image-model MODEL_DATA into the delivered DATA column."""
    try:
        from casacore.tables import table
    except Exception as exc:
        raise RuntimeError(
            "python-casacore is required to merge MODEL_DATA into DATA."
        ) from exc

    with table(str(visibility_path), readonly=False, ack=False) as ms_table:
        columns = set(ms_table.colnames())
        if "DATA" not in columns or "MODEL_DATA" not in columns:
            raise ValueError(
                f"{visibility_path} must contain DATA and MODEL_DATA columns "
                "after CASA ft injection."
            )
        data = ms_table.getcol("DATA")
        model_data = ms_table.getcol("MODEL_DATA")
        ms_table.putcol("DATA", data + model_data)
