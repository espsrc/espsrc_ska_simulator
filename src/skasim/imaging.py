"""imaging.py — dirty (OSKAR) and cleaned (WSClean) imaging wrappers."""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from loguru import logger

from .config import ImgConfig, has_spectral_cube_model
from .manifest import RunContext
from .runtime import require_karabo_module
from .utils import mapping_unit

SKY_MODEL_CMAP = "viridis_r"

# --------------------------------------------------------------------------- #
# spectral-cube imaging helpers
# --------------------------------------------------------------------------- #


def _resolve_spectral_cube_wsclean_config(
    img_config: ImgConfig,
    n_channels: int,
) -> ImgConfig:
    """Return an imaging config adjusted for per-channel cube output.

    - join_channels forced to False
    - channels_out forced to n_channels if the user did not explicitly set it
    - multiscale disabled (line cubes rarely benefit from multiscale)
    """
    kwargs = img_config.model_dump()
    kwargs["join_channels"] = False
    kwargs["multiscale"] = False
    kwargs["multiscale_scales"] = None
    if img_config.channels_out is None:
        kwargs["channels_out"] = n_channels
    return ImgConfig(**kwargs)


# --------------------------------------------------------------------------- #
# dirty imaging (OSKAR)
# --------------------------------------------------------------------------- #


def run_dirty_imaging(
    ctx: RunContext,
    visibility_path: Path,
    fov: u.Quantity,
    center: SkyCoord,
    img_config: ImgConfig,  # << NEW
    sub_dir: Path,  # << NEW — work_dir/{tag}
) -> None:
    """produce dirty image via OSKAR."""
    imager_module = require_karabo_module("karabo.imaging.imager_oskar")
    visibility_module = require_karabo_module("karabo.simulation.visibility")
    work_dir = sub_dir
    vis = visibility_module.Visibility(str(visibility_path))
    imaging_cellsize = fov / img_config.pixels
    cfg = imager_module.OskarDirtyImagerConfig(
        imaging_npixel=img_config.pixels,
        imaging_cellsize=imaging_cellsize.to(u.rad).value,
        combine_across_frequencies=True,
        imaging_phase_centre=center.icrs,
    )
    imager = imager_module.OskarDirtyImager(config=cfg)
    dirty_image = imager.create_dirty_image(vis)

    dirty_png = work_dir / f"{img_config.tag}_dirty.png"
    dirty_fits = work_dir / f"{img_config.tag}_dirty.fits"
    dirty_image.write_to_file(str(dirty_fits), overwrite=True)
    try:
        write_fits_preview(dirty_fits, dirty_png, "OSKAR Dirty Image")
    except Exception as e:
        logger.warning(f"Failed to generate APLpy dirty image preview: {e}")

    logger.debug(f"Dirty PNG: {dirty_png}")
    logger.debug(f"Dirty FITS: {dirty_fits}")

    image_product_id = f"{img_config.tag}_dirty"
    ctx.manifest.add_output(
        "image_product",
        str(dirty_png.relative_to(ctx.work_dir)),
        image_product_id=image_product_id,
        imager="oskar-dirty",
        role="preview",
        metadata={"tag": img_config.tag},
    )
    ctx.manifest.add_output(
        "image_product",
        str(dirty_fits.relative_to(ctx.work_dir)),
        image_product_id=image_product_id,
        imager="oskar-dirty",
        role="dirty",
        metadata={"tag": img_config.tag},
    )


# --------------------------------------------------------------------------- #
# cleaned imaging (WSClean)
# --------------------------------------------------------------------------- #


def build_wsclean_argv(
    img_config: ImgConfig,
    visibility_path: Path,
    fov: u.Quantity,
    output_prefix: str,
    n_channels: int = 1,
) -> list[str]:
    """Build a shell-free WSClean argv list from the resolved imaging config."""
    imaging_cellsize = fov / img_config.pixels
    channels_out = (
        img_config.channels_out if img_config.channels_out is not None else n_channels
    )

    argv = shlex.split(img_config.wsclean_command) + [
        "-weight",
        "briggs",
        str(img_config.robust),
        "-size",
        str(img_config.pixels),
        str(img_config.pixels),
        "-scale",
        f"{imaging_cellsize.to(u.arcsec).value:.6f}asec",
        "-niter",
        str(img_config.clean_iterations),
        "-mgain",
        str(img_config.mgain if img_config.mgain is not None else 0.8),
        "-auto-threshold",
        str(
            img_config.auto_threshold if img_config.auto_threshold is not None else 0.3
        ),
        "-auto-mask",
        str(img_config.auto_mask if img_config.auto_mask is not None else 3.0),
        "-channels-out",
        str(channels_out),
        "-name",
        output_prefix,
    ]

    if img_config.multiscale is not False:
        argv.append("-multiscale")
        if img_config.multiscale_scales:
            argv += [
                "-multiscale-scales",
                ",".join(str(s) for s in img_config.multiscale_scales),
            ]
    if img_config.local_rms is not False:
        argv.append("-local-rms")
    if img_config.join_channels is not False:
        argv.append("-join-channels")
    if img_config.padding is not None:
        argv += ["-padding", str(img_config.padding)]
    if img_config.threads is not None:
        argv += ["-j", str(img_config.threads)]

    argv.append(str(visibility_path))
    return argv


def run_wsclean_command(argv: list[str], work_dir: Path):
    """Run WSClean with argv and an explicit working directory.

    Streams stdout/stderr line-by-line through loguru so that WSClean
    progress appears in the skasim logs in real time.
    """
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    lines: list[str] = []
    with subprocess.Popen(
        argv,
        shell=False,
        cwd=str(work_dir),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    ) as proc:
        assert proc.stdout is not None  # guaranteed by PIPE
        for line in proc.stdout:
            stripped = line.rstrip("\n")
            logger.info("[wsclean] {}", stripped)
            lines.append(line)
        proc.wait()
    combined = "".join(lines)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode,
            argv,
            output=combined,
            stderr=None,
        )
    return subprocess.CompletedProcess(
        argv, proc.returncode, stdout=combined, stderr=""
    )


def wsclean_output_prefix(ctx: RunContext) -> str:
    """Return the stable WSClean output prefix for this run."""
    return f"{ctx.work_dir.name}_wsclean"


def collect_wsclean_outputs(work_dir: Path, output_prefix: str) -> list[Path]:
    """Collect WSClean FITS outputs for one configured output prefix."""
    return sorted(work_dir.glob(f"{output_prefix}*.fits"))


def run_wsclean_imaging(
    ctx: RunContext,
    visibility_path: Path,
    fov: u.Quantity,
    img_config: ImgConfig,  # << NEW
    sub_dir: Path,  # << NEW — work_dir/{tag}
    n_channels: int = 1,
) -> None:
    """produce cleaned image via external WSClean binary."""
    wsclean_module = require_karabo_module("karabo.imaging.imager_wsclean")
    file_handler_module = require_karabo_module("karabo.util.file_handler")
    work_dir = sub_dir

    output_prefix = f"{img_config.tag}_wsclean"

    # spectral-cube mode: force per-channel imaging and no joined-channel fit
    spectral_cube_present = has_spectral_cube_model(ctx.config)
    if spectral_cube_present:
        img_config = _resolve_spectral_cube_wsclean_config(
            img_config,
            n_channels=n_channels,
        )

    argv = build_wsclean_argv(
        img_config,
        visibility_path,
        fov,
        output_prefix=output_prefix,
        n_channels=n_channels,
    )
    logger.info(f"WSClean command: {argv}")

    file_handler_module.FileHandler().get_tmp_dir(
        prefix=wsclean_module.TMP_PREFIX_CUSTOM,
        purpose=wsclean_module.TMP_PURPOSE_CUSTOM,
    )
    run_wsclean_command(argv, work_dir)

    # remove the temporary files created by WSClean
    for tmp in work_dir.glob("wsclean-00*.fits"):
        if tmp.name.startswith(output_prefix):
            continue
        try:
            tmp.unlink()
        except Exception as exc:
            logger.exception("Failed to clean up temp file %s", tmp)

    wsclean_outputs = collect_wsclean_outputs(work_dir, output_prefix)
    mfs_files = [p.name for p in wsclean_outputs if "-MFS-" in p.name]
    logger.info(f"MFS files: {mfs_files}")

    for img_path in wsclean_outputs:
        png_name = img_path.with_suffix(".png").name
        png_path = work_dir / png_name

        # infer image type from filename for correct plot title
        title = "Imaging output (WSClean)"
        if "MFS-image" in img_path.name:
            title = "Cleaned image (WSClean)"
        elif "MFS-model" in img_path.name:
            title = "Component model (WSClean)"
        elif "MFS-residual" in img_path.name:
            title = "Residual (WSClean)"
        elif "MFS-dirty" in img_path.name:
            title = "Dirty image (WSClean)"
        elif "MFS-psf" in img_path.name:
            title = "Point spread function (WSClean)"

        write_fits_preview(img_path, png_path, title)
        role = "image"
        lower_name = img_path.name.lower()
        if "model" in lower_name:
            role = "model"
        elif "residual" in lower_name:
            role = "residual"
        elif "dirty" in lower_name:
            role = "dirty"
        elif "psf" in lower_name:
            role = "psf"
        ctx.manifest.add_output(
            "image_product",
            str(png_path.relative_to(ctx.work_dir)),
            image_product_id=output_prefix,
            imager="wsclean",
            role=f"{role}_preview",
            metadata={"tag": img_config.tag},
        )
        ctx.manifest.add_output(
            "image_product",
            str(img_path.relative_to(ctx.work_dir)),
            image_product_id=output_prefix,
            imager="wsclean",
            role=role,
            metadata={"tag": img_config.tag},
        )

    # collect per-channel FITS clean images and stack them into a single 3D cube.
    # WSClean writes per-channel products as <prefix>-<dddd>-<product>.fits.
    per_channel_images = sorted(work_dir.glob(f"{output_prefix}-????-image.fits"))
    if per_channel_images:
        cube_fits = work_dir / f"{output_prefix}-cube-image.fits"
        if cube_fits.exists():
            cube_fits.unlink()
        try:
            stack_channels(per_channel_images, cube_fits)
        except Exception as exc:
            logger.warning(f"Failed to stack clean cube: {exc}")
        else:
            png_name = cube_fits.with_suffix(".png").name
            png_path = work_dir / png_name
            write_fits_preview(cube_fits, png_path, "Spectral cube (clean)")
            ctx.manifest.add_output(
                "image_product",
                str(png_path.relative_to(ctx.work_dir)),
                image_product_id=output_prefix,
                imager="wsclean",
                role="cube_image_preview",
                metadata={"tag": img_config.tag},
            )
            ctx.manifest.add_output(
                "image_product",
                str(cube_fits.relative_to(ctx.work_dir)),
                image_product_id=output_prefix,
                imager="wsclean",
                role="cube_image",
                metadata={"tag": img_config.tag},
            )

    # produce only moment 8 for spectral-cube clean outputs
    if spectral_cube_present and img_config.imager == "wsclean":
        from .loaders.image_models import run_moment8_for_spectral_cube

        run_moment8_for_spectral_cube(ctx, work_dir, output_prefix, img_config.tag)


def stack_channels(channel_paths: list[Path], output_path: Path) -> None:
    """Stack WSClean per-channel FITS images into a single 3D spectral cube.

    Uses the same logic as ``scripts/wsclean_channels_to_cube.py`` but inlined
    here so it works without putting ``scripts/`` on ``PYTHONPATH``.
    """
    import re

    channels: list[tuple[float, int, np.ndarray, fits.Header]] = []
    for path in channel_paths:
        match = re.search(r"-([0-9]{4,})-[a-zA-Z0-9_]+\.fits$", path.name)
        if not match:
            continue
        ch = int(match.group(1))
        with fits.open(path) as hdul:
            hdu = hdul[0]
            data = np.squeeze(np.asarray(hdu.data, dtype=np.float32))
            if data.ndim != 2:
                raise ValueError(
                    f"{path}: expected 2D image after squeezing, got shape {hdu.data.shape}"
                )
            header = hdu.header.copy()
        crpix3 = header.get("CRPIX3", 1.0)
        crval3 = header.get("CRVAL3")
        cdelt3 = header.get("CDELT3")
        cunit3 = (header.get("CUNIT3") or "").strip().lower()
        if crval3 is None or cdelt3 is None:
            raise ValueError(f"{path}: missing CRVAL3/CDELT3")
        freq_hz = crval3 + (crpix3 - 1.0) * cdelt3
        if cunit3 == "mhz":
            freq_hz *= 1e6
        channels.append((float(freq_hz), ch, data, header))

    if not channels:
        raise ValueError("no valid WSClean per-channel images found")

    channels.sort(key=lambda x: (x[0], x[1]))
    first_data = channels[0][2]
    cube = np.empty((len(channels), *first_data.shape), dtype=np.float32)
    for i, (_, ch, data, _) in enumerate(channels):
        if data.shape != first_data.shape:
            raise ValueError(
                f"channel {ch} has shape {data.shape}, expected {first_data.shape}"
            )
        cube[i] = data

    out_header = channels[0][3].copy()
    # force 3D header
    for key in list(out_header.keys()):
        if re.match(r"NAXIS\d+", key):
            out_header.remove(key, ignore_missing=True)
    out_header["NAXIS"] = 3
    out_header["NAXIS1"] = channels[0][3].get("NAXIS1")
    out_header["NAXIS2"] = channels[0][3].get("NAXIS2")
    out_header["NAXIS3"] = len(channels)
    freqs = np.array([c[0] for c in channels])
    df = freqs[1] - freqs[0] if len(freqs) > 1 else 1.0
    out_header["CRPIX3"] = 1.0
    out_header["CRVAL3"] = float(freqs[0])
    out_header["CDELT3"] = float(df)
    out_header["CTYPE3"] = "FREQ"
    out_header["CUNIT3"] = "Hz"
    out_header["HISTORY"] = "stacked by skasim.imaging.stack_channels"

    fits.writeto(output_path, cube, out_header, overwrite=True)


# --------------------------------------------------------------------------- #
# image previews (FITS -> PNG via APLpy)
# --------------------------------------------------------------------------- #


def write_fits_preview(
    img_path: Path,
    png_path: Path,
    title: str,
    recenter: tuple[float, float, float] | None = None,
    scale_factor: float = 1000.0,
    bunit: str = "mJy/beam",
    colorbar_label: str = "mJy/beam",
) -> None:
    """Write a publication-style PNG preview for a WSClean FITS image, optionally recentered."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import aplpy
    import cmasher as cmr
    import matplotlib.pyplot as plt

    with fits.open(img_path) as source_hdul:
        source_hdu = source_hdul[0]
        data = np.asarray(source_hdu.data).squeeze()
        while data.ndim > 2:
            data = data[0]

        display_data = data * scale_factor
        finite = display_data[np.isfinite(display_data)]
        if finite.size:
            rms = float(np.nanstd(finite))
            vmin = float(np.nanmin(finite))
            vmax = float(np.nanmax(finite))
            if rms > 0:
                vmin = max(vmin, -2.0 * rms)
                vmax = min(vmax, 20.0 * rms)
        else:
            rms = 0.0
            vmin = vmax = None

        hdu = _make_2d_preview_hdu(display_data, source_hdu.header, bunit=bunit)
        hdul = fits.HDUList([hdu])
        fig = plt.figure(figsize=(8, 7))
        ffig = aplpy.FITSFigure(hdul, figure=fig)

        if recenter:
            ra_deg, dec_deg, fov_deg = recenter
            ffig.recenter(ra_deg, dec_deg, width=fov_deg, height=fov_deg)

        cmap = cmr.get_sub_cmap("cmr.rainforest", 0.30, 0.85)
        ffig.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
        if rms > 0 and finite.size:
            levels = 5.0 * rms * np.sqrt(3.0) ** np.arange(1, 25)
            drawable_levels = levels[levels <= np.nanmax(finite)]
            if drawable_levels.size:
                try:
                    ffig.show_contour(
                        hdul,
                        levels=drawable_levels,
                        colors="white",
                        linewidths=0.45,
                    )
                except AttributeError as exc:
                    logger.warning(f"Skipping FITS preview contours: {exc}")
        if "BMAJ" in hdu.header and "BMIN" in hdu.header:
            ffig.add_beam()
            ffig.beam.set_color("white")
            ffig.beam.set_edgecolor("black")
        ffig.axis_labels.set_xtext("RA")
        ffig.axis_labels.set_ytext("Dec")
        ffig.add_colorbar()
        ffig.colorbar.set_axis_label_text(colorbar_label)
        ffig.savefig(str(png_path), dpi=150)
        plt.close(fig)


def _make_2d_preview_hdu(
    data: np.ndarray,
    header: fits.Header,
    bunit: str | None = None,
) -> fits.PrimaryHDU:
    """Build an in-memory 2D celestial HDU suitable for APLpy rendering."""
    try:
        preview_header = WCS(header).celestial.to_header()
    except Exception:
        preview_header = fits.Header()
    for key in ("BMAJ", "BMIN", "BPA", "BUNIT", "OBJECT", "TELESCOP", "INSTRUME"):
        if key in header:
            preview_header[key] = header[key]
    if bunit is not None:
        preview_header["BUNIT"] = bunit
    return fits.PrimaryHDU(data=data, header=preview_header)


def write_sky_model_previews(
    sky_model,
    center: SkyCoord,
    fov: u.Quantity,
    work_dir: Path,
    run_id: str,
) -> list[tuple[str, str]]:
    """Write full and FoV sky-model source previews."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    sources = sky_model.to_json()
    if not sources:
        return []

    ra = np.asarray([src["ra"] for src in sources], dtype=float)
    dec = np.asarray([src["dec"] for src in sources], dtype=float)
    flux = np.asarray([src["I"] for src in sources], dtype=float)
    major_axis = np.asarray(
        [src.get("major_axis", 0.0) or 0.0 for src in sources],
        dtype=float,
    )
    minor_axis = np.asarray(
        [src.get("minor_axis", 0.0) or 0.0 for src in sources],
        dtype=float,
    )
    position_angle = np.asarray(
        [src.get("pa", 0.0) or 0.0 for src in sources], dtype=float
    )
    positive_flux = flux[flux > 0]
    norm = None
    if positive_flux.size:
        norm = LogNorm(
            vmin=float(np.nanmin(positive_flux)), vmax=float(np.nanmax(positive_flux))
        )

    full_name = f"{run_id}_sky_model.png"
    fov_name = f"{run_id}_sky_model_fov.png"
    _plot_sky_model_sources(
        work_dir / full_name,
        ra,
        dec,
        flux,
        major_axis,
        minor_axis,
        position_angle,
        norm,
        title=f"Sky model ({len(sources)} sources)",
    )
    half_fov = fov.to(u.deg).value / 2.0
    _plot_sky_model_sources(
        work_dir / fov_name,
        ra,
        dec,
        flux,
        major_axis,
        minor_axis,
        position_angle,
        norm,
        title=f"Sky model FoV ({fov.to(u.deg).value:.2f} deg)",
        xlim=(center.ra.deg + half_fov, center.ra.deg - half_fov),
        ylim=(center.dec.deg - half_fov, center.dec.deg + half_fov),
        fov_circle=(center.ra.deg, center.dec.deg, half_fov),
    )
    return [(full_name, "sky_model"), (fov_name, "sky_model_fov")]


def _plot_sky_model_sources(
    png_path: Path,
    ra: np.ndarray,
    dec: np.ndarray,
    flux: np.ndarray,
    major_axis: np.ndarray,
    minor_axis: np.ndarray,
    position_angle: np.ndarray,
    norm,
    title: str,
    xlim: tuple[float, float] | None = None,
    ylim: tuple[float, float] | None = None,
    fov_circle: tuple[float, float, float] | None = None,
) -> None:
    """Plot source positions as ellipses with astronomical RA orientation."""
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import PatchCollection

    fig, ax = plt.subplots(figsize=(7, 6), facecolor="white")
    ax.set_title(title)
    ax.set_xlabel("RA (deg)")
    ax.set_ylabel("Dec (deg)")
    ax.grid(True, color="0.85", linestyle=":", linewidth=0.8)
    if xlim is not None:
        ax.set_xlim(*xlim)
        plot_width_deg = abs(xlim[1] - xlim[0])
    else:
        ra_min, ra_max = _padded_limits(ra)
        ax.set_xlim(ra_max, ra_min)
        plot_width_deg = ra_max - ra_min
    if ylim is not None:
        ax.set_ylim(*ylim)
    else:
        ax.set_ylim(*_padded_limits(dec))
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_box_aspect(1)

    compact = _compact_source_mask(major_axis, plot_width_deg)
    resolved = ~compact
    if np.any(resolved):
        ellipses = _sky_model_ellipses(
            ra[resolved],
            dec[resolved],
            major_axis[resolved],
            minor_axis[resolved],
            position_angle[resolved],
        )
        ellipse_collection = PatchCollection(
            ellipses,
            cmap=SKY_MODEL_CMAP,
            norm=norm,
            alpha=0.82,
            edgecolor="black",
            linewidth=0.35,
        )
        ellipse_collection.set_array(flux[resolved])
        ax.add_collection(ellipse_collection)
    if np.any(compact):
        ax.scatter(
            ra[compact],
            dec[compact],
            s=_flux_marker_sizes(flux[compact]),
            c=flux[compact],
            cmap=SKY_MODEL_CMAP,
            norm=norm,
            marker="+",
            linewidths=1.2,
            alpha=0.9,
        )
    if fov_circle is not None:
        from matplotlib.patches import Circle

        ax.add_patch(
            Circle(
                (fov_circle[0], fov_circle[1]),
                fov_circle[2],
                fill=False,
                color="tab:red",
                linestyle="--",
                linewidth=1.2,
            )
        )
    scalar = ScalarMappable(norm=norm, cmap=SKY_MODEL_CMAP)
    scalar.set_array(flux)
    cbar = fig.colorbar(scalar, ax=ax)
    cbar.set_label("Stokes I (Jy)")
    fig.tight_layout()
    fig.savefig(png_path, dpi=140)
    plt.close(fig)


def _padded_limits(
    values: np.ndarray, pad_fraction: float = 0.05
) -> tuple[float, float]:
    """Return finite min/max limits with a small visual padding."""
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return (0.0, 1.0)
    lower = float(np.nanmin(finite))
    upper = float(np.nanmax(finite))
    span = upper - lower
    if span <= 0:
        span = max(abs(lower) * 0.01, 1.0 / 3600.0)
    pad = span * pad_fraction
    return (lower - pad, upper + pad)


def _compact_source_mask(
    major_axis_arcsec: np.ndarray,
    plot_width_deg: float,
) -> np.ndarray:
    """Return sources too small to read as ellipses at the plotted FoV."""
    threshold_arcsec = max(3.0, abs(plot_width_deg) * 3600.0 * 0.01)
    return major_axis_arcsec < threshold_arcsec


def _flux_marker_sizes(flux: np.ndarray) -> np.ndarray:
    """Map source flux densities to visible cross marker areas."""
    positive = flux[np.isfinite(flux) & (flux > 0)]
    if positive.size == 0:
        return np.full(flux.shape, 45.0)
    lo = float(np.nanmin(positive))
    hi = float(np.nanmax(positive))
    safe_flux = np.clip(flux, lo, hi)
    if hi <= lo:
        scaled = np.ones_like(safe_flux)
    else:
        scaled = (np.log10(safe_flux) - np.log10(lo)) / (np.log10(hi) - np.log10(lo))
    return 35.0 + scaled * 140.0


def _sky_model_position_angle(pa_deg: float) -> float:
    """Convert astronomical PA east of north to Matplotlib angle from +x."""
    return 90.0 - pa_deg


def _sky_model_ellipses(
    ra: np.ndarray,
    dec: np.ndarray,
    major_axis_arcsec: np.ndarray,
    minor_axis_arcsec: np.ndarray,
    position_angle_deg: np.ndarray,
) -> list:
    """Convert source shape metadata to Matplotlib ellipses in degree units."""
    from matplotlib.patches import Ellipse

    ellipses = []
    fallback_arcsec = 8.0
    for ra_deg, dec_deg, major, minor, pa in zip(
        ra,
        dec,
        major_axis_arcsec,
        minor_axis_arcsec,
        position_angle_deg,
    ):
        major = float(major) if np.isfinite(major) and major > 0 else fallback_arcsec
        minor = float(minor) if np.isfinite(minor) and minor > 0 else major
        ellipses.append(
            Ellipse(
                (float(ra_deg), float(dec_deg)),
                width=major / 3600.0,
                height=minor / 3600.0,
                angle=_sky_model_position_angle(float(pa)) if np.isfinite(pa) else 90.0,
            )
        )
    return ellipses
