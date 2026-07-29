"""WSClean -predict injection for spectral-cube models.

WSClean predicts visibilities from a set of per-channel model images and writes
them into the MODEL_DATA column of a Measurement Set.  This module splits a
resampled 3D spectral cube into the naming convention WSClean expects and
calls ``wsclean -predict``.
"""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.io import fits
from loguru import logger

from ..config import ImgConfig


def write_per_channel_model_fits(
    cube_data: np.ndarray,
    header_template: fits.Header,
    output_dir: Path,
    prefix: str,
    freq_axis: int = 3,
) -> list[Path]:
    """Split a 3D spectral cube into WSClean per-channel model images.

    Parameters
    ----------
    cube_data
        Numpy array in axis order (freq, dec, ra) after any reordering.
    header_template
        FITS header describing the full cube.  Must contain CTYPE1=RA,
        CTYPE2=DEC and a frequency axis.
    output_dir
        Directory where ``<prefix>-NNNN-model.fits`` files will be written.
    prefix
        WSClean prefix, e.g. ``model_entry_01_spectral_cube``.
    freq_axis
        1-based FITS axis index for frequency in ``header_template``.

    Returns
    -------
    paths
        Ordered list of per-channel FITS paths.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    n_freq = cube_data.shape[0]
    crpix = float(header_template.get(f"CRPIX{freq_axis}", 1.0))  # type: ignore[arg-type]
    crval_key = f"CRVAL{freq_axis}"
    cdelt_key = f"CDELT{freq_axis}"
    cunit_key = f"CUNIT{freq_axis}"

    if crval_key not in header_template or cdelt_key not in header_template:
        raise ValueError(
            f"header_template missing {crval_key} or {cdelt_key}; "
            "cannot build per-channel WSClean model headers."
        )

    crval = float(header_template[crval_key])  # type: ignore[arg-type]
    cdelt = float(header_template[cdelt_key])  # type: ignore[arg-type]
    cunit = str(header_template.get(cunit_key) or "").strip().lower()
    if cunit == "mhz":
        crval *= 1e6
        cdelt *= 1e6

    out_paths: list[Path] = []
    for i in range(n_freq):
        hdr = header_template.copy()
        # strip old NAXIS* keywords and rebuild a 3D single-plane header
        for key in list(hdr.keys()):
            if key.startswith("NAXIS") and key != "NAXIS":
                hdr.remove(key)
        hdr["NAXIS"] = 3
        hdr["NAXIS1"] = cube_data.shape[2]
        hdr["NAXIS2"] = cube_data.shape[1]
        hdr["NAXIS3"] = 1

        for k in (
            "CTYPE1", "CRPIX1", "CRVAL1", "CDELT1", "CUNIT1",
            "CTYPE2", "CRPIX2", "CRVAL2", "CDELT2", "CUNIT2",
        ):
            if k in header_template:
                hdr[k] = header_template[k]

        hdr["CTYPE3"] = "FREQ"
        hdr["CRPIX3"] = 1.0
        hdr["CRVAL3"] = crval + (i - (crpix - 1.0)) * cdelt
        hdr["CDELT3"] = cdelt
        hdr["CUNIT3"] = "Hz"
        hdr["BUNIT"] = header_template.get("BUNIT", "Jy/pixel")

        plane = cube_data[i, :, :][np.newaxis, :, :]  # shape (1, dec, ra)
        path = output_dir / f"{prefix}-{i:04d}-model.fits"
        fits.writeto(path, plane, hdr, overwrite=True)
        out_paths.append(path)

    logger.info(
        f"wrote {len(out_paths)} per-channel model FITS for WSClean predict "
        f"under {output_dir}"
    )
    return out_paths


def _imaging_cellsize(fov_deg: float | None, pixels: int) -> u.Quantity:
    """Return angular pixel size from FoV and number of pixels."""
    fov = float(fov_deg if fov_deg is not None else 1.0) * u.deg
    return fov / pixels


def build_wsclean_predict_argv(
    wsclean_command: str,
    visibility_path: Path,
    img_config: ImgConfig,
    prefix: str,
    n_channels: int,
    pixel_size_arcsec: float | None = None,
    n_pixels: int | None = None,
) -> list[str]:
    """Build argv for ``wsclean -predict``.

    Only the executable from ``wsclean_command`` is kept; imaging/deconvolution
    flags (``-niter``, ``-mgain``, ``-auto-mask``, etc.) are stripped because
    they are invalid in predict mode.  Use ``wsclean_predict_command`` if you
    need custom predict-only flags.
    """
    if pixel_size_arcsec is None or n_pixels is None:
        imaging_cellsize = _imaging_cellsize(img_config.fov_deg, img_config.pixels)
        pixel_size_arcsec = imaging_cellsize.to(u.arcsec).value
        n_pixels = img_config.pixels

    # Keep only the executable from the supplied command line; imaging flags
    # are not applicable to ``wsclean -predict`` and make it fail.
    tokens = shlex.split(wsclean_command)
    executable = tokens[0] if tokens else "wsclean"

    argv = [executable] + [
        "-predict",
        "-gridder",
        "wgridder",
        "-wgridder-accuracy",
        "1e-5",
        "-size",
        str(n_pixels),
        str(n_pixels),
        "-scale",
        f"{pixel_size_arcsec:.6f}asec",
        "-channels-out",
        str(n_channels),
        "-name",
        prefix,
    ]
    if img_config.threads is not None:
        argv += ["-j", str(img_config.threads)]
    argv.append(str(visibility_path))
    return argv


def run_wsclean_predict(
    wsclean_command: str,
    visibility_path: Path,
    img_config: ImgConfig,
    prefix: str,
    n_channels: int,
    work_dir: Path,
    pixel_size_arcsec: float | None = None,
    n_pixels: int | None = None,
) -> None:
    """Run WSClean in predict mode to fill MODEL_DATA of the MS."""
    argv = build_wsclean_predict_argv(
        wsclean_command,
        visibility_path,
        img_config,
        prefix,
        n_channels,
        pixel_size_arcsec=pixel_size_arcsec,
        n_pixels=n_pixels,
    )
    logger.info(f"WSClean predict command: {argv}")

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
        if proc.stdout is None:
            raise RuntimeError("subprocess.Popen returned None stdout with PIPE")
        for line in proc.stdout:
            stripped = line.rstrip("\n")
            logger.info("[wsclean-predict] {}", stripped)
            lines.append(stripped)
        proc.wait()
    combined = "\n".join(lines)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode,
            argv,
            output=combined,
            stderr=None,
        )


def merge_model_data_into_data(visibility_path: Path) -> None:
    """Add image-model MODEL_DATA into the delivered DATA column.

    This is a thin compatibility wrapper around the canonical implementation in
    ``src/skasim/loaders/image_models.py``.  It exists so that legacy imports
    from ``skasim.loaders.wsclean_predict`` keep working; new code should import
    from ``skasim.loaders`` or ``skasim.loaders.image_models`` instead.
    """
    from .image_models import merge_model_data_into_data as _canonical_merge

    return _canonical_merge(visibility_path)


def inject_spectral_cube_with_wsclean_predict(
    ctx,
    entry,
    index: int,
    visibility_path: Path,
    img_config: ImgConfig,
    cube_data: np.ndarray,
    header: fits.Header,
    freq_axis: int,
) -> dict:
    """Write per-channel FITS and run ``wsclean -predict`` for a spectral cube."""
    prefix = f"model_entry_{index + 1:02d}_spectral_cube"
    model_dir = ctx.work_dir / f"{prefix}_predict_models"

    model_paths = write_per_channel_model_fits(
        cube_data, header, model_dir, prefix, freq_axis=freq_axis
    )

    run_wsclean_predict(
        wsclean_command=img_config.wsclean_command,
        visibility_path=visibility_path,
        img_config=img_config,
        prefix=prefix,
        n_channels=cube_data.shape[0],
        work_dir=model_dir,
    )

    return {
        "model_entry_index": index,
        "model_type": entry.type,
        "backend": "wsclean_predict",
        "n_channels": len(model_paths),
        "prefix": prefix,
        "model_dir": str(model_dir),
    }


def inject_continuum_i_alpha_with_wsclean_predict(
    ctx,
    entry,
    index: int,
    visibility_path: Path,
    img_config: ImgConfig,
    product,
) -> dict:
    """Inject a continuum I+alpha model using ``wsclean -predict``.

    The CASA Taylor-term images (tt0, tt1) produced for the model are exported to
    per-channel FITS images following WSClean's ``<prefix>-NNNN-model.fits``
    convention.  Each channel is computed from the first-order Taylor expansion
    around the observation reference frequency:

        I(ν) = tt0 + tt1 * (ν - ν0) / ν0

    where ν0 is the observation band centre.  ``wsclean -predict`` then fills
    the MODEL_DATA column of the visibility MS.
    """
    if product.nterms != 2:
        raise ValueError(
            "inject_continuum_i_alpha_with_wsclean_predict requires a 2-term "
            "Taylor product (tt0 and tt1)."
        )

    prefix = f"model_entry_{index + 1:02d}_continuum_predict"
    model_dir = ctx.work_dir / f"{prefix}_models"
    model_dir.mkdir(parents=True, exist_ok=True)

    # Use the intermediate FITS files (4D with degenerate STOKES/FREQ axes) so we
    # can read tt0/tt1 without requiring a CASA image export step.
    fits_paths = getattr(product, "intermediates", None)
    if not fits_paths or len(fits_paths) != 2:
        raise ValueError(
            "continuum_i_alpha product must expose two intermediate FITS paths "
            "(tt0.fits and tt1.fits)."
        )

    tt0_fits = fits_paths[0]
    tt1_fits = fits_paths[1]

    with fits.open(tt0_fits) as hdul0, fits.open(tt1_fits) as hdul1:
        tt0_data = np.asarray(hdul0[0].data, dtype=float).squeeze()
        tt1_data = np.asarray(hdul1[0].data, dtype=float).squeeze()
        header_template = hdul0[0].header.copy()

    if tt0_data.ndim != 2 or tt1_data.ndim != 2:
        raise ValueError(
            f"intermediate tt0 and tt1 must squeeze to 2D spatial arrays; "
            f"got tt0 {tt0_data.shape}, tt1 {tt1_data.shape}"
        )

    if tt0_data.shape != tt1_data.shape:
        raise ValueError(
            f"tt0 and tt1 shapes differ: {tt0_data.shape} vs {tt1_data.shape}"
        )

    # Derive the pixel geometry from the model header so WSClean predict matches.
    cdelt1 = abs(float(header_template.get("CDELT1", 0.0)))
    cdelt2 = abs(float(header_template.get("CDELT2", 0.0)))
    if cdelt1 <= 0.0 or cdelt2 <= 0.0:
        raise ValueError(
            f"model header missing CDELT1/CDELT2; cannot determine pixel size"
        )
    # Average pixel size in arcseconds (FITS CDELT is in degrees).
    pixel_size_deg = 0.5 * (cdelt1 + cdelt2)
    pixel_size_arcsec = pixel_size_deg * 3600.0

    ny, nx = tt0_data.shape
    if ny != nx:
        logger.warning(
            "continuum_i_alpha model image is not square ({}×{}); WSClean predict "
            "will use {} pixels for both dimensions, which may distort the model.",
            nx, ny, nx,
        )
    n_pixels = nx

    obs = ctx.config.observation
    n_channels = obs.n_channels
    center_hz = obs.frequency_mhz * 1e6
    bandwidth_hz = obs.bandwidth_mhz * 1e6
    chan_width_hz = bandwidth_hz / n_channels
    start_hz = center_hz - bandwidth_hz / 2.0

    model_paths: list[Path] = []
    for i in range(n_channels):
        freq_hz = start_hz + (i + 0.5) * chan_width_hz
        # First-order Taylor expansion around the observation reference.
        relative = (freq_hz - center_hz) / center_hz
        plane = tt0_data + tt1_data * relative

        hdr = header_template.copy()
        # Strip old NAXIS* keywords and rebuild a 3D single-plane header.
        for key in list(hdr.keys()):
            if key.startswith("NAXIS") and key != "NAXIS":
                hdr.remove(key)
        hdr["NAXIS"] = 3
        hdr["NAXIS1"] = plane.shape[-1]
        hdr["NAXIS2"] = plane.shape[-2]
        hdr["NAXIS3"] = 1

        for k in (
            "CTYPE1", "CRPIX1", "CRVAL1", "CDELT1", "CUNIT1",
            "CTYPE2", "CRPIX2", "CRVAL2", "CDELT2", "CUNIT2",
        ):
            if k in header_template:
                hdr[k] = header_template[k]

        hdr["CTYPE3"] = "FREQ"
        hdr["CRPIX3"] = 1.0
        hdr["CRVAL3"] = freq_hz
        hdr["CDELT3"] = chan_width_hz
        hdr["CUNIT3"] = "Hz"
        hdr["BUNIT"] = header_template.get("BUNIT", "Jy/pixel")

        path = model_dir / f"{prefix}-{i:04d}-model.fits"
        fits.writeto(path, plane[np.newaxis, :, :], hdr, overwrite=True)
        model_paths.append(path)

    logger.info(
        f"wrote {len(model_paths)} per-channel model FITS for continuum I+alpha "
        f"WSClean predict under {model_dir}"
    )

    wsclean_command = img_config.wsclean_predict_command or img_config.wsclean_command

    run_wsclean_predict(
        wsclean_command=wsclean_command,
        visibility_path=visibility_path,
        img_config=img_config,
        prefix=prefix,
        n_channels=n_channels,
        work_dir=model_dir,
        pixel_size_arcsec=pixel_size_arcsec,
        n_pixels=n_pixels,
    )

    return {
        "model_entry_index": index,
        "model_type": entry.type,
        "backend": "wsclean_predict",
        "n_channels": len(model_paths),
        "prefix": prefix,
        "model_dir": str(model_dir),
        "pixel_size_arcsec": pixel_size_arcsec,
        "n_pixels": n_pixels,
        "wsclean_command": wsclean_command,
    }


__all__ = [
    "build_wsclean_predict_argv",
    "inject_continuum_i_alpha_with_wsclean_predict",
    "inject_spectral_cube_with_wsclean_predict",
    "merge_model_data_into_data",
    "run_wsclean_predict",
    "write_per_channel_model_fits",
]
