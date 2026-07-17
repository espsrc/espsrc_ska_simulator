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
) -> list[str]:
    """Build argv for ``wsclean -predict``."""
    imaging_cellsize = _imaging_cellsize(img_config.fov_deg, img_config.pixels)

    argv = shlex.split(wsclean_command) + [
        "-predict",
        "-gridder",
        "wgridder",
        "-wgridder-accuracy",
        "1e-5",
        "-size",
        str(img_config.pixels),
        str(img_config.pixels),
        "-scale",
        f"{imaging_cellsize.to(u.arcsec).value:.6f}asec",
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
) -> None:
    """Run WSClean in predict mode to fill MODEL_DATA of the MS."""
    argv = build_wsclean_predict_argv(
        wsclean_command,
        visibility_path,
        img_config,
        prefix,
        n_channels,
    )
    logger.info(f"WSClean predict command: {argv}")

    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"

    with subprocess.Popen(
        argv,
        shell=False,
        cwd=str(work_dir),
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    ) as proc:
        assert proc.stdout is not None
        for line in proc.stdout:
            logger.info("[wsclean-predict] {}", line.rstrip("\n"))
        proc.wait()
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(
            proc.returncode,
            argv,
            output="",
            stderr=None,
        )


def merge_model_data_into_data(visibility_path: Path) -> None:
    """Add MODEL_DATA into the delivered DATA column."""
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
                f"{visibility_path} must contain DATA and MODEL_DATA columns."
            )
        data = ms_table.getcol("DATA")
        model_data = ms_table.getcol("MODEL_DATA")
        ms_table.putcol("DATA", data + model_data)


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

__all__ = [
    "build_wsclean_predict_argv",
    "inject_spectral_cube_with_wsclean_predict",
    "merge_model_data_into_data",
    "run_wsclean_predict",
    "write_per_channel_model_fits",
]
