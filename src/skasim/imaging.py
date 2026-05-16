"""imaging.py — dirty (OSKAR) and cleaned (WSClean) imaging wrappers."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

import astropy.units as u
from astropy.coordinates import SkyCoord
from karabo.imaging.image import Image
from karabo.imaging.imager_oskar import OskarDirtyImager, OskarDirtyImagerConfig
from karabo.imaging.imager_wsclean import (
    _WSCLEAN_BINARY,
    TMP_PREFIX_CUSTOM,
    TMP_PURPOSE_CUSTOM,
)
from karabo.simulation.visibility import Visibility
from karabo.util.file_handler import FileHandler
from loguru import logger

from .config import SimConfig
from .utils import show_exc

# --------------------------------------------------------------------------- #
# dirty imaging (OSKAR)
# --------------------------------------------------------------------------- #


def run_dirty_imaging(
    config: SimConfig,
    visibility_path: Path,
    fov: u.Quantity,
    center: SkyCoord,
    work_dir: Path,
) -> None:
    """produce dirty image via OSKAR."""
    vis = Visibility(str(visibility_path))
    imaging_cellsize = fov / config.imaging.pixels
    cfg = OskarDirtyImagerConfig(
        imaging_npixel=config.imaging.pixels,
        imaging_cellsize=imaging_cellsize.to(u.rad).value,
        combine_across_frequencies=True,
        imaging_phase_centre=center,
    )
    imager = OskarDirtyImager(config=cfg)
    dirty_image = imager.create_dirty_image(vis)

    dirty_png = work_dir / f"{work_dir.name}_dirty.png"
    dirty_fits = work_dir / f"{work_dir.name}_dirty.fits"
    dirty_image.plot(
        title="Dirty image OSKAR",
        filename=str(dirty_png),
        wcs_enabled=True,
        xlabel="RA",
        ylabel="DEC",
    )
    dirty_image.write_to_file(str(dirty_fits), overwrite=True)
    logger.debug(f"Dirty PNG: {dirty_png}")
    logger.debug(f"Dirty FITS: {dirty_fits}")


# --------------------------------------------------------------------------- #
# cleaned imaging (WSClean)
# --------------------------------------------------------------------------- #


def run_wsclean_imaging(
    config: SimConfig,
    visibility_path: Path,
    fov: u.Quantity,
    work_dir: Path,
) -> None:
    """produce cleaned image via external WSClean binary."""
    imaging_cellsize = fov / config.imaging.pixels
    threshold = "--auto-threshold 0.3"
    custom_command = (
        f"wsclean -weight briggs {config.imaging.robust} -multiscale "
        f"-size {config.imaging.pixels} {config.imaging.pixels} "
        f"-scale {imaging_cellsize.to(u.arcsec).value:.6f}asec "
        f"-niter {config.niter} -mgain 0.8 {threshold} "
        f"-auto-mask 3 -channels-out 8 -join-channels -local-rms "
        f"{visibility_path}"
    )
    logger.info(f"WSClean command: {custom_command}")

    FileHandler().get_tmp_dir(
        prefix=TMP_PREFIX_CUSTOM,
        purpose=TMP_PURPOSE_CUSTOM,
    )
    expected_prefix = f"{_WSCLEAN_BINARY} "
    if not custom_command.startswith(expected_prefix):
        raise ValueError(f"Command must start with '{expected_prefix}'")

    cmd = f"OPENBLAS_NUM_THREADS=1 {custom_command}"
    proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    logger.info(f"WSClean stdout: {proc.stdout}")

    # remove the temporary files created by WSClea
    for tmp in glob.glob("wsclean-00*.fits"):
        try:
            os.remove(tmp)
        except Exception as exc:
            logger.error(show_exc(exc))

    mfs_files = glob.glob("*-MFS-*.fits")
    logger.info(f"MFS files: {mfs_files}")

    from matplotlib.colors import PowerNorm

    gamma = 0.3
    for img_path in glob.glob("wsclean-*.fits"):
        img = Image(path=img_path)
        png_name = f"{work_dir.name}_{img_path.replace('.fits', '.png')}"
        img.plot(
            title="Cleaned image (WSCLEAN)",
            filename=png_name,
            wcs_enabled=True,
            xlabel="RA",
            ylabel="DEC",
            norm=PowerNorm(gamma),
        )
        new_name = (
            f"{work_dir.name}_bw{config.observation.bandwidth_mhz:.0f}_"
            f"ch{config.observation.n_channels}_fr{config.observation.freq_mhz:.0f}_"
            f"sec{config.observation.seconds}{img_path.replace('wsclean-', '')}"
        )
        shutil.move(img_path, new_name)
        logger.debug(f"Renamed {img_path} -> {new_name}")
