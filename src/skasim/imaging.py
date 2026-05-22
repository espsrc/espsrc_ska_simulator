"""imaging.py — dirty (OSKAR) and cleaned (WSClean) imaging wrappers."""

from __future__ import annotations

import glob
import os
import shutil
import subprocess
from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
from loguru import logger

from .config import SimConfig
from .manifest import RunContext
from .runtime import require_karabo_module
from .utils import show_exc

# --------------------------------------------------------------------------- #
# dirty imaging (OSKAR)
# --------------------------------------------------------------------------- #


def run_dirty_imaging(
    ctx: RunContext,
    visibility_path: Path,
    fov: u.Quantity,
    center: SkyCoord,
) -> None:
    """produce dirty image via OSKAR."""
    imager_module = require_karabo_module("karabo.imaging.imager_oskar")
    visibility_module = require_karabo_module("karabo.simulation.visibility")
    config = ctx.config
    work_dir = ctx.work_dir
    vis = visibility_module.Visibility(str(visibility_path))
    imaging_cellsize = fov / config.imaging.pixels
    cfg = imager_module.OskarDirtyImagerConfig(
        imaging_npixel=config.imaging.pixels,
        imaging_cellsize=imaging_cellsize.to(u.rad).value,
        combine_across_frequencies=True,
        imaging_phase_centre=center,
    )
    imager = imager_module.OskarDirtyImager(config=cfg)
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

    ctx.manifest.outputs.extend([
        str(dirty_png.relative_to(work_dir)),
        str(dirty_fits.relative_to(work_dir)),
    ])


# --------------------------------------------------------------------------- #
# cleaned imaging (WSClean)
# --------------------------------------------------------------------------- #


def run_wsclean_imaging(
    ctx: RunContext,
    visibility_path: Path,
    fov: u.Quantity,
) -> None:
    """produce cleaned image via external WSClean binary."""
    image_module = require_karabo_module("karabo.imaging.image")
    wsclean_module = require_karabo_module("karabo.imaging.imager_wsclean")
    file_handler_module = require_karabo_module("karabo.util.file_handler")
    config = ctx.config
    work_dir = ctx.work_dir

    # switch to work_dir so WSClean outputs and glob-based file ops
    # resolve relative to the job directory, not the caller's CWD
    orig_cwd = Path.cwd()
    os.chdir(str(work_dir))

    try:
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

        file_handler_module.FileHandler().get_tmp_dir(
            prefix=wsclean_module.TMP_PREFIX_CUSTOM,
            purpose=wsclean_module.TMP_PURPOSE_CUSTOM,
        )
        expected_prefix = f"{wsclean_module._WSCLEAN_BINARY} "
        if not custom_command.startswith(expected_prefix):
            raise ValueError(f"Command must start with '{expected_prefix}'")

        cmd = f"OPENBLAS_NUM_THREADS=1 {custom_command}"
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
        logger.info(f"WSClean stdout: {proc.stdout}")

        # remove the temporary files created by WSClean
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
            img = image_module.Image(path=img_path)
            png_name = f"{work_dir.name}_{img_path.replace('.fits', '.png')}"

            # infer image type from filename for correct plot title
            title = "Imaging output (WSClean)"
            if "MFS-image" in img_path:
                title = "Cleaned image (WSClean)"
            elif "MFS-model" in img_path:
                title = "Component model (WSClean)"
            elif "MFS-residual" in img_path:
                title = "Residual (WSClean)"
            elif "MFS-dirty" in img_path:
                title = "Dirty image (WSClean)"
            elif "MFS-psf" in img_path:
                title = "Point spread function (WSClean)"

            img.plot(
                title=title,
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

            ctx.manifest.outputs.extend([png_name, new_name])

    finally:
        os.chdir(str(orig_cwd))
