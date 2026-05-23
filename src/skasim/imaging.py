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
    try:
        dirty_image.plot(
            title="Dirty image OSKAR",
            filename=str(dirty_png),
            wcs_enabled=True,
            xlabel="RA",
            ylabel="DEC",
            block=False,
        )
    finally:
        import matplotlib.pyplot as plt

        plt.close("all")
    dirty_image.write_to_file(str(dirty_fits), overwrite=True)
    logger.debug(f"Dirty PNG: {dirty_png}")
    logger.debug(f"Dirty FITS: {dirty_fits}")

    image_product_id = f"{work_dir.name}_dirty"
    ctx.manifest.add_output(
        "image_product",
        str(dirty_png.relative_to(work_dir)),
        image_product_id=image_product_id,
        imager="oskar-dirty",
        role="preview",
    )
    ctx.manifest.add_output(
        "image_product",
        str(dirty_fits.relative_to(work_dir)),
        image_product_id=image_product_id,
        imager="oskar-dirty",
        role="image",
    )


# --------------------------------------------------------------------------- #
# cleaned imaging (WSClean)
# --------------------------------------------------------------------------- #


def build_wsclean_argv(
    config: SimConfig,
    visibility_path: Path,
    fov: u.Quantity,
    output_prefix: str,
) -> list[str]:
    """Build a shell-free WSClean argv list from the resolved imaging config."""
    imaging_cellsize = fov / config.imaging.pixels
    channels_out = min(config.observation.n_channels or 1, 8)
    return shlex.split(config.imaging.wsclean_command) + [
        "-weight",
        "briggs",
        str(config.imaging.robust),
        "-multiscale",
        "-size",
        str(config.imaging.pixels),
        str(config.imaging.pixels),
        "-scale",
        f"{imaging_cellsize.to(u.arcsec).value:.6f}asec",
        "-niter",
        str(config.clean_iterations),
        "-mgain",
        "0.8",
        "-auto-threshold",
        "0.3",
        "-auto-mask",
        "3",
        "-channels-out",
        str(channels_out),
        "-join-channels",
        "-local-rms",
        "-name",
        output_prefix,
        str(visibility_path),
    ]


def run_wsclean_command(argv: list[str], work_dir: Path):
    """Run WSClean with argv and an explicit working directory."""
    env = os.environ.copy()
    env["OPENBLAS_NUM_THREADS"] = "1"
    return subprocess.run(
        argv,
        shell=False,
        cwd=str(work_dir),
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )


def collect_wsclean_outputs(work_dir: Path, output_prefix: str) -> list[Path]:
    """Collect WSClean FITS outputs for one configured output prefix."""
    return sorted(work_dir.glob(f"{output_prefix}*.fits"))


def wsclean_output_prefix(ctx: RunContext) -> str:
    """Return the stable WSClean output prefix for this run."""
    return f"{ctx.work_dir.name}_wsclean"


def write_fits_preview(img_path: Path, png_path: Path, title: str) -> None:
    """Write a non-interactive PNG preview for a WSClean FITS image."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt
    from matplotlib.colors import PowerNorm

    data = fits.getdata(img_path)
    data = data.squeeze()
    while data.ndim > 2:
        data = data[0]

    fig, ax = plt.subplots(figsize=(7, 6))
    finite = data[np.isfinite(data)]
    if finite.size:
        vmin, vmax = np.nanpercentile(finite, [1, 99])
        if vmin == vmax:
            vmin, vmax = None, None
    else:
        vmin, vmax = None, None
    image = ax.imshow(
        data,
        origin="lower",
        cmap="viridis",
        norm=PowerNorm(0.3, vmin=vmin, vmax=vmax),
    )
    ax.set_title(title)
    ax.set_xlabel("RA")
    ax.set_ylabel("DEC")
    fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    fig.savefig(png_path, dpi=130)
    plt.close(fig)


def run_wsclean_imaging(
    ctx: RunContext,
    visibility_path: Path,
    fov: u.Quantity,
) -> None:
    """produce cleaned image via external WSClean binary."""
    wsclean_module = require_karabo_module("karabo.imaging.imager_wsclean")
    file_handler_module = require_karabo_module("karabo.util.file_handler")
    config = ctx.config
    work_dir = ctx.work_dir

    output_prefix = wsclean_output_prefix(ctx)
    argv = build_wsclean_argv(config, visibility_path, fov, output_prefix=output_prefix)
    logger.info(f"WSClean command: {argv}")

    file_handler_module.FileHandler().get_tmp_dir(
        prefix=wsclean_module.TMP_PREFIX_CUSTOM,
        purpose=wsclean_module.TMP_PURPOSE_CUSTOM,
    )
    proc = run_wsclean_command(argv, work_dir)
    logger.info(f"WSClean stdout: {proc.stdout}")

    # remove the temporary files created by WSClean
    for tmp in work_dir.glob("wsclean-00*.fits"):
        if tmp.name.startswith(output_prefix):
            continue
        try:
            tmp.unlink()
        except Exception as exc:
            logger.error(show_exc(exc))

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
            png_path.name,
            image_product_id=output_prefix,
            imager="wsclean",
            role=f"{role}_preview",
        )
        ctx.manifest.add_output(
            "image_product",
            img_path.name,
            image_product_id=output_prefix,
            imager="wsclean",
            role=role,
        )
