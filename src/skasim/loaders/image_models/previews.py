"""Image-model preview generation for the weblog sky-model section."""

from __future__ import annotations

import warnings
from pathlib import Path

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from loguru import logger

from ...config import CasaTaylorTermsModelEntry, SpectralCubeModelEntry
from ...imaging import write_fits_preview
from ...manifest import RunContext
from .casa_interop import run_casa_exportfits
from .fits_io import (
    _find_frequency_axis,
    _fits_axis_to_numpy,
    _freq_axis_centres,
    _squeeze_degenerate_axes,
    _strip_spectral_axis_from_header,
    image_model_entries,
    primary_model_fits_path,
    read_fits_image_info,
)

# suppress fits formatting fixes
warnings.simplefilter("ignore", category=FITSFixedWarning)
# suppress polar motion fallback warnings
warnings.filterwarnings("ignore", message=".*polar motions.*")


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
    freq_axis = _freq_axis_centres(header, nchan, axis=3)
    cunit3 = (header.get("CUNIT3") or "Hz").strip()
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

    base_header = _strip_spectral_axis_from_header(header)

    out_fits = work_dir / f"{output_prefix}-moment8.fits"
    base_header["BUNIT"] = header.get("BUNIT") or "Jy/beam"
    base_header["MOMENT"] = 8
    base_header["HISTORY"] = "produced by skasim.run_moment8_for_spectral_cube"
    if out_fits.exists():
        out_fits.unlink()
    fits.writeto(
        out_fits, np.asarray(moment8, dtype=np.float32), base_header, overwrite=True
    )
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
        out_header = _strip_spectral_axis_from_header(header)
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
