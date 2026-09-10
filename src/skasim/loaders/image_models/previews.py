"""Image-model preview generation for the weblog sky-model section."""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Dict, Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from loguru import logger

from ...config import CasaTaylorTermsModelEntry, SpectralCubeModelEntry
from ...imaging import _make_2d_preview_hdu
from ...manifest import RunContext
from .casa_interop import run_casa_exportfits
from .fits_io import (
    FitsImageInfo,
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
        # keep the existing single-panel moment-8 preview untouched
        from ...imaging import write_fits_preview

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

        from ...imaging import write_fits_preview

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
    center: Optional[SkyCoord] = None,
) -> None:
    """Write two-panel FITS model previews for the weblog sky-model section.

    Each continuum/Stokes-map/CASA-Taylor-term entry produces a side-by-side
    PNG: LEFT shows the FITS at its natural image extent so small models are
    legible; RIGHT shows the simulation FoV contextual view with the primary
    beam as a circle and the image-model footprint as a rectangle. Spectral
    cubes keep the existing single-panel moment-8 preview.

    Parameters
    ----------
    ctx
        Run context with the model configuration and manifest.
    fov
        Simulation field of view used for the right-hand contextual panel.
    center
        Simulation phase centre. When ``None``, the contextual panel falls back
        to the FITS model's own centre if it has a usable WCS.
    """
    entries = image_model_entries(ctx.config)
    if not entries:
        return

    for index, entry in enumerate(entries, start=1):
        # spectral cubes keep their existing moment-8 preview
        if isinstance(entry, SpectralCubeModelEntry):
            continue

        image_path = primary_model_fits_path(entry)
        export_path = None
        if image_path is None and isinstance(entry, CasaTaylorTermsModelEntry):
            image_path = Path(entry.tt0).expanduser().resolve()
            export_path = ctx.work_dir / f"model_entry_{index:02d}_casa_taylor.tt0.fits"
        if image_path is None:
            continue

        info = None
        try:
            info = read_fits_image_info(image_path)
        except Exception as exc:
            logger.debug(
                f"write_image_model_previews: failed to read {image_path}: {exc}"
            )
            # for casa taylor terms, the input may be a casa table; still try the export path
            if export_path is not None:
                pass
            else:
                continue

        suffix = "" if len(entries) == 1 else f"_{index:02d}"
        png_name = f"{ctx.work_dir.name}_fits_model{suffix}.png"
        png_path = ctx.work_dir / png_name
        preview_source = image_path
        if export_path is not None:
            run_casa_exportfits(ctx.work_dir, image_path, export_path)
            preview_source = export_path

        # if reading the original failed, read the exported FITS for metadata only
        if info is None and export_path is not None:
            try:
                info = read_fits_image_info(export_path)
            except Exception as exc:
                logger.debug(
                    f"write_image_model_previews: failed to read exported {export_path}: {exc}"
                )

        _write_image_model_two_panel_preview(
            preview_source,
            png_path,
            fov,
            title=f"FITS Model ({entry.type})",
            center=center,
        )
        footprint = _footprint_stats(preview_source, info, center, fov)
        ctx.manifest.add_output(
            "plot",
            png_name,
            role="fits_model",
            metadata={
                "model_entry_index": index - 1,
                "model_type": entry.type,
                "source_fits": str(image_path),
                "preview_fits": str(preview_source),
                **footprint,
            },
        )

    # For spectral-cube inputs, also render a moment-8 (peak) preview of the raw cube.
    write_spectral_cube_input_preview(ctx, fov)


def _write_image_model_two_panel_preview(
    img_path: Path,
    png_path: Path,
    fov: u.Quantity,
    title: str = "FITS Model",
    center: Optional[SkyCoord] = None,
    scale_factor: float = 1000.0,
    bunit: str = "mJy/pixel",
    colorbar_label: str = "mJy/pixel",
) -> None:
    """Render a two-panel preview: natural FITS extent (left) and simulation FoV context (right)."""
    import matplotlib

    matplotlib.use("Agg", force=True)
    import aplpy
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

        model_extent = _fits_image_extent_deg(source_hdu.header)
        model_center = _fits_image_center(source_hdu.header)
        if model_center is None and center is not None:
            # when the FITS has no WCS, the contextual panel degrades to the
            # simulation centre plus a small fallback footprint
            model_center = center

        fig = plt.figure(figsize=(14, 6))

        # left panel: natural image extent
        left = aplpy.FITSFigure(hdul, figure=fig, subplot=[0.06, 0.12, 0.43, 0.78])
        if model_center is not None and model_extent is not None:
            width_deg = model_extent["width_deg"] * 1.05
            height_deg = model_extent["height_deg"] * 1.05
            left.recenter(
                model_center.ra.deg,
                model_center.dec.deg,
                width=width_deg,
                height=height_deg,
            )
        cmap = _preview_cmap()
        left.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax)
        left.axis_labels.set_xtext("RA")
        left.axis_labels.set_ytext("Dec")
        left.set_title("Natural extent")
        left.add_colorbar()
        left.colorbar.set_axis_label_text(colorbar_label)
        if rms > 0 and finite.size:
            levels = 5.0 * rms * np.sqrt(3.0) ** np.arange(1, 25)
            drawable_levels = levels[levels <= np.nanmax(finite)]
            if drawable_levels.size:
                try:
                    left.show_contour(
                        hdul,
                        levels=drawable_levels,
                        colors="white",
                        linewidths=0.45,
                    )
                except AttributeError as exc:
                    logger.warning(f"Skipping FITS preview contours: {exc}")

        # right panel: contextual FoV view in world coordinates. This deliberately
        # does not draw the FITS raster: it is a geometry panel, not a second crop.
        right = fig.add_axes([0.54, 0.12, 0.43, 0.78])
        ctx_center = center if center is not None else model_center
        ctx_fov_deg = fov.to(u.deg).value
        if ctx_center is not None:
            _draw_fov_annotations(
                right,
                ctx_center.ra.deg,
                ctx_center.dec.deg,
                ctx_fov_deg,
                model_center,
                model_extent,
            )
            half_fov = ctx_fov_deg / 2.0
            right.set_xlim(ctx_center.ra.deg + half_fov, ctx_center.ra.deg - half_fov)
            right.set_ylim(ctx_center.dec.deg - half_fov, ctx_center.dec.deg + half_fov)
        else:
            right.text(
                0.5,
                0.5,
                "No WCS: FoV context unavailable",
                ha="center",
                va="center",
                transform=right.transAxes,
            )
        right.set_aspect("equal", adjustable="box")
        right.set_xlabel("RA (deg)")
        right.set_ylabel("Dec (deg)")
        right.set_title(f"FoV context = {ctx_fov_deg:.4f}°")
        right.grid(True, color="0.85", linestyle=":", linewidth=0.8)
        handles, labels = right.get_legend_handles_labels()
        if handles:
            right.legend(loc="upper right", fontsize="small")

        fig.suptitle(title, fontsize=12)
        fig.savefig(str(png_path), dpi=150)
        plt.close(fig)


def _fits_image_extent_deg(header: fits.Header) -> Optional[Dict[str, float]]:
    """Return approximate RA/Dec extent of a 2D FITS image in degrees.

    Falls back to pixel scale when a celestial WCS is not available.
    """
    naxis1 = int(header.get("NAXIS1") or 1)
    naxis2 = int(header.get("NAXIS2") or 1)
    try:
        wcs = WCS(header).celestial
        pix = np.array([[0.0, 0.0], [naxis1 - 1.0, naxis2 - 1.0]])
        sky = wcs.pixel_to_world(pix[:, 0], pix[:, 1])
        ra_span = float(abs(sky[0].ra.deg - sky[1].ra.deg))
        dec_span = float(abs(sky[0].dec.deg - sky[1].dec.deg))
        return {"width_deg": ra_span, "height_deg": dec_span}
    except Exception:
        cdelt1 = abs(float(header.get("CDELT1") or header.get("CD1_1") or 0.0))
        cdelt2 = abs(float(header.get("CDELT2") or header.get("CD2_2") or 0.0))
        if cdelt1 <= 0 or cdelt2 <= 0:
            return None
        return {
            "width_deg": cdelt1 * naxis1,
            "height_deg": cdelt2 * naxis2,
        }


def _fits_image_center(header: fits.Header) -> Optional[SkyCoord]:
    """Return the celestial centre of a 2D FITS image, or None if WCS is unavailable."""
    naxis1 = int(header.get("NAXIS1") or 1)
    naxis2 = int(header.get("NAXIS2") or 1)
    try:
        wcs = WCS(header).celestial
        center = wcs.pixel_to_world((naxis1 - 1) / 2.0, (naxis2 - 1) / 2.0)
        if isinstance(center, SkyCoord):
            return center
    except Exception as exc:
        logger.debug(f"_fits_image_center: WCS failed: {exc}")
    return None


def _footprint_stats(
    preview_source: Path,
    info: Optional[FitsImageInfo],
    center: Optional[SkyCoord],
    fov: u.Quantity,
) -> Dict[str, object]:
    """Summarize a FITS model's on-sky footprint for the weblog's Sky Model panel.

    Computed from the same header the two-panel preview renders from, so the
    reported numbers always match what the plot shows.
    """
    stats: Dict[str, object] = {}
    try:
        header = fits.getheader(preview_source)
    except Exception as exc:
        logger.debug(f"_footprint_stats: failed to read {preview_source}: {exc}")
        return stats

    extent = _fits_image_extent_deg(header)
    model_center = _fits_image_center(header)

    if info is not None:
        ny, nx = info.spatial_shape
        stats["image_size"] = [nx, ny]
        if extent is not None and nx and ny:
            stats["pixel_scale_arcsec"] = (
                (extent["width_deg"] / nx + extent["height_deg"] / ny) / 2.0 * 3600.0
            )

    if extent is not None:
        stats["extent_width_deg"] = extent["width_deg"]
        stats["extent_height_deg"] = extent["height_deg"]
        fov_deg = fov.to(u.deg).value
        if fov_deg > 0:
            stats["fov_fraction"] = (
                max(extent["width_deg"], extent["height_deg"]) / fov_deg
            )

    if model_center is not None:
        stats["center_ra_deg"] = model_center.ra.deg
        stats["center_dec_deg"] = model_center.dec.deg
        if center is not None:
            stats["offset_from_center_arcsec"] = model_center.separation(center).arcsec

    return stats


def _draw_fov_annotations(
    ax,
    center_ra: float,
    center_dec: float,
    fov_deg: float,
    model_center: Optional[SkyCoord],
    model_extent: Optional[Dict[str, float]],
) -> None:
    """Draw the primary beam circle and the model footprint rectangle on a Matplotlib axes."""
    from matplotlib.patches import Circle, Rectangle

    half_fov = fov_deg / 2.0
    # primary beam circle
    ax.add_patch(
        Circle(
            (center_ra, center_dec),
            half_fov,
            fill=False,
            color="tab:red",
            linestyle="--",
            linewidth=1.4,
            label="primary beam",
        )
    )

    # model footprint rectangle; degrade to a tiny marker when WCS is unavailable
    if model_center is not None and model_extent is not None:
        half_w = model_extent["width_deg"] / 2.0
        half_h = model_extent["height_deg"] / 2.0
        ax.add_patch(
            Rectangle(
                (model_center.ra.deg - half_w, model_center.dec.deg - half_h),
                2.0 * half_w,
                2.0 * half_h,
                fill=False,
                color="tab:cyan",
                linestyle="-",
                linewidth=1.2,
                label="model footprint",
            )
        )
    elif model_center is not None:
        ax.scatter(
            [model_center.ra.deg],
            [model_center.dec.deg],
            marker="x",
            color="tab:cyan",
            s=80,
            label="model centre (no extent)",
        )


def _preview_cmap():
    """Return the same rainforest sub-cmap used by the single-panel previews."""
    import cmasher as cmr

    return cmr.get_sub_cmap("cmr.rainforest", 0.30, 0.85)
