"""weblog.py — render weblog.html from a RunManifest via Jinja2."""

from __future__ import annotations

import base64
from pathlib import Path

import numpy as np
from astropy.io import fits
from jinja2 import Environment, PackageLoader, select_autoescape

from .manifest import RunManifest


def _humanize_seconds(total_s: float) -> str:
    """convert seconds to '2m 34s' or '1h 2m 3s'."""
    if total_s < 60:
        return f"{total_s:.1f}s"
    minutes, seconds = divmod(int(total_s), 60)
    if minutes < 60:
        return f"{minutes}m {seconds}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h {minutes}m {seconds}s"


def _find_image_outputs(manifest: RunManifest, work_dir: Path) -> list[dict]:
    """extract image paths, captions, and base64 data for inline embedding."""
    results: list[dict] = []
    for o in manifest.outputs:
        if o.kind == "image_product":
            continue
        output_path = o.path if hasattr(o, "path") else o
        if "telescope" in output_path.lower() or o.role == "telescope":
            continue
        if "sky_model" in output_path.lower() or o.role in ("sky_model", "sky_model_fov"):
            continue
        if not output_path.endswith((".png", ".jpg", ".jpeg")):
            continue
        fpath = work_dir / output_path
        if not fpath.exists():
            continue  # skip missing files (e.g. during failed runs)
        caption = _infer_image_caption(Path(output_path).name)
        data_b64 = _file_to_base64_data_uri(fpath)
        results.append({"path": output_path, "caption": caption, "data": data_b64})
    return results


def _find_science_products(manifest: RunManifest, work_dir: Path) -> list[dict]:
    """group MFS and dirty image products into model, clean, residual, PSF, and dirty entries."""
    products: dict[str, dict] = {}
    for output in manifest.outputs:
        if output.kind != "image_product":
            continue
        if output.imager == "wsclean" and "-MFS-" not in output.path:
            continue

        product_id = output.image_product_id or "imaging"
        product = products.setdefault(
            product_id,
            {
                "product_id": product_id,
                "imager": output.imager,
                "model": None,
                "clean": None,
                "residual": None,
                "psf": None,
                "dirty": None,
                "stats": None,
            },
        )
        role = output.role or ""
        key = None
        if role == "model":
            key = "model"
        elif role == "image":
            if output.imager == "oskar-dirty":
                key = "dirty"
            else:
                key = "clean"
        elif role == "dirty":
            key = "dirty"
        elif role == "residual":
            key = "residual"
        elif role == "psf":
            key = "psf"
        elif role.endswith("_preview"):
            continue

        if key is None:
            continue
        fpath = work_dir / output.path
        preview = fpath.with_suffix(".png")
        product[key] = {
            "fits": output.path,
            "preview": preview.name if preview.exists() else None,
            "data": _file_to_base64_data_uri(preview) if preview.exists() else None,
        }

    for product in products.values():
        clean = product["clean"]
        residual = product["residual"]
        if clean and residual:
            product["stats"] = _calculate_image_stats(
                work_dir / clean["fits"],
                work_dir / residual["fits"],
            )

    return [product for product in products.values() if any(product.get(k) for k in ("model", "clean", "residual", "dirty"))]


def _calculate_image_stats(clean_path: Path, residual_path: Path) -> dict | None:
    """calculate peak, residual RMS, and peak/RMS from FITS image data."""
    if not clean_path.exists() or not residual_path.exists():
        return None
    try:
        clean_data = _read_fits_image(clean_path)
        residual_data = _read_fits_image(residual_path)
    except Exception:
        return None
    finite_clean = clean_data[np.isfinite(clean_data)]
    finite_residual = residual_data[np.isfinite(residual_data)]
    if finite_clean.size == 0 or finite_residual.size == 0:
        return None

    peak_jy = float(np.nanmax(finite_clean))
    rms_jy = float(np.nanstd(finite_residual))
    snr = peak_jy / rms_jy if rms_jy > 0 else None
    return {
        "peak": _format_flux(peak_jy),
        "rms": _format_flux(rms_jy),
        "snr": f"{snr:.1f}" if snr is not None else None,
    }


def _read_fits_image(path: Path) -> np.ndarray:
    """read a FITS image and reduce singleton axes to a 2D array."""
    with fits.open(path) as hdul:
        data = hdul[0].data
        if data is None:
            raise ValueError(f"{path} has no image data")
        data = np.asarray(data).squeeze()
    while data.ndim > 2:
        data = data[0]
    return data


def _format_flux(val_jy: float) -> str:
    """Format flux value dynamically using the most suitable unit (Jy, mJy, uJy, nJy)."""
    val_abs = abs(val_jy)
    if val_abs >= 1.0:
        return f"{val_jy:.3f} Jy/beam"
    elif val_abs >= 1e-3:
        return f"{val_jy * 1e3:.3f} mJy/beam"
    elif val_abs >= 1e-6:
        return f"{val_jy * 1e6:.2f} μJy/beam"
    else:
        return f"{val_jy * 1e9:.2f} nJy/beam"


def _file_to_base64_data_uri(fpath: Path) -> str:
    """read a binary image and return a base64 data URI string."""
    data = fpath.read_bytes()
    b64 = base64.b64encode(data).decode("ascii")
    # infer mime type from extension
    ext = fpath.suffix.lower()
    mime = "image/png" if ext == ".png" else "image/jpeg"
    return f"data:{mime};base64,{b64}"


def _infer_image_caption(filename: str) -> str:
    """return a human-readable caption from an output filename."""
    fname = Path(filename).stem
    lower = fname.lower()
    if "dirty" in lower:
        if "oskar" in lower or "-dirty" in lower:
            return "Dirty image (OSKAR)"
        return "Dirty image"
    if "psf" in lower:
        return "Point spread function"
    if "residual" in lower:
        return "Residual"
    if "model" in lower:
        return "Component model"
    if "clean" in lower or "image" in lower:
        return "Cleaned image (WSClean)"
    if "telescope" in lower:
        return "Telescope layout"
    return fname  # fallback to filename


def render_weblog(manifest: RunManifest, work_dir: Path) -> str:
    """generate weblog.html content and write it to disk. returns the html string."""
    env = Environment(
        loader=PackageLoader("skasim", "templates"),
        autoescape=select_autoescape(["html", "xml"]),
    )
    template = env.get_template("weblog.html.j2")

    # compute summary stats
    total_elapsed = None
    if manifest.completed_at and manifest.started_at:
        total_elapsed = (manifest.completed_at - manifest.started_at).total_seconds()

    milestone_lookup = {m.name: m for m in manifest.milestones}
    simulation_duration = None
    if (
        "simulation_started" in milestone_lookup
        and "simulation_completed" in milestone_lookup
    ):
        simulation_duration = (
            milestone_lookup["simulation_completed"].timestamp_utc
            - milestone_lookup["simulation_started"].timestamp_utc
        ).total_seconds()

    imaging_duration = None
    if (
        "imaging_started" in milestone_lookup
        and "imaging_completed" in milestone_lookup
    ):
        imaging_duration = (
            milestone_lookup["imaging_completed"].timestamp_utc
            - milestone_lookup["imaging_started"].timestamp_utc
        ).total_seconds()

    # Derived telescope properties
    from .utils import get_diameter
    dish_diameter = None
    derived_fov = None
    telescope_name = manifest.config.telescope
    freq_mhz = manifest.config.observation.frequency_mhz
    
    try:
        diam = get_diameter(telescope_name)
        dish_diameter = f"{diam.value:.1f} {diam.unit}"
        if freq_mhz:
            wavelength_m = 299.792458 / freq_mhz
            fov_rad = 1.25 * wavelength_m / diam.value
            fov_deg_val = fov_rad * 180.0 / 3.141592653589793
            derived_fov = f"{fov_deg_val:.2f}°"
    except Exception:
        pass
        
    n_stations = None
    if "telescope_built" in milestone_lookup and milestone_lookup["telescope_built"].details:
        n_stations = milestone_lookup["telescope_built"].details.get("n_stations")

    # Resolve telescope plot if present
    telescope_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        if "telescope" in output_path.lower() or o.role == "telescope":
            fpath = work_dir / output_path
            if fpath.exists():
                telescope_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath)
                }
            break

    # Resolve sky model plot if present
    sky_model_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role == "sky_model" or ("sky_model" in output_path.lower() and "_fov" not in output_path.lower()):
            fpath = work_dir / output_path
            if fpath.exists():
                sky_model_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath)
                }
            break

    # Resolve sky model fov plot if present
    sky_model_fov_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role == "sky_model_fov" or ("sky_model" in output_path.lower() and "_fov" in output_path.lower()):
            fpath = work_dir / output_path
            if fpath.exists():
                sky_model_fov_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath)
                }
            break

    html = template.render(
        manifest=manifest,
        total_elapsed=_humanize_seconds(total_elapsed) if total_elapsed else None,
        simulation_duration=(
            _humanize_seconds(simulation_duration) if simulation_duration else None
        ),
        imaging_duration=(
            _humanize_seconds(imaging_duration) if imaging_duration else None
        ),
        images=_find_image_outputs(manifest, work_dir),
        science_products=_find_science_products(manifest, work_dir),
        telescope_plot=telescope_plot,
        sky_model_plot=sky_model_plot,
        sky_model_fov_plot=sky_model_fov_plot,
        dish_diameter=dish_diameter,
        derived_fov=derived_fov,
        n_stations=n_stations,
    )

    weblog_path = work_dir / "weblog.html"
    weblog_path.write_text(html, encoding="utf-8")
    return html
