"""weblog.py — render weblog.html from a RunManifest via Jinja2."""

from __future__ import annotations

import base64
import json
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
from astropy.io import fits
from jinja2 import Environment, PackageLoader, select_autoescape

from .manifest import RunManifest

KNOWN_ANTENNA_COUNTS = {
    "MEERKAT": 64,
}


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
        if o.role == "uv_coverage":
            continue
        if "sky_model" in output_path.lower() or o.role in (
            "sky_model",
            "sky_model_fov",
        ):
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
        tag = (output.metadata or {}).get("tag", "default")
        product = products.setdefault(
            product_id,
            {
                "product_id": product_id,
                "tag": tag,
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
            "beam": _read_fits_beam(fpath),
        }

    for product in products.values():
        clean = product["clean"]
        residual = product["residual"]
        if clean and residual:
            product["stats"] = _calculate_image_stats(
                work_dir / clean["fits"],
                work_dir / residual["fits"],
            )
        product["beam"] = _science_product_beam(product)

    return [
        product
        for product in products.values()
        if any(product.get(k) for k in ("model", "clean", "residual", "dirty"))
    ]


def _science_product_beam(product: dict) -> dict | None:
    """Return one representative beam for a science product."""
    for key in ("clean", "dirty", "residual", "model", "psf"):
        item = product.get(key)
        if item and item.get("beam"):
            return item["beam"]
    return None


def _read_fits_beam(path: Path) -> dict | None:
    """Return FITS beam metadata in readable units when present."""
    if not path.exists():
        return None
    try:
        with fits.open(path) as hdul:
            header = hdul[0].header
            bmaj = header.get("BMAJ")
            bmin = header.get("BMIN")
            bpa = header.get("BPA")
    except Exception:
        return None
    if bmaj is None and bmin is None and bpa is None:
        return None
    return {
        "bmaj": _format_angle_deg(bmaj) if bmaj is not None else None,
        "bmin": _format_angle_deg(bmin) if bmin is not None else None,
        "bpa": f"{float(bpa):.2f} deg" if bpa is not None else None,
    }


def _format_angle_deg(value_deg: float) -> str:
    """Format an angular size provided in degrees."""
    value = float(value_deg)
    abs_value = abs(value)
    if abs_value >= 1.0:
        return f"{value:.4f} deg"
    if abs_value >= 1.0 / 60.0:
        return f"{value * 60.0:.3f} arcmin"
    return f"{value * 3600.0:.3f} arcsec"


def _format_float(value: float, digits: int = 3) -> str:
    """Format a number without distracting trailing zeroes."""
    return f"{float(value):.{digits}f}".rstrip("0").rstrip(".")


def _observation_summary(manifest: RunManifest) -> dict:
    """Build display rows for the resolved spectral and time setup."""
    obs = manifest.config.observation
    details = _find_milestone_details(manifest, "observation_configured")
    centre = float(obs.frequency_mhz)
    n_channels = int(obs.n_channels or 1)
    channel_width = float(obs.channel_width_mhz or 0.0)
    total_bandwidth = float(obs.bandwidth_mhz or n_channels * channel_width)
    min_freq = float(details.get("min_frequency_mhz", centre - total_bandwidth / 2.0))
    max_freq = float(details.get("max_frequency_mhz", centre + total_bandwidth / 2.0))
    n_timesteps = details.get("n_timesteps")
    integration_time = None
    if n_timesteps:
        integration_time = float(obs.observation_time_s) / float(n_timesteps)
    return {
        "frequency_min_mhz": _format_float(min_freq),
        "frequency_max_mhz": _format_float(max_freq),
        "central_frequency_mhz": _format_float(centre),
        "n_channels": n_channels,
        "channel_bandwidth_mhz": _format_float(channel_width),
        "total_bandwidth_mhz": _format_float(total_bandwidth),
        "observation_time_s": _format_float(obs.observation_time_s),
        "integration_time_s": (
            _format_float(integration_time) if integration_time is not None else None
        ),
        "n_timesteps": n_timesteps,
    }


def _imaging_summary_for_block(
    img_config,
    fov_deg: float | None,
    beam: dict | None,
) -> dict:
    """Build display rows for imaging geometry and pixelization for a single ImgConfig block."""
    pixel_size_arcsec = None
    if fov_deg is not None:
        pixel_size_arcsec = fov_deg * 3600.0 / img_config.pixels
    return {
        "imager": img_config.imager,
        "pixels": img_config.pixels,
        "image_size": f"{img_config.pixels} x {img_config.pixels}",
        "fov_deg": _format_float(fov_deg, 4) if fov_deg is not None else None,
        "pixel_size_arcsec": (
            _format_float(pixel_size_arcsec, 4)
            if pixel_size_arcsec is not None
            else None
        ),
        "beam": _format_beam(beam),
    }


def _format_beam(beam: dict | None) -> str | None:
    """Format one beam dictionary for compact display."""
    if not beam:
        return None
    parts = []
    if beam.get("bmaj"):
        parts.append(f"BMAJ {beam['bmaj']}")
    if beam.get("bmin"):
        parts.append(f"BMIN {beam['bmin']}")
    if beam.get("bpa"):
        parts.append(f"BPA {beam['bpa']}")
    return " · ".join(parts) if parts else None


def _imager_parameter_rows_for_block(
    img_config,
    manifest: RunManifest,
    fov_deg: float | None,
    center: tuple[float, float] | None,
) -> list[tuple[str, str]]:
    """Return the effective OSKAR or WSClean parameters for a single ImgConfig block."""
    config = manifest.config
    n_channels = int(config.observation.n_channels or 1)
    channels_out = min(n_channels, 8)
    output_prefix = _first_image_product_id(manifest, img_config.imager)
    visibility_input = _first_output_path(manifest, "visibility")
    pixel_size_arcsec = None
    if fov_deg is not None:
        pixel_size_arcsec = fov_deg * 3600.0 / img_config.pixels

    if img_config.imager == "wsclean":
        rows = [
            ("Imager", "WSClean"),
            ("Command", img_config.wsclean_command),
            ("Weighting", f"Briggs robust {_format_float(img_config.robust)}"),
            ("Multiscale", "enabled"),
            ("Image size", f"{img_config.pixels} x {img_config.pixels} pixels"),
            (
                "Pixel scale",
                f"{_format_float(pixel_size_arcsec, 4)} arcsec"
                if pixel_size_arcsec is not None
                else "unknown",
            ),
            ("Clean iterations", str(img_config.clean_iterations)),
            ("Major-cycle gain", "0.8"),
            ("Auto threshold", "0.3"),
            ("Auto mask", "3"),
            ("Channels out", str(channels_out)),
            ("Join channels", "enabled"),
            ("Local RMS", "enabled"),
        ]
        if output_prefix is not None:
            rows.append(("Output prefix", output_prefix))
        if visibility_input is not None:
            rows.append(("Visibility input", visibility_input))
        return rows

    rows = [
        ("Imager", "OSKAR dirty"),
        ("Image size", f"{img_config.pixels} x {img_config.pixels} pixels"),
        (
            "Pixel scale",
            f"{_format_float(pixel_size_arcsec, 4)} arcsec"
            if pixel_size_arcsec is not None
            else "unknown",
        ),
        ("Combine frequencies", "enabled"),
    ]
    if center is not None:
        rows.append(
            (
                "Phase centre",
                f"RA {_format_float(center[0], 6)} deg, "
                f"Dec {_format_float(center[1], 6)} deg",
            )
        )
    if visibility_input is not None:
        rows.append(("Visibility input", visibility_input))
    return rows


def _build_imaging_tabs(
    manifest: RunManifest,
    science_products: list[dict],
    fov_deg_default: float | None,
    center: tuple[float, float] | None,
) -> list[dict]:
    """Build one tab dict per imaging block with config summary, parameters, and products."""
    config = manifest.config
    if not config.imaging:
        return []

    tabs: list[dict] = []

    for img_config in config.imaging:
        tag = img_config.tag

        # Compute FoV for this block
        fov_deg = img_config.fov_deg
        if fov_deg is None:
            fov_deg = fov_deg_default

        # Find beam from this tag's products
        tag_products = [p for p in science_products if p.get("tag") == tag]
        beam = None
        for p in tag_products:
            if p.get("beam"):
                beam = p["beam"]
                break

        tabs.append(
            {
                "tag": tag,
                "imaging_summary": _imaging_summary_for_block(
                    img_config, fov_deg, beam
                ),
                "imager_parameter_rows": _imager_parameter_rows_for_block(
                    img_config,
                    manifest,
                    fov_deg,
                    center,
                ),
                "products": tag_products,
            }
        )

    return tabs


def _first_image_product_id(manifest: RunManifest, imager: str) -> str | None:
    """Return the first image product ID for an imager."""
    for output in manifest.outputs:
        if output.kind == "image_product" and output.imager == imager:
            return output.image_product_id
    return None


def _first_output_path(manifest: RunManifest, kind: str) -> str | None:
    """Return the first output path for a manifest output kind."""
    for output in manifest.outputs:
        if output.kind == kind:
            return output.path
    return None


def _find_milestone_details(manifest: RunManifest, name: str) -> dict:
    """Return details for the first milestone with a matching name."""
    for milestone in manifest.milestones:
        if milestone.name == name:
            return milestone.details or {}
    return {}


def _antenna_count(manifest: RunManifest) -> int | None:
    """Return antenna/station count from the manifest or known telescope metadata."""
    details = _find_milestone_details(manifest, "telescope_built")
    count = details.get("n_stations")
    if count is not None:
        return int(count)
    return KNOWN_ANTENNA_COUNTS.get(manifest.config.telescope.upper())


def _software_versions() -> list[tuple[str, str]]:
    """Return concise runtime software versions for weblog reproducibility."""
    return [
        ("skasim", _package_version("skasim")),
        ("Karabo", _package_version("karabo-pipeline")),
        ("OSKAR", _conda_package_version("oskarpy") or _package_version("oskarpy")),
        ("WSClean", _conda_package_version("wsclean") or _package_version("wsclean")),
        ("Python", sys.version.split()[0]),
    ]


def _package_version(package_name: str) -> str:
    """Return Python package metadata version or an explicit unknown marker."""
    try:
        return version(package_name)
    except PackageNotFoundError:
        return "unknown"


def _conda_package_version(package_name: str) -> str | None:
    """Return Conda package version from the active prefix when available."""
    conda_meta = Path(sys.prefix) / "conda-meta"
    if not conda_meta.exists():
        return None
    for path in sorted(conda_meta.glob(f"{package_name}-*.json")):
        try:
            metadata = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if metadata.get("name") == package_name and metadata.get("version"):
            return str(metadata["version"])
    return None


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

    # aggregate imaging duration across all tags
    imaging_starts = [
        m.timestamp_utc
        for m in manifest.milestones
        if m.name.endswith("_started") and m.name.startswith("imaging_")
    ]
    imaging_ends = [
        m.timestamp_utc
        for m in manifest.milestones
        if m.name.endswith("_completed") and m.name.startswith("imaging_")
    ]
    imaging_duration = None
    if imaging_starts and imaging_ends:
        imaging_duration = (max(imaging_ends) - min(imaging_starts)).total_seconds()

    # Derived telescope properties
    from .utils import get_diameter

    dish_diameter = None
    derived_fov = None
    fov_deg_value = (
        manifest.config.imaging[0].fov_deg if manifest.config.imaging else None
    )
    telescope_name = manifest.config.telescope
    freq_mhz = manifest.config.observation.frequency_mhz

    try:
        diam = get_diameter(telescope_name)
        dish_diameter = f"{diam.value:.1f} {diam.unit}"
        if freq_mhz:
            wavelength_m = 299.792458 / freq_mhz
            fov_rad = 1.25 * wavelength_m / diam.value
            fov_deg_val = fov_rad * 180.0 / 3.141592653589793
            if fov_deg_value is None:
                fov_deg_value = fov_deg_val
            derived_fov = f"{fov_deg_val:.2f}°"
    except Exception:
        pass

    n_stations = _antenna_count(manifest)

    # Resolve telescope plot if present
    telescope_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        if "telescope" in output_path.lower() or o.role == "telescope":
            fpath = work_dir / output_path
            if fpath.exists():
                telescope_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath),
                }
            break

    # Resolve UV coverage plot if present
    uv_coverage_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role == "uv_coverage" or "uvcoverage" in output_path.lower():
            fpath = work_dir / output_path
            if fpath.exists() and output_path.endswith((".png", ".jpg", ".jpeg")):
                uv_coverage_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath),
                }
                break

    # Resolve sky model plot if present
    sky_model_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role == "sky_model" or (
            "sky_model" in output_path.lower() and "_fov" not in output_path.lower()
        ):
            fpath = work_dir / output_path
            if fpath.exists():
                sky_model_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath),
                }
            break

    # Resolve sky model fov plot if present
    sky_model_fov_plot = None
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role == "sky_model_fov" or (
            "sky_model" in output_path.lower() and "_fov" in output_path.lower()
        ):
            fpath = work_dir / output_path
            if fpath.exists():
                sky_model_fov_plot = {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath),
                }
            break

    # Resolve FITS image-model preview if present
    fits_model_plots = []
    for o in manifest.outputs:
        output_path = o.path if hasattr(o, "path") else o
        role = getattr(o, "role", None)
        if role != "fits_model":
            continue
        fpath = work_dir / output_path
        if fpath.exists() and output_path.endswith((".png", ".jpg", ".jpeg")):
            fits_model_plots.append(
                {
                    "path": output_path,
                    "data": _file_to_base64_data_uri(fpath),
                    "model_type": o.metadata.get("model_type")
                    if hasattr(o, "metadata")
                    else None,
                }
            )

    observation_details = _find_milestone_details(manifest, "observation_configured")
    center = None
    if (
        "phase_center_ra_deg" in observation_details
        and "phase_center_dec_deg" in observation_details
    ):
        center = (
            float(observation_details["phase_center_ra_deg"]),
            float(observation_details["phase_center_dec_deg"]),
        )

    science_products = _find_science_products(manifest, work_dir)
    imaging_tabs = _build_imaging_tabs(
        manifest,
        science_products,
        fov_deg_value,
        center,
    )

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
        imaging_tabs=imaging_tabs,
        telescope_plot=telescope_plot,
        uv_coverage_plot=uv_coverage_plot,
        sky_model_plot=sky_model_plot,
        sky_model_fov_plot=sky_model_fov_plot,
        fits_model_plots=fits_model_plots,
        observation_summary=_observation_summary(manifest),
        software_versions=_software_versions(),
        dish_diameter=dish_diameter,
        derived_fov=derived_fov,
        n_stations=n_stations,
    )

    weblog_path = work_dir / "weblog.html"
    weblog_path.write_text(html, encoding="utf-8")
    return html
