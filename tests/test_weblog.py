"""Weblog rendering behavior."""

from datetime import datetime, timedelta, timezone

import numpy as np
from astropy.io import fits

from skasim.config import SimConfig
from skasim.manifest import RunManifest
from skasim.weblog import render_weblog


def test_weblog_renders_structured_outputs(tmp_path):
    """The weblog displays structured output kinds and paths."""
    manifest = RunManifest(
        run_id="example",
        started_at=datetime(2026, 5, 22, 17, 30, 0),
        config=SimConfig(),
    )
    manifest.add_output("visibility", "visibilities.MS")
    manifest.add_output("log", "example.log")

    html = render_weblog(manifest, tmp_path)

    assert "visibility" in html
    assert "visibilities.MS" in html
    assert "log" in html
    assert "example.log" in html


def test_weblog_skips_missing_image_outputs(tmp_path):
    """Failed runs can render even if an image output record points to a missing file."""
    manifest = RunManifest(
        run_id="missing-image",
        started_at=datetime(2026, 5, 22, 17, 30, 0),
        config=SimConfig(),
    )
    manifest.add_output("image_product", "missing.png", imager="wsclean", role="image")
    manifest.mark_failed("imaging failed before writing image")

    html = render_weblog(manifest, tmp_path)

    assert "imaging failed before writing image" in html
    assert "missing.png" in html


def test_weblog_uses_pipeline_milestone_names_for_durations(tmp_path):
    """Weblog duration cards are derived from actual pipeline milestone names."""
    started = datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc)
    manifest = RunManifest(
        run_id="durations",
        started_at=started,
        config=SimConfig(),
    )
    manifest.add_milestone("simulation_started", "started")
    manifest.milestones[-1].timestamp_utc = started
    manifest.add_milestone("simulation_completed", "completed")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=75)
    manifest.add_milestone("imaging_started", "started")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=80)
    manifest.add_milestone("imaging_completed", "completed")
    manifest.milestones[-1].timestamp_utc = started + timedelta(seconds=90)

    html = render_weblog(manifest, tmp_path)

    assert "simulation 1m 15s" in html
    assert "imaging 10.0s" in html


def test_weblog_groups_wsclean_mfs_products_and_stats(tmp_path):
    """MFS products render as model/clean/residual with peak and RMS stats."""
    manifest = RunManifest(
        run_id="science-products",
        started_at=datetime(2026, 5, 22, 17, 30, 0, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    prefix = "science-products_wsclean"
    for role in ("model", "image", "residual", "psf"):
        stem = f"{prefix}-MFS-{role if role != 'image' else 'image'}"
        fits_path = tmp_path / f"{stem}.fits"
        png_path = tmp_path / f"{stem}.png"
        data = np.array([[1.0e-3, 2.0e-3], [3.0e-3, 4.0e-3]])
        if role == "residual":
            data = np.array([[0.0, 1.0e-6], [-1.0e-6, 0.0]])
        fits.writeto(fits_path, data, overwrite=True)
        png_path.write_bytes(
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
            b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02"
            b"\x00\x00\x00\x90wS\xde\x00\x00\x00\x00IEND\xaeB`\x82"
        )
        preview_role = role
        manifest.add_output(
            "image_product",
            fits_path.name,
            image_product_id=prefix,
            imager="wsclean",
            role=preview_role,
        )
        manifest.add_output(
            "image_product",
            png_path.name,
            image_product_id=prefix,
            imager="wsclean",
            role=f"{preview_role}_preview",
        )

    html = render_weblog(manifest, tmp_path)

    assert "Science Products" in html
    assert "Model" in html
    assert "Clean" in html
    assert "Residual" in html
    assert "View PSF" in html
    assert "Peak:" in html
    assert "4.000 mJy/beam" in html
    assert "RMS:" in html
