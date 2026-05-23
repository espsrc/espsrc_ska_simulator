"""Weblog rendering behavior."""

from datetime import datetime, timedelta, timezone

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

    assert "<h3>simulation</h3>" in html
    assert "1m 15s" in html
    assert "<h3>imaging</h3>" in html
    assert "10.0s" in html
