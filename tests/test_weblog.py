"""Weblog rendering behavior."""

from datetime import datetime

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
