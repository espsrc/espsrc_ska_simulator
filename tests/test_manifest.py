"""Run manifest behavior."""

import json
from datetime import datetime

from skasim.config import SimConfig
from skasim.manifest import RunManifest, create_run_context


def test_manifest_serializes_structured_outputs():
    """Manifest outputs include kind, path, and image product metadata."""
    manifest = RunManifest(
        run_id="example",
        started_at=datetime(2026, 5, 22, 17, 30, 0),
        config=SimConfig(),
    )

    manifest.add_output(kind="visibility", path="visibilities.MS")
    manifest.add_output(
        kind="image_product",
        path="example_wsclean-MFS-image.fits",
        image_product_id="example_wsclean",
        imager="wsclean",
        role="image",
    )

    data = json.loads(manifest.model_dump_json())

    assert data["outputs"][0]["kind"] == "visibility"
    assert data["outputs"][0]["path"] == "visibilities.MS"
    assert data["outputs"][1]["kind"] == "image_product"
    assert data["outputs"][1]["imager"] == "wsclean"
    assert data["outputs"][1]["role"] == "image"
    assert "example_wsclean-MFS-image.fits" in manifest.model_dump_json()


def test_create_run_context_records_log_and_manifest_outputs(tmp_path):
    """Run contexts identify log and manifest outputs by kind."""
    config = SimConfig(output_prefix=str(tmp_path / "example"))
    ctx = create_run_context(config)

    outputs = {output.kind: output.path for output in ctx.manifest.outputs}

    assert outputs["log"].endswith(".log")
    assert outputs["manifest"] == "run_manifest.json"
