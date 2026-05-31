"""Image-model validation and reporting behavior."""

import sys
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from skasim.config import SimConfig
from skasim.image_models import (
    inject_image_models,
    image_model_center,
    merge_model_data_into_data,
    validate_continuum_i_alpha,
    write_image_model_previews,
)
from skasim.manifest import RunManifest, create_run_context
from skasim.weblog import render_weblog


def _write_model_pair(tmp_path, shape=(16, 16), alpha_shape=None):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10.0, 2.0]
    wcs.wcs.crpix = [shape[1] / 2.0, shape[0] / 2.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    stokes_header = wcs.to_header()
    stokes_header["BUNIT"] = "Jy/pixel"
    alpha_header = wcs.to_header()
    alpha_header["BUNIT"] = "1"

    stokes_i = tmp_path / "stokes_i.fits"
    alpha = tmp_path / "alpha.fits"
    fits.writeto(stokes_i, np.ones(shape), stokes_header, overwrite=True)
    fits.writeto(
        alpha,
        np.zeros(alpha_shape or shape),
        alpha_header,
        overwrite=True,
    )
    return stokes_i, alpha


def test_validate_continuum_i_alpha_accepts_matching_fits(tmp_path):
    """Continuum I+alpha validation reports model metadata."""
    stokes_i, alpha = _write_model_pair(tmp_path)
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ]
    )

    report = validate_continuum_i_alpha(cfg.models[0])

    assert report["spatial_shape"] == [16, 16]
    assert report["unit"] == "Jy/pixel"
    assert report["reference_frequency_hz"] == pytest.approx(1.4e9)


def test_validate_continuum_i_alpha_rejects_shape_mismatch(tmp_path):
    """Stokes I and alpha maps must share a spatial grid."""
    stokes_i, alpha = _write_model_pair(tmp_path, alpha_shape=(8, 8))
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ]
    )

    with pytest.raises(ValueError, match="matching spatial dimensions"):
        validate_continuum_i_alpha(cfg.models[0])


def test_validate_continuum_i_alpha_accepts_real_model_fixture(tmp_path):
    """Real WSClean model FITS headers are accepted as Jy/pixel image models."""
    fixture = (
        Path(__file__).resolve().parents[1]
        / "models_for_testing"
        / "CY4223_L_004_20180322_avg_target_3C320_vla_loop_01-0000-model.fits"
    )
    if not fixture.exists():
        pytest.skip("models_for_testing fixture directory is not present")
    alpha = tmp_path / "alpha.fits"
    with fits.open(fixture) as hdul:
        header = hdul[0].header.copy()
        header["BUNIT"] = "1"
        fits.writeto(alpha, np.zeros_like(hdul[0].data), header, overwrite=True)
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(fixture),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ]
    )

    report = validate_continuum_i_alpha(cfg.models[0])

    assert report["spatial_shape"] == [1600, 1600]
    assert report["unit"] == "JY/PIXEL"


def test_image_model_center_uses_first_model_wcs(tmp_path):
    """Image-only runs can derive the phase centre from FITS WCS."""
    stokes_i, alpha = _write_model_pair(tmp_path)
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ]
    )

    center = image_model_center(cfg.models)

    assert center is not None
    assert center.ra.deg == pytest.approx(10.0, abs=1e-3)
    assert center.dec.deg == pytest.approx(2.0, abs=1e-3)


def test_weblog_renders_fits_model_preview(tmp_path):
    """FITS model previews appear in the weblog Sky Model section."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")
    stokes_i, alpha = _write_model_pair(tmp_path, shape=(32, 32))
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)

    write_image_model_previews(
        ctx,
        SkyCoord(10.0 * u.deg, 2.0 * u.deg),
        0.05 * u.deg,
    )
    html = render_weblog(ctx.manifest, ctx.work_dir)

    assert "FITS Model" in html
    assert "continuum_i_alpha" in html
    assert (ctx.work_dir / "run_fits_model.png").exists()


def test_weblog_renders_existing_fits_model_output(tmp_path):
    """The weblog can render a pre-recorded FITS model preview output."""
    manifest = RunManifest(
        run_id="preview",
        started_at=datetime(2026, 5, 31, tzinfo=timezone.utc),
        config=SimConfig(),
    )
    png_path = tmp_path / "preview_fits_model.png"
    png_path.write_bytes(
        b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
        b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02"
        b"\x00\x00\x00\x90wS\xde\x00\x00\x00\x00IEND\xaeB`\x82"
    )
    manifest.add_output(
        "plot",
        png_path.name,
        role="fits_model",
        metadata={"model_type": "continuum_i_alpha"},
    )

    html = render_weblog(manifest, tmp_path)

    assert "FITS Model" in html
    assert "continuum_i_alpha" in html


def test_inject_image_models_runs_continuum_entries_in_order(tmp_path, monkeypatch):
    """Image injection validates entries, predicts with CASA ft, and merges once."""
    stokes_i, alpha = _write_model_pair(tmp_path)
    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            },
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes_i),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            },
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    calls = []

    class Product:
        nterms = 2
        reffreq = "1400000000.0Hz"

        def __init__(self, index):
            self.model_paths = [ctx.work_dir / f"model-{index}.image"]
            self.intermediates = []

    def fake_prepare(ctx_arg, entry, index):
        assert ctx_arg is ctx
        return Product(index)

    def fake_run_casa_ft(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr("skasim.image_models.prepare_continuum_i_alpha_for_casa", fake_prepare)
    monkeypatch.setattr("skasim.image_models.run_casa_ft", fake_run_casa_ft)
    monkeypatch.setattr(
        "skasim.image_models.merge_model_data_into_data",
        lambda visibility_path: calls.append({"merged": visibility_path}),
    )

    inject_image_models(ctx, tmp_path / "visibilities.MS")

    assert [call.get("incremental") for call in calls[:2]] == [False, True]
    assert calls[2]["merged"] == tmp_path / "visibilities.MS"
    assert [m.name for m in ctx.manifest.milestones].count("image_model_injected") == 2
    assert ctx.manifest.milestones[-1].name == "image_injection_completed"


def test_merge_model_data_into_data_adds_model_column(monkeypatch, tmp_path):
    """The final delivered DATA column includes image MODEL_DATA."""
    data = np.array([[[1 + 0j, 2 + 0j]]])
    model = np.array([[[3 + 0j, 4 + 0j]]])
    written = {}

    class FakeTable:
        def __init__(self, *_args, **_kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *_exc):
            return False

        def colnames(self):
            return ["DATA", "MODEL_DATA"]

        def getcol(self, name):
            return {"DATA": data, "MODEL_DATA": model}[name]

        def putcol(self, name, value):
            written[name] = value

    casacore_module = ModuleType("casacore")
    tables_module = ModuleType("casacore.tables")
    tables_module.table = FakeTable
    monkeypatch.setitem(sys.modules, "casacore", casacore_module)
    monkeypatch.setitem(sys.modules, "casacore.tables", tables_module)

    merge_model_data_into_data(tmp_path / "visibilities.MS")

    assert np.array_equal(written["DATA"], data + model)
