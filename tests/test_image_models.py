"""Image-model validation and reporting behavior."""

from datetime import datetime, timezone
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from skasim.config import SimConfig
from skasim.loaders.image_models import (
    image_model_center,
    inject_image_models,
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
    fits.writeto(
        stokes_i, np.ones(shape, dtype=np.float32), stokes_header, overwrite=True
    )
    fits.writeto(
        alpha,
        np.zeros(alpha_shape or shape, dtype=np.float32),
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
    from pathlib import Path

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

    monkeypatch.setattr(
        "skasim.loaders.image_models.prepare_continuum_i_alpha_for_casa", fake_prepare
    )
    monkeypatch.setattr("skasim.loaders.image_models.run_casa_ft", fake_run_casa_ft)
    monkeypatch.setattr(
        "skasim.loaders.image_models.merge_model_data_into_data",
        lambda visibility_path: calls.append({"merged": str(visibility_path)}),
    )

    inject_image_models(ctx, visibility_path=Path("fake.ms"))

    # two entries -> two ft calls + one merge
    assert len(calls) == 3
    ft_calls = [c for c in calls if "model" in c]
    assert ft_calls[0].get("incremental") is False
    assert ft_calls[1].get("incremental") is True
    assert any(c.get("merged") == "fake.ms" for c in calls)


def test_merge_model_data_into_data():
    """MODEL_DATA is accumulated into DATA."""
    calls = []

    class FakeTable:
        def __init__(self, *args, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def colnames(self):
            return ["DATA", "MODEL_DATA"]

        def getcol(self, name):
            if name == "DATA":
                return np.ones((4, 1, 1))
            if name == "MODEL_DATA":
                return np.ones((4, 1, 1)) * 0.5
            return None

        def putcol(self, name, value):
            calls.append((name, value))

    import skasim.loaders.image_models as im_mod

    orig_table = getattr(im_mod, "table", None)
    im_mod.table = lambda *a, **k: FakeTable(*a, **k)
    try:
        merge_model_data_into_data(Path("fake.ms"))
    finally:
        if orig_table is not None:
            im_mod.table = orig_table
        else:
            if hasattr(im_mod, "table"):
                delattr(im_mod, "table")

    assert len(calls) == 1
    assert calls[0][0] == "DATA"
    assert np.allclose(calls[0][1], 1.5)
