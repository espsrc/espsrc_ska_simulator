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
from skasim.loaders.image_models import (
    image_model_center,
    import_casa_tasks,
    inject_image_models,
    merge_model_data_into_data,
    prepare_casa_taylor_terms,
    run_casa_exportfits,
    run_casa_ft,
    run_casa_importfits,
    run_casa_set_spectral_coordinate,
    validate_casa_taylor_terms,
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


def test_weblog_renders_casa_taylor_term_preview(tmp_path, monkeypatch):
    """CASA Taylor-term previews are exported to FITS before rendering."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")
    tt0 = tmp_path / "model.tt0"
    tt1 = tmp_path / "model.tt1"
    tt0.mkdir()
    tt1.mkdir()
    (tt0 / "table.dat").write_text("table", encoding="utf-8")
    (tt1 / "table.dat").write_text("table", encoding="utf-8")
    cfg = SimConfig(
        models=[
            {
                "type": "casa_taylor_terms",
                "tt0": str(tt0),
                "tt1": str(tt1),
                "reference_frequency_hz": 1.5e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    calls = []

    def fake_export(work_dir, imagename, fitsimage):
        calls.append((work_dir, imagename, fitsimage))
        stokes_i, _ = _write_model_pair(tmp_path, shape=(32, 32))
        fitsimage.write_bytes(stokes_i.read_bytes())

    monkeypatch.setattr("skasim.loaders.image_models.run_casa_exportfits", fake_export)

    write_image_model_previews(
        ctx,
        SkyCoord(10.0 * u.deg, 2.0 * u.deg),
        0.05 * u.deg,
    )
    html = render_weblog(ctx.manifest, ctx.work_dir)

    assert calls[0][1] == tt0.resolve()
    assert "FITS Model" in html
    assert "casa_taylor_terms" in html
    assert (ctx.work_dir / "run_fits_model.png").exists()


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
    ft_calls = [c for c in calls if "visibility_path" in c]
    assert ft_calls[0].get("incremental") is False
    assert ft_calls[1].get("incremental") is True
    assert any(c.get("merged") == "fake.ms" for c in calls)


def test_merge_model_data_into_data(monkeypatch, tmp_path):
    """MODEL_DATA is accumulated into DATA."""
    calls = []

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
            return {
                "DATA": np.ones((4, 1, 1)),
                "MODEL_DATA": np.ones((4, 1, 1)) * 0.5,
            }[name]

        def putcol(self, name, value):
            calls.append((name, value))

    # inject fake casacore.tables.table before the function imports it
    fake_tables = ModuleType("casacore.tables")
    fake_tables.table = FakeTable
    monkeypatch.setitem(sys.modules, "casacore.tables", fake_tables)

    merge_model_data_into_data(tmp_path / "visibilities.MS")

    assert len(calls) == 1
    assert calls[0][0] == "DATA"
    assert np.allclose(calls[0][1], 1.5)


def test_casa_taylor_terms_validate_and_prepare_existing_images(tmp_path, monkeypatch):
    """Existing CASA Taylor-term image directories are copied and aligned."""
    tt0 = tmp_path / "model.tt0"
    tt1 = tmp_path / "model.tt1"
    tt0.mkdir()
    tt1.mkdir()
    (tt0 / "table.dat").write_text("table", encoding="utf-8")
    (tt1 / "table.dat").write_text("table", encoding="utf-8")
    cfg = SimConfig(
        models=[
            {
                "type": "casa_taylor_terms",
                "tt0": str(tt0),
                "tt1": str(tt1),
                "reference_frequency_hz": 1.5e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    calls = []

    # fake casacore table to satisfy require_casacore + pixel reads
    fake_map_data = np.ones((4, 4), dtype=np.float32)

    class FakeCasacoreTable:
        def __init__(self, path, readonly=True, ack=False):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def getcol(self, name):
            return fake_map_data

        def putcol(self, name, value):
            pass

        def getcolkeywords(self, name):
            return {"dimnames": ["frequency"]}

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casacore",
        lambda: FakeCasacoreTable,
    )
    monkeypatch.setattr(
        "skasim.loaders.image_models.run_casa_set_spectral_coordinate",
        lambda work_dir, image_paths, frequency_hz: calls.append(
            (work_dir, image_paths, frequency_hz)
        ),
    )

    report = validate_casa_taylor_terms(cfg.models[0])
    product = prepare_casa_taylor_terms(ctx, cfg.models[0], 0)

    assert report["nterms"] == 2
    assert product.model_paths == [
        ctx.work_dir / "model_entry_01_casa_taylor.tt0.image",
        ctx.work_dir / "model_entry_01_casa_taylor.tt1.image",
    ]
    assert product.nterms == 2
    assert product.reffreq == "700000000.0Hz"
    assert (product.model_paths[0] / "table.dat").exists()
    # adjust_spectral_reference calls _set_crval4_via_script for each term
    assert len(calls) == 2
    assert calls[0] == (ctx.work_dir, [product.model_paths[0]], 700_000_000.0)
    assert calls[1] == (ctx.work_dir, [product.model_paths[1]], 700_000_000.0)


def test_inject_image_models_runs_casa_taylor_terms(tmp_path, monkeypatch):
    """CASA Taylor-term entries are predicted with ft and merged into DATA."""
    tt0 = tmp_path / "model.tt0"
    tt1 = tmp_path / "model.tt1"
    tt0.mkdir()
    tt1.mkdir()
    (tt0 / "table.dat").write_text("table", encoding="utf-8")
    (tt1 / "table.dat").write_text("table", encoding="utf-8")
    cfg = SimConfig(
        models=[
            {
                "type": "casa_taylor_terms",
                "tt0": str(tt0),
                "tt1": str(tt1),
                "reference_frequency_hz": 1.5e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    calls = []

    # fake casacore table to satisfy require_casacore + pixel reads
    fake_map_data = np.ones((4, 4), dtype=np.float32)

    class FakeCasacoreTable:
        def __init__(self, path, readonly=True, ack=False):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def getcol(self, name):
            return fake_map_data

        def putcol(self, name, value):
            pass

        def getcolkeywords(self, name):
            return {"dimnames": ["frequency"]}

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casacore",
        lambda: FakeCasacoreTable,
    )
    monkeypatch.setattr(
        "skasim.loaders.image_models.run_casa_set_spectral_coordinate",
        lambda work_dir, image_paths, frequency_hz: calls.append(
            {"spectral_copy": image_paths, "frequency_hz": frequency_hz}
        ),
    )
    monkeypatch.setattr(
        "skasim.loaders.image_models.run_casa_ft",
        lambda **kwargs: calls.append(kwargs),
    )
    monkeypatch.setattr(
        "skasim.loaders.image_models.merge_model_data_into_data",
        lambda visibility_path: calls.append({"merged": visibility_path}),
    )

    inject_image_models(ctx, tmp_path / "visibilities.MS")

    copied_paths = [
        ctx.work_dir / "model_entry_01_casa_taylor.tt0.image",
        ctx.work_dir / "model_entry_01_casa_taylor.tt1.image",
    ]
    # prepare_casa_taylor_terms now adjusts spectral reference inline:
    #   calls[0]: adjust_spectral_reference for tt0 (via run_casa_set_spectral_coordinate)
    #   calls[1]: adjust_spectral_reference for tt1 (via run_casa_set_spectral_coordinate)
    #   calls[2]: run_casa_ft
    #   calls[3]: merge_model_data_into_data
    assert calls[0]["frequency_hz"] == 700_000_000.0
    assert calls[1]["frequency_hz"] == 700_000_000.0
    assert calls[2]["model_paths"] == copied_paths
    assert calls[2]["nterms"] == 2
    assert calls[2]["reffreq"] == "700000000.0Hz"
    assert calls[3]["merged"] == tmp_path / "visibilities.MS"
    assert ctx.manifest.milestones[-1].name == "image_injection_completed"


def test_run_casa_importfits_uses_batch_fallback(tmp_path, monkeypatch):
    """CASA importfits can run through a local CASA executable when casatasks is absent."""
    calls = []

    class Result:
        returncode = 0
        stdout = "ok"

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        # create the side-effects that the real CASA importfits would produce
        for _fitsimage, imagename in [
            (tmp_path / "model.tt0.fits", tmp_path / "model.tt0.image"),
        ]:
            table_dir = (
                Path(imagename) if not isinstance(imagename, Path) else imagename
            )
            table_dir.mkdir(parents=True, exist_ok=True)
            (table_dir / "table.dat").write_text("fake", encoding="utf-8")
        return Result()

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casa_executable",
        lambda: Path("/opt/casa/bin/casa"),
    )
    monkeypatch.setattr("skasim.loaders.image_models.subprocess.run", fake_run)

    run_casa_importfits(
        tmp_path,
        [(tmp_path / "model.tt0.fits", tmp_path / "model.tt0.image")],
    )

    script = (tmp_path / "skasim_casa_importfits.py").read_text(encoding="utf-8")
    assert "from casatasks import importfits" in script
    assert "importfits(" in script
    assert "model.tt0.fits" in script
    assert calls[0][0][:5] == [
        "/opt/casa/bin/casa",
        "--nologger",
        "--nogui",
        "--log2term",
        "-c",
    ]
    assert calls[0][1]["cwd"] == str(tmp_path)


def test_run_casa_ft_uses_batch_fallback(tmp_path, monkeypatch):
    """CASA ft falls back to a batch CASA process when casatasks is absent."""
    import sys
    from types import ModuleType

    calls = []

    class Result:
        returncode = 0
        stdout = "ok"

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return Result()

    # mock casacore.tables.table for the post-batch verification
    class FakeTable:
        def __init__(self, path, **kwargs):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def colnames(self):
            return ["DATA", "MODEL_DATA"]

    fake_tables = ModuleType("casacore.tables")
    fake_tables.table = FakeTable
    monkeypatch.setitem(sys.modules, "casacore.tables", fake_tables)

    monkeypatch.setattr("skasim.loaders.image_models.import_casa_tasks", lambda: None)
    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casa_executable",
        lambda: Path("/opt/casa/bin/casa"),
    )
    monkeypatch.setattr("skasim.loaders.image_models.subprocess.run", fake_run)

    run_casa_ft(
        visibility_path=tmp_path / "visibilities.MS",
        model_paths=[tmp_path / "model.tt0.image", tmp_path / "model.tt1.image"],
        nterms=2,
        reffreq="1400000000.0Hz",
        incremental=True,
    )

    script = (tmp_path / "skasim_casa_ft.py").read_text(encoding="utf-8")
    assert "from casatasks import ft" in script
    assert "model.tt0.image" in script
    assert "model.tt1.image" in script
    assert "incremental=True" in script
    assert calls[0][0][-1] == str(tmp_path / "skasim_casa_ft.py")


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
