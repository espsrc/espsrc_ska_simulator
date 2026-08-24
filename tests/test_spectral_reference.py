"""Tests for spectral reference adjustment in image models."""

import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from skasim.config import CasaTaylorTermsModelEntry, ContinuumIAlphaModelEntry, ObsConfig, SimConfig
from skasim.loaders.image_models import adjust_spectral_reference
from skasim.manifest import create_run_context
from skasim.runtime import CasacoreRuntimeError


# ---------------------------------------------------------------------------
# adjust_spectral_reference tests
# ---------------------------------------------------------------------------

def _make_casa_image_dir(tmp_path, name):
    """Create a minimal CASA image directory with table.dat."""
    image_dir = tmp_path / name
    image_dir.mkdir()
    (image_dir / "table.dat").write_text("fake", encoding="utf-8")
    return image_dir


def _patch_require_casacore_and_table(monkeypatch, pixel_data_map):
    """Patch require_casacore and fake_table to return pixel data from a dict.

    pixel_data_map: {image_path_str: ndarray} — maps image path to pixel data.
    Returns a dict that collects putcol calls: {image_path_str: ndarray}.
    """
    written = {}

    class FakeTable:
        def __init__(self, path, **kwargs):
            self.path = str(path)

        def __enter__(self):
            return self

        def __exit__(self, *args):
            pass

        def getcol(self, name):
            assert name == "map"
            return pixel_data_map[self.path].copy()

        def putcol(self, name, value):
            written[self.path] = value

        def getcolkeywords(self, name):
            return {"dimnames": ["frequency"]}

    def fake_require_casacore():
        return FakeTable

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casacore", fake_require_casacore,
    )
    return written


def test_adjust_spectral_reference_nterms1_crval_only(tmp_path, monkeypatch):
    """nterms=1 (alpha_map=None): only CRVAL4 is set, no pixel data modification."""
    image_dir = _make_casa_image_dir(tmp_path, "tt0.image")

    set_crval_calls = []
    monkeypatch.setattr(
        "skasim.loaders.image_models._set_crval4_via_script",
        lambda work_dir, paths, freq: set_crval_calls.append((paths, freq)),
    )

    result = adjust_spectral_reference(
        image_dir, old_ref_hz=1.5e9, new_ref_hz=1.284e9, alpha_map=None,
    )

    assert result == pytest.approx(1.284e9)
    assert len(set_crval_calls) == 1
    assert set_crval_calls[0][1] == pytest.approx(1.284e9)


def test_adjust_spectral_reference_nterms2_scales_pixel_data(tmp_path, monkeypatch):
    """nterms≥2 (alpha provided): tt0 pixel data is scaled by (ν_new/ν_old)^α."""
    image_dir = _make_casa_image_dir(tmp_path, "tt0.image")
    pixel_data = np.ones((4, 4), dtype=np.float64) * 100.0

    # mock _set_crval4_via_script (called after pixel data)
    monkeypatch.setattr(
        "skasim.loaders.image_models._set_crval4_via_script",
        lambda work_dir, paths, freq: None,
    )

    written = _patch_require_casacore_and_table(
        monkeypatch, {str(image_dir): pixel_data},
    )

    old_ref = 1.5e9
    new_ref = 1.284e9
    alpha = -0.7
    expected_factor = (new_ref / old_ref) ** alpha

    result = adjust_spectral_reference(
        image_dir, old_ref_hz=old_ref, new_ref_hz=new_ref, alpha_map=np.full_like(pixel_data, alpha),
    )

    assert result == pytest.approx(new_ref)
    assert str(image_dir) in written
    np.testing.assert_allclose(written[str(image_dir)], 100.0 * expected_factor, rtol=1e-6)


def test_adjust_spectral_reference_idempotent(tmp_path, monkeypatch):
    """When old_ref == new_ref and alpha is provided, pixel data is unchanged."""
    image_dir = _make_casa_image_dir(tmp_path, "tt0.image")
    pixel_data = np.ones((4, 4), dtype=np.float64) * 50.0

    monkeypatch.setattr(
        "skasim.loaders.image_models._set_crval4_via_script",
        lambda work_dir, paths, freq: None,
    )

    written = _patch_require_casacore_and_table(
        monkeypatch, {str(image_dir): pixel_data},
    )

    same_freq = 1.4e9
    result = adjust_spectral_reference(
        image_dir, old_ref_hz=same_freq, new_ref_hz=same_freq, alpha_map=np.full_like(pixel_data, -0.7),
    )

    assert result == pytest.approx(same_freq)
    # factor = (1.4e9 / 1.4e9) ^ -0.7 = 1.0
    np.testing.assert_allclose(written[str(image_dir)], 50.0, rtol=1e-10)


def test_adjust_spectral_reference_alpha_zero(tmp_path, monkeypatch):
    """When alpha=0, the correction factor is (ratio)^0 = 1 — no change to pixels."""
    image_dir = _make_casa_image_dir(tmp_path, "tt0.image")
    pixel_data = np.ones((4, 4), dtype=np.float64) * 77.0

    monkeypatch.setattr(
        "skasim.loaders.image_models._set_crval4_via_script",
        lambda work_dir, paths, freq: None,
    )

    written = _patch_require_casacore_and_table(
        monkeypatch, {str(image_dir): pixel_data},
    )

    result = adjust_spectral_reference(
        image_dir, old_ref_hz=5.99e9, new_ref_hz=1.284e9, alpha_map=np.zeros_like(pixel_data),
    )

    assert result == pytest.approx(1.284e9)
    np.testing.assert_allclose(written[str(image_dir)], 77.0, rtol=1e-10)


def test_adjust_spectral_reference_positive_alpha(tmp_path, monkeypatch):
    """Positive alpha (spectral index > 0) scales tt0 correctly."""
    image_dir = _make_casa_image_dir(tmp_path, "tt0.image")
    pixel_data = np.ones((2, 2), dtype=np.float64) * 10.0

    monkeypatch.setattr(
        "skasim.loaders.image_models._set_crval4_via_script",
        lambda work_dir, paths, freq: None,
    )

    written = _patch_require_casacore_and_table(
        monkeypatch, {str(image_dir): pixel_data},
    )

    # alpha = 0.5, old = 1.0 GHz, new = 2.0 GHz
    # factor = (2.0/1.0)^0.5 = sqrt(2)
    # tt0' = 10.0 * sqrt(2) ≈ 14.142
    result = adjust_spectral_reference(
        image_dir, old_ref_hz=1.0e9, new_ref_hz=2.0e9, alpha_map=np.full_like(pixel_data, 0.5),
    )

    expected = 10.0 * (2.0) ** 0.5
    np.testing.assert_allclose(written[str(image_dir)], expected, rtol=1e-6)


# ---------------------------------------------------------------------------
# prepare_casa_taylor_terms integration tests
# ---------------------------------------------------------------------------

def test_prepare_casa_taylor_terms_adjusts_spectral_reference(tmp_path, monkeypatch):
    """prepare_casa_taylor_terms adjusts reffreq and records a milestone."""
    tt0 = tmp_path / "model.tt0"
    tt1 = tmp_path / "model.tt1"
    tt0.mkdir()
    tt1.mkdir()
    (tt0 / "table.dat").write_text("table", encoding="utf-8")
    (tt1 / "table.dat").write_text("table", encoding="utf-8")

    cfg = SimConfig(
        models=[
            CasaTaylorTermsModelEntry(
                type="casa_taylor_terms",
                tt0=str(tt0),
                tt1=str(tt1),
                reference_frequency_hz=5.991892258e9,
            )
        ],
        observation=ObsConfig(frequency_mhz=1284.0, observation_time_s=1),
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)

    adjust_calls = []

    def fake_adjust(image_path, old_ref_hz, new_ref_hz, alpha_map=None):
        adjust_calls.append({
            "image_path": image_path,
            "old_ref_hz": old_ref_hz,
            "new_ref_hz": new_ref_hz,
            "alpha_map": alpha_map,
        })
        return new_ref_hz

    monkeypatch.setattr(
        "skasim.loaders.image_models.adjust_spectral_reference",
        fake_adjust,
    )

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casacore",
        lambda: MagicMock(),
    )

    from skasim.loaders.image_models import prepare_casa_taylor_terms
    product = prepare_casa_taylor_terms(ctx, cfg.models[0], 0)

    # reffreq should be observation centre frequency, not config reference
    assert product.reffreq == f"{1284.0e6}Hz"
    # two adjust calls: tt0 with element-wise alpha, tt1 scaled by the same factor
    assert len(adjust_calls) == 2
    assert adjust_calls[0]["alpha_map"] is not None  # tt0 gets pixel data correction
    assert adjust_calls[1]["alpha_map"] is not None  # tt1 also scaled to preserve tt1/tt0 ratio
    # milestone recorded
    milestone = [m for m in ctx.manifest.milestones if m.name == "adjusted_spectral_reference"]
    assert len(milestone) == 1
    assert milestone[0].details["old_reference_frequency_hz"] == pytest.approx(5.991892258e9)
    assert milestone[0].details["new_reference_frequency_hz"] == pytest.approx(1284.0e6)


def test_prepare_casa_taylor_terms_nterms1(tmp_path, monkeypatch):
    """casa_taylor_terms with tt0 only (nterms=1) gets CRVAL4-only adjustment."""
    tt0 = tmp_path / "model.tt0"
    tt0.mkdir()
    (tt0 / "table.dat").write_text("table", encoding="utf-8")

    cfg = SimConfig(
        models=[
            CasaTaylorTermsModelEntry(
                type="casa_taylor_terms",
                tt0=str(tt0),
                reference_frequency_hz=5.991892258e9,
            )
        ],
        observation=ObsConfig(frequency_mhz=1284.0, observation_time_s=1),
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)

    adjust_calls = []

    def fake_adjust(image_path, old_ref_hz, new_ref_hz, alpha_map=None):
        adjust_calls.append({"alpha_map": alpha_map})
        return new_ref_hz

    monkeypatch.setattr(
        "skasim.loaders.image_models.adjust_spectral_reference",
        fake_adjust,
    )

    from skasim.loaders.image_models import prepare_casa_taylor_terms
    product = prepare_casa_taylor_terms(ctx, cfg.models[0], 0)

    assert product.reffreq == f"{1284.0e6}Hz"
    assert len(adjust_calls) == 1
    assert adjust_calls[0]["alpha_map"] is None  # nterms=1, CRVAL4 only


# ---------------------------------------------------------------------------
# prepare_continuum_i_alpha_for_casa tests
# ---------------------------------------------------------------------------

def test_prepare_continuum_i_alpha_adjusts_spectral_reference(tmp_path, monkeypatch):
    """prepare_continuum_i_alpha_for_casa adjusts reffreq to observation centre."""
    from astropy.io import fits
    from astropy.wcs import WCS

    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10.0, 2.0]
    wcs.wcs.crpix = [8.0, 8.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    header = wcs.to_header()
    header["BUNIT"] = "Jy/pixel"

    stokes_path = tmp_path / "stokes_i.fits"
    alpha_path = tmp_path / "alpha.fits"
    fits.writeto(stokes_path, np.ones((16, 16), dtype=np.float32), header, overwrite=True)
    alpha_header = header.copy()
    alpha_header["BUNIT"] = "1"
    fits.writeto(alpha_path, np.full((16, 16), -0.7, dtype=np.float32), alpha_header, overwrite=True)

    cfg = SimConfig(
        models=[
            ContinuumIAlphaModelEntry(
                type="continuum_i_alpha",
                stokes_i=str(stokes_path),
                alpha=str(alpha_path),
                reference_frequency_hz=1.4e9,
            )
        ],
        observation=ObsConfig(frequency_mhz=1284.0, observation_time_s=1),
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)

    adjust_calls = []

    def fake_adjust(image_path, old_ref_hz, new_ref_hz, alpha_map=None):
        adjust_calls.append({
            "old_ref_hz": old_ref_hz,
            "new_ref_hz": new_ref_hz,
            "alpha_map": alpha_map,
        })
        return new_ref_hz

    monkeypatch.setattr(
        "skasim.loaders.image_models.adjust_spectral_reference",
        fake_adjust,
    )
    monkeypatch.setattr(
        "skasim.loaders.image_models.run_casa_importfits",
        lambda work_dir, images: None,
    )

    monkeypatch.setattr(
        "skasim.loaders.image_models.require_casacore",
        lambda: MagicMock(),
    )

    from skasim.loaders.image_models import prepare_continuum_i_alpha_for_casa
    product = prepare_continuum_i_alpha_for_casa(ctx, cfg.models[0], 0)

    assert product.reffreq == f"{1284.0e6}Hz"
    # two adjust calls; both use alpha_map=None because pixel scaling is done
    # directly on the intermediate FITS before importfits.
    assert len(adjust_calls) == 2
    assert adjust_calls[0]["alpha_map"] is None
    assert adjust_calls[0]["old_ref_hz"] == pytest.approx(1.4e9)
    assert adjust_calls[0]["new_ref_hz"] == pytest.approx(1284.0e6)
    assert adjust_calls[1]["alpha_map"] is None
    # milestone still records the mean spectral index used for pixel scaling
    milestone = [m for m in ctx.manifest.milestones if m.name == "adjusted_spectral_reference"]
    assert len(milestone) == 1
    assert milestone[0].details["alpha_mean"] == pytest.approx(-0.7)
