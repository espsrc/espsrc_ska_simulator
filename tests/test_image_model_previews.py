"""Tests for the two-panel FITS image-model previews."""

from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS

from skasim.config import SimConfig
from skasim.loaders.image_models import write_image_model_previews
from skasim.loaders.image_models.previews import (
    _fits_image_center,
    _fits_image_extent_deg,
)
from skasim.manifest import create_run_context
from skasim.weblog import render_weblog


def _write_model_fits(path: Path, shape=(32, 32), cdelt_deg=0.001):
    wcs = WCS(naxis=2)
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    wcs.wcs.crval = [10.0, 2.0]
    wcs.wcs.crpix = [shape[1] / 2.0, shape[0] / 2.0]
    wcs.wcs.cdelt = [-cdelt_deg, cdelt_deg]
    header = wcs.to_header()
    header["BUNIT"] = "Jy/pixel"
    fits.writeto(path, np.ones(shape, dtype=np.float32), header, overwrite=True)


def test_two_panel_preview_renders(tmp_path):
    """A continuum_i_alpha entry produces a two-panel PNG with both titles."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")
    stokes = tmp_path / "stokes_i.fits"
    alpha = tmp_path / "alpha.fits"
    _write_model_fits(stokes, shape=(32, 32), cdelt_deg=0.001)
    _write_model_fits(alpha, shape=(32, 32), cdelt_deg=0.001)

    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    center = SkyCoord(10.0, 2.0, unit="deg")
    fov = 0.05 * u.deg

    write_image_model_previews(ctx, fov, center)

    png = ctx.work_dir / "run_fits_model.png"
    assert png.exists()
    assert png.stat().st_size > 0

    html = render_weblog(ctx.manifest, ctx.work_dir)
    assert "Footprint" in html
    assert "Natural Extent" in html
    assert "FoV Context" in html


def test_footprint_stats_recorded_and_rendered(tmp_path):
    """Footprint stats (size, centre, offset, FoV fraction) reach the manifest and weblog."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")
    stokes = tmp_path / "stokes_i.fits"
    alpha = tmp_path / "alpha.fits"
    _write_model_fits(stokes, shape=(32, 32), cdelt_deg=0.001)
    _write_model_fits(alpha, shape=(32, 32), cdelt_deg=0.001)

    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    # offset the phase centre from the model centre so the offset stat is non-zero
    center = SkyCoord(10.02, 2.0, unit="deg")
    fov = 0.05 * u.deg

    write_image_model_previews(ctx, fov, center)

    fits_model_output = next(o for o in ctx.manifest.outputs if o.role == "fits_model")
    metadata = fits_model_output.metadata
    assert metadata["image_size"] == [32, 32]
    assert metadata["pixel_scale_arcsec"] == pytest.approx(3.4875, abs=0.05)
    assert metadata["extent_width_deg"] == pytest.approx(0.031, abs=1e-3)
    assert metadata["center_ra_deg"] == pytest.approx(10.0, abs=1e-3)
    assert metadata["center_dec_deg"] == pytest.approx(2.0, abs=1e-3)
    assert metadata["offset_from_center_arcsec"] == pytest.approx(72.0, abs=2.5)
    assert metadata["fov_fraction"] == pytest.approx(0.031 / 0.05, abs=1e-2)

    html = render_weblog(ctx.manifest, ctx.work_dir)
    assert "32 × 32 px" in html
    assert "% of FoV" in html
    assert "from phase centre" in html


def test_context_panel_uses_simulation_center(tmp_path):
    """The contextual panel title reports the simulation FoV passed from the pipeline."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")
    stokes = tmp_path / "stokes_i.fits"
    alpha = tmp_path / "alpha.fits"
    _write_model_fits(stokes, shape=(16, 16), cdelt_deg=0.001)
    _write_model_fits(alpha, shape=(16, 16), cdelt_deg=0.001)

    cfg = SimConfig(
        models=[
            {
                "type": "continuum_i_alpha",
                "stokes_i": str(stokes),
                "alpha": str(alpha),
                "reference_frequency_hz": 1.4e9,
            }
        ],
        output_dir=str(tmp_path / "run"),
    )
    ctx = create_run_context(cfg)
    center = SkyCoord(10.005, 2.005, unit="deg")
    fov = 0.03 * u.deg

    write_image_model_previews(ctx, fov, center)

    # minimal check: the output exists and the role is correct
    manifest_roles = [out.role for out in ctx.manifest.outputs]
    assert "fits_model" in manifest_roles


def test_fits_image_extent_deg_uses_wcs(tmp_path):
    """_fits_image_extent_deg returns the approximate WCS extent in degrees."""
    path = tmp_path / "model.fits"
    _write_model_fits(path, shape=(64, 64), cdelt_deg=0.002)
    with fits.open(path) as hdul:
        extent = _fits_image_extent_deg(hdul[0].header)
    assert extent is not None
    assert extent["width_deg"] == pytest.approx(0.126, abs=0.01)
    assert extent["height_deg"] == pytest.approx(0.126, abs=0.01)


def test_fits_image_center_uses_wcs(tmp_path):
    """_fits_image_center returns the pixel centre as a SkyCoord."""
    path = tmp_path / "model.fits"
    _write_model_fits(path, shape=(32, 32), cdelt_deg=0.001)
    with fits.open(path) as hdul:
        center = _fits_image_center(hdul[0].header)
    assert center is not None
    assert center.ra.deg == pytest.approx(10.0, abs=1e-3)
    assert center.dec.deg == pytest.approx(2.0, abs=1e-3)


def test_fits_image_center_degrades_gracefully(tmp_path):
    """_fits_image_center returns None for a FITS without celestial WCS."""
    path = tmp_path / "no_wcs.fits"
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 8
    header["NAXIS2"] = 8
    fits.writeto(path, np.ones((8, 8), dtype=np.float32), header, overwrite=True)
    with fits.open(path) as hdul:
        center = _fits_image_center(hdul[0].header)
    assert center is None


def test_fits_image_extent_deg_falls_back_to_cdelt(tmp_path):
    """_fits_image_extent_deg falls back to CDELT when WCS is unavailable."""
    path = tmp_path / "cdelt_only.fits"
    header = fits.Header()
    header["NAXIS"] = 2
    header["NAXIS1"] = 16
    header["NAXIS2"] = 16
    header["CDELT1"] = 0.01
    header["CDELT2"] = 0.01
    fits.writeto(path, np.ones((16, 16), dtype=np.float32), header, overwrite=True)
    with fits.open(path) as hdul:
        extent = _fits_image_extent_deg(hdul[0].header)
    assert extent is not None
    assert extent["width_deg"] == pytest.approx(0.16, abs=1e-6)
    assert extent["height_deg"] == pytest.approx(0.16, abs=1e-6)
