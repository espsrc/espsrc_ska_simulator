"""Tests for FitsImageLoader — FITS image ingestion via OSKAR."""

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.io import fits

from skasim.loaders import FitsImageLoader
from skasim.runtime import OskarRuntimeError


class _FakeOskarSky:
    """Minimal mock of oskar.Sky for unit tests."""

    def __init__(self, array):
        self._array = array

    def to_array(self):
        return self._array

    @classmethod
    def from_fits_file(cls, fpath, **kwargs):
        return cls(kwargs.pop("_test_array", np.array([[10.0, 20.0, 1.0]])))


def _make_fake_fits_image(tmp_path, naxis1=64, naxis2=64, crval3=None):
    """Write a minimal 2D FITS image to disk and return its path."""
    data = np.zeros((naxis2, naxis1), dtype=np.float32)
    data[30:34, 30:34] = 1.0  # small square source
    header = {
        "NAXIS": 2,
        "NAXIS1": naxis1,
        "NAXIS2": naxis2,
        "CRPIX1": naxis1 / 2,
        "CRPIX2": naxis2 / 2,
        "CRVAL1": 10.0,
        "CRVAL2": 20.0,
        "CDELT1": 0.01,
        "CDELT2": 0.01,
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CUNIT1": "deg",
        "CUNIT2": "deg",
    }
    if crval3 is not None:
        header["NAXIS"] = 3
        header["NAXIS3"] = 1
        header["CRPIX3"] = 1
        header["CRVAL3"] = crval3
        header["CDELT3"] = 1.0
        header["CTYPE3"] = "FREQ"
        header["CUNIT3"] = "Hz"
        data = data.reshape(1, naxis2, naxis1)

    hdu = fits.PrimaryHDU(data, header=fits.Header(header))
    path = tmp_path / "test_image.fits"
    hdu.writeto(str(path), overwrite=True)
    return path


# -----------------------------------------------------------------------------
# Happy path
# -----------------------------------------------------------------------------


def test_fits_image_loader_happy_path(monkeypatch, tmp_path):
    """FITS image loads via mocked OSKAR and produces a SkyModel."""
    fpath = _make_fake_fits_image(tmp_path)

    fake_array = np.array(
        [
            [10.0, 20.0, 1.0],
            [10.01, 20.01, 0.5],
        ]
    )

    class _MockSky:
        def __init__(self, array):
            self._array = array

        def to_array(self):
            return self._array

        @classmethod
        def from_fits_file(cls, fpath, **kwargs):
            return cls(fake_array)

    def _mock_oskar():
        return type("MockOskar", (), {"Sky": _MockSky})()

    monkeypatch.setattr(
        "skasim.loaders.fits_image.require_oskar_module",
        _mock_oskar,
    )

    loader = FitsImageLoader(fpath, fallback_freq_mhz=700.0)
    sky_model = loader.load()

    assert sky_model is not None
    assert sky_model.phase_center is not None
    # Center should be near the image center from WCS
    assert pytest.approx(sky_model.phase_center.ra.to(u.deg).value, abs=1e-3) == 10.00
    assert pytest.approx(sky_model.phase_center.dec.to(u.deg).value, abs=1e-3) == 20.00


# -----------------------------------------------------------------------------
# Frequency extraction
# -----------------------------------------------------------------------------


def test_extracts_freq_from_header(monkeypatch, tmp_path):
    """When CTYPE3=FREQ exists, loader notes the frequency."""
    fpath = _make_fake_fits_image(tmp_path, crval3=1.4e9)

    monkeypatch.setattr(
        "skasim.loaders.fits_image.require_oskar_module",
        lambda: type(
            "MockOskar",
            (),
            {
                "Sky": type(
                    "Sky",
                    (_FakeOskarSky,),
                    {
                        "from_fits_file": lambda cls, fp, **kw: cls(
                            np.array([[10.0, 20.0, 1.0]])
                        ),
                    },
                )
            },
        ),
    )

    loader = FitsImageLoader(fpath, fallback_freq_mhz=700.0)
    freq_hz = loader._extract_or_fallback_freq_hz()
    assert freq_hz == pytest.approx(1.4e9)


def test_fallback_when_no_freq_axis(monkeypatch, tmp_path):
    """When no FREQ axis, fallback frequency is used."""
    fpath = _make_fake_fits_image(tmp_path)

    monkeypatch.setattr(
        "skasim.loaders.fits_image.require_oskar_module",
        lambda: type(
            "MockOskar",
            (),
            {
                "Sky": type(
                    "Sky",
                    (_FakeOskarSky,),
                    {
                        "from_fits_file": lambda cls, fp, **kw: cls(
                            np.array([[10.0, 20.0, 1.0]])
                        ),
                    },
                )
            },
        ),
    )

    loader = FitsImageLoader(fpath, fallback_freq_mhz=1400.0)
    freq_hz = loader._extract_or_fallback_freq_hz()
    assert freq_hz == pytest.approx(1.4e9)


# -----------------------------------------------------------------------------
# Errors
# -----------------------------------------------------------------------------


def test_missing_file_raises():
    """Non-existent FITS raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        FitsImageLoader("/nonexistent/image.fits").load()


def test_missing_oskar_raises(monkeypatch, tmp_path):
    """When oskar is absent, OskarRuntimeError is raised."""
    fpath = _make_fake_fits_image(tmp_path)
    monkeypatch.setattr(
        "skasim.loaders.fits_image.require_oskar_module",
        lambda: (_ for _ in ()).throw(OskarRuntimeError("missing SDK")),
    )
    with pytest.raises(OskarRuntimeError, match="missing SDK"):
        FitsImageLoader(fpath).load()


# -----------------------------------------------------------------------------
# Phase center
# -----------------------------------------------------------------------------


def test_phase_center_from_wcs_list_return(monkeypatch, tmp_path):
    """Phase centre works when pixel_to_world returns a list instead of tuple."""
    fpath = _make_fake_fits_image(tmp_path, naxis1=100, naxis2=100)

    mock_coord = SkyCoord(ra=283.51816667 * u.deg, dec=3.59527778 * u.deg, frame="fk5")
    mock_list = [mock_coord]

    def fake_pixel_to_world(*args, **kwargs):
        return mock_list

    monkeypatch.setattr(
        "skasim.loaders.fits_image.WCS",
        lambda h: type(
            "FakeWCS",
            (),
            {
                "pixel_n_dim": 2,
                "pixel_to_world": fake_pixel_to_world,
            },
        )(),
    )

    loader = FitsImageLoader(fpath)
    center = loader._compute_phase_center()

    assert isinstance(center, SkyCoord)
    assert pytest.approx(center.ra.to(u.deg).value, abs=1e-6) == 283.51816667
    assert pytest.approx(center.dec.to(u.deg).value, abs=1e-6) == 3.59527778
