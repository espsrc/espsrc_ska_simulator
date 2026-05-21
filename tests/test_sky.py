"""tests/test_sky.py for Source and SkyModel classes"""

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord
from astropy.table import Table

from skasim.sky import SkyModel, Source

# -----------------------------------------------------------------------------
# Source creation
# -----------------------------------------------------------------------------


def test_source_creation_basic():
    """basic constructor."""
    s = Source(ra=10.0, dec=20.0, I=1.5)
    assert s.ra.value == 10.0
    assert s.dec.value == 20.0
    assert s.I.value == 1.5
    assert s.flux.value == 1.5
    assert s.spec_index == 0


def test_source_creation_with_units():
    """constructor with explicit units."""
    s = Source(
        ra=10 * u.deg,
        dec=20 * u.deg,
        I=5 * u.Jy,
        Q=0.1 * u.Jy,
        ref_freq=1.4e9 * u.Hz,
        spec_index=-0.7,
    )
    assert s.I.unit == u.Jy
    assert s.Q.value == pytest.approx(0.1)
    assert s.ref_freq.to(u.GHz).value == pytest.approx(1.4)
    assert s.spec_index == pytest.approx(-0.7)


def test_source_creation_wrong_units_raises():
    """incompatible units"""
    with pytest.raises(ValueError):
        Source(ra=10.0, dec=20.0, I=5 * u.m)  # should be Jy


def test_source_creation_all_parameters():
    """optional parameters."""
    s = Source(
        ra=150.0,
        dec=2.5,
        I=3.0,
        Q=0.2,
        U=0.1,
        V=0.05,
        ref_freq=200e6,
        spec_index=-0.7,
        rot_meas=0.5,
        major_axis=5.0,
        minor_axis=3.0,
        pa=45.0,
        true_redshift=0.1,
        obs_redshift=0.05,
        obj_id="custom_id",
        resolved=True,
        isl_rms=0.01,
    )
    assert s.obj_id == "custom_id"
    assert s.resolved is True
    assert s.true_redshift == 0.1
    assert s.obs_redshift == 0.05
    assert s.isl_rms.to(u.Jy).value == pytest.approx(0.01)


# -----------------------------------------------------------------------------
# Source defaults
# -----------------------------------------------------------------------------


def test_source_defaults():
    """optional fields with defaults."""
    s = Source(ra=0, dec=0, I=1.0)
    assert s.Q.value == 0.0
    assert s.U.value == 0.0
    assert s.V.value == 0.0
    assert s.rot_meas.value == 0.0
    assert s.isl_rms.value == 0.0
    assert s.major_axis.value == 0.0
    assert s.minor_axis.value == 0.0
    assert s.pa.value == 0.0


# -----------------------------------------------------------------------------
# obj_id
# -----------------------------------------------------------------------------


def test_source_auto_obj_id():
    """obj_id autogeneration from coords."""
    s = Source(ra=10.5, dec=20.5, I=1.0)
    assert "10.5" in s.obj_id
    assert "20.5" in s.obj_id


def test_source_custom_obj_id():
    """custom obj_id."""
    s = Source(ra=0, dec=0, I=1.0, obj_id="my_source")
    assert s.obj_id == "my_source"


# -----------------------------------------------------------------------------
# resolved flag
# -----------------------------------------------------------------------------


def test_source_resolved_true():
    """resolved source with axes."""
    s = Source(ra=0, dec=0, I=1.0, resolved=True, major_axis=5.0, minor_axis=3.0)
    assert s.resolved is True
    assert s.major_axis.to(u.arcsec).value == 5.0
    assert s.minor_axis.to(u.arcsec).value == 3.0


# -----------------------------------------------------------------------------
# to_json / from_json roundtrip
# -----------------------------------------------------------------------------
def test_source_json_roundtrip():
    """safe serialization"""
    original = Source(
        ra=150.0,
        dec=2.5,
        I=3.0,
        Q=0.2,
        U=0.1,
        V=0.05,
        ref_freq=200e6,
        spec_index=-0.7,
        major_axis=5.0,
        minor_axis=3.0,
        pa=45.0,
        true_redshift=0.1,
        obs_redshift=0.05,
    )
    d = original.to_json()
    restored = Source.from_json(d)
    assert restored.ra.value == pytest.approx(150.0)
    assert restored.dec.value == pytest.approx(2.5)
    assert restored.I.value == pytest.approx(3.0)
    assert restored.spec_index == pytest.approx(-0.7)
    assert restored.major_axis.to(u.arcsec).value == pytest.approx(5.0)
    assert restored.true_redshift == pytest.approx(0.1)
    assert restored.obs_redshift == pytest.approx(0.05)


def test_source_to_json_hmsdms_format():
    """verify to_json supports hmsdms coordinate formatting."""
    s = Source(ra=150.0, dec=2.5, I=1.0)
    d = s.to_json(coords_fmt="hmsdms")
    assert isinstance(d["ra"], str)
    assert ":" in d["ra"]


# -----------------------------------------------------------------------------
# from_array
# -----------------------------------------------------------------------------


def test_source_from_array_3_elements():
    """source constructor from minimal array: ra, dec, I."""
    arr = [150.0, 2.5, 3.0]
    s = Source.from_array(arr)
    assert s.ra.value == 150.0
    assert s.I.value == 3.0
    assert s.Q.value == 0.0


def test_source_from_array_6_elements():
    """source constructor with stoke parameters."""
    arr = [150.0, 2.5, 3.0, 0.2, 0.1, 0.05]
    s = Source.from_array(arr)
    assert s.Q.value == 0.2
    assert s.U.value == 0.1
    assert s.V.value == 0.05


def test_source_from_array_12_elements():
    """source constructor with full parameters."""
    arr = [150.0, 2.5, 3.0, 0.2, 0.1, 0.05, 200e6, -0.7, 0.0, 5.0, 3.0, 45.0]
    s = Source.from_array(arr)
    assert s.ref_freq.to(u.Hz).value == pytest.approx(200e6)
    assert s.spec_index == -0.7
    assert s.major_axis.to(u.arcsec).value == 5.0


def test_source_from_array_14_elements():
    """same, with redshift"""
    arr = [150.0, 2.5, 3.0, 0.2, 0.1, 0.05, 200e6, -0.7, 0.0, 5.0, 3.0, 45.0, 0.1, 0.05]
    s = Source.from_array(arr)
    assert s.true_redshift == 0.1
    assert s.obs_redshift == 0.05


def test_source_from_array_16_elements():
    """same, with resolved flag and rms."""
    arr = [
        150.0,
        2.5,
        3.0,
        0.2,
        0.1,
        0.05,
        200e6,
        -0.7,
        0.0,
        5.0,
        3.0,
        45.0,
        0.1,
        0.05,
        True,
        0.01,
    ]
    s = Source.from_array(arr)
    assert s.resolved is True
    assert s.isl_rms.to(u.Jy).value == pytest.approx(0.01)


def test_source_from_array_too_short_raises():
    """insuficient info for source creation."""
    with pytest.raises(ValueError, match="at least 3 elements"):
        Source.from_array([1.0, 2.0])


def test_source_from_array_unknown_length_raises():
    """invalid array for source creation."""
    with pytest.raises(ValueError, match="must have 3, 6, or 12"):
        Source.from_array([1.0] * 5)


def test_source_from_array_with_colnames():
    """source creation with explicit column mapping"""
    arr = [150.0, 2.5, 3.0]
    s = Source.from_array(arr, colnames=["ra", "dec", "I"])
    assert s.ra.value == 150.0
    assert s.Q.value == 0.0


# -----------------------------------------------------------------------------
# from_table_in_fits
# -----------------------------------------------------------------------------


def test_source_from_table_in_fits_basic():
    """standard astropy table."""
    t = Table()
    t["RA"] = [10.0, 11.0]
    t["DEC"] = [20.0, 21.0]
    t["STK_I"] = [1.0, 2.0]
    sources = Source.from_table_in_fits(t)
    assert len(sources) == 2
    assert sources[0].ra.value == 10.0
    assert sources[0].I.value == 1.0


def test_source_from_table_in_fits_with_id():
    """astropy table with column ID."""
    t = Table()
    t["RA"] = [10.0]
    t["DEC"] = [20.0]
    t["STK_I"] = [1.0]
    t["ID"] = ["src_A"]
    sources = Source.from_table_in_fits(t)
    assert sources[0].obj_id == "src_A"


def test_source_from_table_in_fits_alternative_columns():
    """resolution of non-standard colnames."""
    t = Table()
    t["RA"] = [10.0]
    t["DEC"] = [20.0]
    t["S_INT"] = [1.0]  # alias for STK_I
    t["NU_EFF"] = [100e6]  # alias for REFFREQ
    t["SPECIDX"] = [-0.7]
    sources = Source.from_table_in_fits(t)
    assert sources[0].I.value == 1.0
    assert sources[0].ref_freq.to(u.Hz).value == pytest.approx(100e6)
    assert sources[0].spec_index == -0.7


# -----------------------------------------------------------------------------
# to_sky_model
# -----------------------------------------------------------------------------


def test_source_to_sky_model_reduced():
    """from source to sky model."""
    s = Source(ra=10.0, dec=20.0, I=2.0)
    tup = s.to_sky_model(reduced_form=True)
    assert tup == (10.0, 20.0, 2.0)


def test_source_to_sky_model_full():
    """same with many parameters."""
    s = Source(ra=10.0, dec=20.0, I=2.0, Q=0.1, spec_index=-0.5, true_redshift=0.1)
    tup = s.to_sky_model(reduced_form=False)
    assert len(tup) == 14
    assert tup[2] == pytest.approx(2.0)
    assert tup[7] == pytest.approx(-0.5)


# -----------------------------------------------------------------------------
# get_flux / flux property
# -----------------------------------------------------------------------------


def test_source_get_flux_no_frequency():
    """without freq argument."""
    s = Source(ra=0, dec=0, I=2.0, ref_freq=100e6)
    assert s.get_flux().value == pytest.approx(2.0)


def test_source_get_flux_at_frequency():
    """spectral index scaling is applied correctly."""
    s = Source(ra=0, dec=0, I=1.0, ref_freq=100e6, spec_index=-0.7)
    flux_200 = s.get_flux(freq=200e6)
    expected = 1.0 * (200e6 / 100e6) ** -0.7
    assert flux_200.value == pytest.approx(expected)


def test_source_get_flux_custom_alpha():
    """alpha parameter overrides the source spectral index."""
    s = Source(ra=0, dec=0, I=1.0, ref_freq=100e6, spec_index=-0.7)
    flux = s.get_flux(freq=200e6, alpha=-0.5)
    expected = 1.0 * (200e6 / 100e6) ** -0.5
    assert flux.value == pytest.approx(expected)


def test_source_flux_property():
    """flux property is an alias for Stokes I."""
    s = Source(ra=0, dec=0, I=3.0)
    assert s.flux.value == 3.0


# -----------------------------------------------------------------------------
# to_fits_fmt
# -----------------------------------------------------------------------------


def test_source_to_fits_fmt():
    """FITS dictionary keys match the expected column names."""
    s = Source(ra=10.0, dec=20.0, I=1.5, Q=0.1, spec_index=-0.7)
    d = s.to_fits_fmt()
    assert d["RA"] == 10.0
    assert d["DEC"] == 20.0
    assert d["STK_I"] == 1.5
    assert d["STK_Q"] == 0.1
    assert d["SPECIDX"] == -0.7


# -----------------------------------------------------------------------------
# SkyModel class
# -----------------------------------------------------------------------------


def test_skymodel_empty():
    """empty SkyModel has no sources."""
    sm = SkyModel()
    assert sm.sources is None or len(sm.sources) == 0


def test_skymodel_from_array():
    """SkyModel can be built from arrays."""
    arr = np.array([[10.0, 20.0, 1.0], [11.0, 21.0, 2.0]])
    sm = SkyModel(arr)
    assert len(sm.sources) == 2
    center = sm.get_center()
    assert pytest.approx(center.ra.value) == 10.5
    assert pytest.approx(center.dec.value) == 20.5


def test_skymodel_get_center_cached():
    """if phase_center is set, get_center returns it."""
    sm = SkyModel()
    sm.phase_center = SkyCoord(ra=15.0 * u.deg, dec=25.0 * u.deg, frame="icrs")
    center = sm.get_center()
    assert center.ra.value == 15.0
    assert center.dec.value == 25.0


def test_skymodel_get_center_no_sources_raises():
    """get_center raises AttributeError when there is nothing to centre on."""
    sm = SkyModel()
    with pytest.raises(AttributeError):
        sm.get_center()


def test_skymodel_to_json():
    """to_json serialises the source list."""
    s1 = Source(ra=10.0, dec=20.0, I=1.0)
    s2 = Source(ra=11.0, dec=21.0, I=2.0)
    arr = np.array(
        [s1.to_sky_model(reduced_form=False), s2.to_sky_model(reduced_form=False)]
    )
    sm = SkyModel(arr)
    json_data = sm.to_json()
    assert len(json_data) == 2
    assert json_data[0]["ra"] == 10.0


def test_skymodel_to_json_empty_sources():
    """empty source list yields empty JSON list."""
    sm = SkyModel(np.empty((0, 14)))
    assert sm.to_json() == []


def test_skymodel_from_json_roundtrip():
    """from_json reconstructs a SkyModel identically."""
    s1 = Source(ra=10.0, dec=20.0, I=1.0)
    s2 = Source(ra=11.0, dec=21.0, I=2.0)
    arr = np.array(
        [s1.to_sky_model(reduced_form=False), s2.to_sky_model(reduced_form=False)]
    )
    sm = SkyModel(arr)
    json_data = sm.to_json()

    sm2 = SkyModel.from_json(json_data)
    assert sm2 is not None
    assert len(sm2.sources) == 2
    center = sm2.get_center()
    assert pytest.approx(center.ra.value) == 10.5


def test_skymodel_from_json_bad_data_returns_none():
    """malformed JSON is caught and returns None."""
    sm = SkyModel.from_json([{"ra": "invalid"}])
    assert sm is None


def test_skymodel_from_json_empty_list_returns_none():
    """Empty list triggers exception handling and returns None."""
    sm = SkyModel.from_json([])
    assert sm is None


def test_skymodel_from_fits_table_missing_file():
    """missing FITS file raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        SkyModel.from_fits_table("/nonexistent/file.fits")
