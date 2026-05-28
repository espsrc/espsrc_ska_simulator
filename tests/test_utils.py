"""tests/test_utils.py"""

import json
import sys
import tempfile
from pathlib import Path

import astropy.units as u
import numpy as np
import pytest

from skasim.utils import (
    DIAMETERS,
    NpEncoder,
    define_extra_units,
    get_diameter,
    mapping_unit,
    # printlog,
    # show_exc,
)

# -----------------------------------------------------------------------------
# mapping_unit
# -----------------------------------------------------------------------------


def test_mapping_unit_none():
    assert mapping_unit(None) is None


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("JY", "Jy"),
        ("DEG", "deg"),
        ("JY/BEAM", "Jy/beam"),
        ("HZ", "Hz"),
    ],
)
def test_mapping_unit_known(raw, expected):
    assert mapping_unit(raw) == expected


def test_mapping_unit_unknown():
    assert mapping_unit("foo") == "foo"


# -----------------------------------------------------------------------------
# NpEncoder
# -----------------------------------------------------------------------------


def test_npencoder_bool():
    assert json.dumps({"ok": np.bool_(True)}, cls=NpEncoder) == '{"ok": true}'


def test_npencoder_integer():
    val = np.int64(42)
    assert json.dumps({"n": val}, cls=NpEncoder) == '{"n": 42}'


def test_npencoder_floating():
    val = np.float64(3.14)
    assert json.dumps({"x": val}, cls=NpEncoder) == '{"x": 3.14}'


def test_npencoder_ndarray():
    arr = np.array([1, 2, 3])
    result = json.loads(json.dumps({"a": arr}, cls=NpEncoder))
    assert result["a"] == [1, 2, 3]


# -----------------------------------------------------------------------------
# get_diameter
# -----------------------------------------------------------------------------


@pytest.mark.parametrize("name", list(DIAMETERS.keys()))
def test_get_diameter_exact(name):
    diam = get_diameter(name)
    assert isinstance(diam, u.Quantity)
    assert diam.unit.physical_type == "length"
    assert diam > 0


def test_get_diameter_lowercase():
    assert get_diameter("ska1mid") == DIAMETERS["SKA1MID"]


def test_get_diameter_ska_fallback_low():
    assert get_diameter("SKA1-LOW") == DIAMETERS["SKA1LOW"]


def test_get_diameter_ska_fallback_mid():
    assert get_diameter("SKA1-MID") == DIAMETERS["SKA1MID"]


def test_get_diameter_unknown_raises():
    with pytest.raises(ValueError, match="Telescope.*not found"):
        get_diameter("FAKE_SCOPE")


# -----------------------------------------------------------------------------
# define_extra_units
# -----------------------------------------------------------------------------


def test_define_extra_units_idempotent():
    # calling twice should not crash
    define_extra_units()
    define_extra_units()
    # TODO: verify behavior
    assert u.Unit("JY").is_equivalent(u.Jy)
