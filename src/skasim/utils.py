"""helpers extracted from legacy scripts/utils.py."""

import json
import os
import sys
from datetime import datetime
from typing import Optional

import astropy.units as u
import numpy as np
from loguru import logger


def init_logger(log_file: Optional[str] = None) -> None:
    """Cconfigure loguru: stderr + optional log file."""
    logger.remove()
    logger.add(
        sys.stderr,
        colorize=True,
        level="INFO",
    )
    if log_file:
        logger.add(
            log_file,
            colorize=False,
            level="INFO",
            enqueue=True,
        )


def define_extra_units() -> None:
    extra_units = [
        u.def_unit("JY", 1 * u.Jy),
        u.def_unit("DEG", 1 * u.deg),
        u.def_unit("JY/BEAM", 1 * u.Jy / u.sr),
        u.def_unit("HZ", 1 * u.Hz),
    ]
    u.add_enabled_units(extra_units)


def mapping_unit(unit_str: Optional[str]) -> Optional[str]:
    if unit_str is None:
        return None
    return {
        "JY": "Jy",
        "DEG": "deg",
        "JY/BEAM": "Jy/beam",
        "HZ": "Hz",
    }.get(unit_str, unit_str)


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.bool_):
            return bool(obj)
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)





# TODO: define elsewhere or check if karabo provides
DIAMETERS = {
    "ALMA": 25 * u.m,
    "APEX": 12 * u.m,
    "ATCA": 22 * u.m,
    "CARMA": 10.4 * u.m,
    "GBT": 100 * u.m,
    "GMRT": 45 * u.m,
    "IRAM30M": 30 * u.m,
    "JCMT": 15 * u.m,
    "LOFAR": 25 * u.m,
    "MEERKAT": 13.5 * u.m,
    "MRT": 30 * u.m,
    "NRAO12M": 12 * u.m,
    "NRAO20M": 20 * u.m,
    "NRAO40M": 40 * u.m,
    "NRAO45M": 45 * u.m,
    "NRAO90M": 90 * u.m,
    "PARKES": 64 * u.m,
    "SMA": 6.5 * u.m,
    "SKA1LOW": 38 * u.m,
    "SKA1MID": 15 * u.m,
}


def get_diameter(telescope_name: str):
    name = telescope_name.upper()
    if name in DIAMETERS:
        return DIAMETERS[name]
    if "SKA" in name or "SKA1" in name:
        if "LOW" in name:
            return DIAMETERS["SKA1LOW"]
        if "MID" in name:
            return DIAMETERS["SKA1MID"]
        raise ValueError(
            f"Telescope {telescope_name} not found. "
            f"Available: {', '.join(DIAMETERS.keys())}"
        )
    raise ValueError(
        f"Telescope {telescope_name} not found. "
        f"Available: {', '.join(DIAMETERS.keys())}"
    )
