"""Runtime dependency helpers."""

from __future__ import annotations

import importlib
from types import ModuleType


KARABO_INSTALL_MESSAGE = (
    "Karabo is required for full simulation execution. Install and activate the "
    "conda skasim environment with `conda env create -f environment.yml` before "
    "running simulations."
)

OSKAR_INSTALL_MESSAGE = (
    "OSKAR is required for FITS image ingestion. Install and activate the "
    "conda skasim environment with `conda env create -f environment.yml` before "
    "using --fits-image."
)


class KaraboRuntimeError(RuntimeError):
    """Raised when full simulation execution needs Karabo but it is unavailable."""


class OskarRuntimeError(RuntimeError):
    """Raised when FITS image ingestion needs OSKAR but it is unavailable."""


def require_karabo_module(module_name: str) -> ModuleType:
    """Import a Karabo module or raise the supported runtime setup message."""
    try:
        return importlib.import_module(module_name)
    except ImportError as exc:
        # Only translate missing Karabo imports; dependency errors inside Karabo
        # modules should keep their original traceback.
        if exc.name == "karabo" or (exc.name or "").startswith("karabo."):
            raise KaraboRuntimeError(KARABO_INSTALL_MESSAGE) from exc
        raise


def require_oskar_module() -> ModuleType:
    """Import the top-level `oskar` module or raise a clear error."""
    try:
        return importlib.import_module("oskar")
    except ImportError as exc:
        raise OskarRuntimeError(OSKAR_INSTALL_MESSAGE) from exc
