"""loaders/fits_image.py — Load a SkyModel from a FITS image via OSKAR."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from loguru import logger
import xarray as xr

from ..runtime import require_oskar_module
from ..sky import SkyModel


class FitsImageLoader:
    """Load a `SkyModel` from a FITS image using `oskar.Sky`.

    The protocol is:
        oskar.Sky.from_fits_file(path) → to_array() → np.ndarray → SkyModel

    This loader encapsulates the entire OSKAR dependency, allowing lazy import
    via `require_oskar_module` so that the module can be imported even when
    OSKAR is not installed.
    """

    def __init__(
        self,
        fpath: str | Path,
        fallback_freq_mhz: Optional[float] = None,
    ) -> None:
        self.fpath = Path(fpath)
        self.fallback_freq_mhz = fallback_freq_mhz

    def load(self) -> SkyModel:
        """Convert the FITS image to a Karabo SkyModel."""
        if not self.fpath.exists():
            raise FileNotFoundError(f"FITS image not found: {self.fpath}")

        oskar = require_oskar_module()
        freq_hz = self._extract_or_fallback_freq_hz()

        logger.info(f"Loading FITS image via OSKAR: {self.fpath}")
        oskar_sky = oskar.Sky.from_fits_file(
            str(self.fpath),
        )
        data_array = oskar_sky.to_array()
        np_array = self._to_numpy_array(data_array)

        sky_model = SkyModel(sources=xr.DataArray(np_array))
        sky_model.phase_center = self._compute_phase_center()
        sky_model.get_center()

        logger.info(
            f"FITS image loaded: {self.fpath.name} — ref_freq={freq_hz:.3e} Hz"
        )
        return sky_model

    def _extract_or_fallback_freq_hz(self) -> float:
        """Read CRVAL3 from the FITS header or fall back to observation frequency.

        In the meantime, this method is informative only and does not
        mutate the loaded sky model — OSKAR already embeds frequencies
        into the array during from_fits_file().
        """
        header = self._read_primary_header()

        freq_hz: Optional[float] = None
        for i in range(1, header.get("NAXIS", 0) + 1):
            ctype = header.get(f"CTYPE{i}", "")
            if "FREQ" in ctype.upper():
                freq_hz = header.get(f"CRVAL{i}", None)
                if freq_hz is not None:
                    break

        if freq_hz is None:
            freq_hz = (self.fallback_freq_mhz or 700.0) * 1e6
            logger.warning(
                f"No FREQ axis in {self.fpath.name}; "
                f"using fallback frequency {freq_hz:.3e} Hz"
            )
        return freq_hz

    def _read_primary_header(self) -> dict:
        """Return the primary HDU header as a dict-like object."""
        with fits.open(str(self.fpath)) as hdul:
            return hdul[0].header

    def _compute_phase_center(self) -> SkyCoord:
        """Derive phase centre from the FITS WCS.

        Handles multi-dimensional WCS by matching the number of pixel
        arguments to the dimensions exposed by the parsed WCS.
        """
        header = self._read_primary_header()
        wcs = WCS(header)
        spatial_pixel_dim = wcs.pixel_n_dim
        # Note: wcs.world_n_dim == 4 (RA, Dec, Freq, Stokes), but
        # pixel_n_dim can be 2 if astropy removed unused axes.

        # centre of the spatial pixel grid (0-based)
        center_x = header.get("CRPIX1", 1) - 1
        center_y = header.get("CRPIX2", 1) - 1

        if spatial_pixel_dim == 2:
            sky = wcs.pixel_to_world(center_x, center_y)
        else:
            # Fill extra pixel axes with their reference pixel (0-based)
            extra_pixels = [
                header.get(f"CRPIX{i}", 1) - 1
                for i in range(3, spatial_pixel_dim + 1)
            ]
            sky = wcs.pixel_to_world(center_x, center_y, *extra_pixels)

        if isinstance(sky, SkyCoord):
            return sky
        # multi-dimensional return: extract celestial component
        if isinstance(sky, (tuple, list)):
            for item in sky:
                if isinstance(item, SkyCoord):
                    return item
            # if first element has .ra/.dec (e.g. SkyCoord), use it
            if hasattr(sky[0], "ra") and hasattr(sky[0], "dec"):
                return SkyCoord(ra=sky[0].ra, dec=sky[0].dec, frame="icrs")
            ra, dec = sky[0], sky[1]
            return SkyCoord(ra=ra, dec=dec, frame="icrs")

        # Guard against unexpected types (e.g. plain list from some WCS)
        if hasattr(sky, "ra") and hasattr(sky, "dec"):
            return SkyCoord(ra=sky.ra, dec=sky.dec, frame="icrs")
        raise ValueError(
            f"Unexpected type returned by pixel_to_world: {type(sky)}. "
            f"Expected SkyCoord or tuple. Contents: {sky!r}"
        )

    @staticmethod
    def _to_numpy_array(data_array) -> np.ndarray:
        """Ensure the OSKAR array is a plain NumPy array.

        OSKAR's to_array() returns a 2D array with columns:
            [ra, dec, I, ...]
        depending on its internal encoding.
        """
        if hasattr(data_array, "to_numpy"):
            return data_array.to_numpy()
        return np.asarray(data_array)
