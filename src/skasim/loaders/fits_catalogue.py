"""loaders/fits_catalogue.py — FITS catalogue loading helpers.

Loads FITS table sky models using an explicit column mapping.
"""

from __future__ import annotations

from typing import Dict, List, Optional

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.units import UnitBase
from loguru import logger

from ..runtime import require_karabo_module
from ..sky import SkyModel, Source
from ..utils import mapping_unit


class FitsCatalogLoader:
    """Load a `SkyModel` from a FITS table using a `column_mapping` string.

    Automatically falls back to a custom loader (no Karabo) when mapped columns lack TUNIT in the FITS header.
    When all mapped columnsprovide TUNIT, the path via `SkyPrefixMapping` + `SkySourcesUnits` is used.
    """

    # default astropy units for each position in column_mapping
    _UNIT_BY_POS: Dict[int, UnitBase] = {
        1: u.deg,  # ra
        2: u.deg,  # dec
        3: u.Jy,  # I
        4: u.Jy,  # Q
        5: u.Jy,  # U
        6: u.Jy,  # V
        7: u.dimensionless_unscaled,  # spectral_index
        8: u.MHz,  # ref_freq
        9: u.rad / u.m**2,  # rot_meas
        10: u.arcsec,  # major_axis
        11: u.arcsec,  # minor_axis
        12: u.deg,  # pa
    }

    def __init__(
        self,
        fpath: str,
        column_mapping: str,
        scale_I: float = 1.0,
        ref_freq_hz: Optional[float] = None,
        frequency: Optional["u.Quantity"] = None,
    ) -> None:
        self.fpath = fpath
        self.cols_mapping = [int(i) for i in column_mapping.split(",")]
        self.scale_I = scale_I
        self.ref_freq_hz = ref_freq_hz
        self.frequency = frequency

    # public methods
    # ------------------------------------------------------------------

    def has_missing_unit(self) -> bool:
        """Return True if any mapped column that requires a unit lacks it.

        Only positions 1-12 are checked, with 0 (id) and 7 (spectral_index)
        ignored as dimensionless
        """
        REQUIRED_UNIT_POS = {1, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12}
        with fits.open(self.fpath) as hdul:
            hdu1 = hdul[1]
            for pos, idx in enumerate(self.cols_mapping):
                if pos not in REQUIRED_UNIT_POS or idx < 0:
                    continue
                if hdu1.columns[idx].unit is None:
                    return True
        return False

    def load(self) -> "SkyModel":
        if self.has_missing_unit():
            logger.warning("Some mapped columns lack TUNIT; using own FITS loader.")
            return self._load_own(Source, SkyModel)
        return self._load_karabo()

    # internal, own loader (no Karabo)
    def _load_own(self, Source: type, SkyModel: type) -> "SkyModel":
        table = Table.read(self.fpath)
        sources = []
        for row in table:
            kwargs = self._row_to_kwargs(row, table.colnames)
            src = Source(**kwargs)
            if self.scale_I != 1.0:
                src.I *= self.scale_I
            if src.ref_freq == 0 and self.ref_freq_hz is not None:
                src.ref_freq = self.ref_freq_hz * u.Hz
            if src.ref_freq == 0 and self.frequency is not None:
                src.ref_freq = self.frequency.to(u.Hz)
            sources.append(src)

        if not sources:
            raise ValueError(f"No sources could be built from {self.fpath}")
        sky_array = np.array([s.to_sky_model() for s in sources])
        sky_model = SkyModel(sky_array)
        sky_model.get_center()
        logger.info(f"Loaded {len(sources)} sources via fallback loader: {self.fpath}")
        return sky_model

    def _row_to_kwargs(self, row, colnames: List[str]) -> Dict[str, object]:
        kwargs: Dict[str, object] = {}
        # mandatory
        for pos, canon in [(1, "ra"), (2, "dec"), (3, "I")]:
            col = colnames[self.cols_mapping[pos]]
            kwargs[canon] = float(row[col]) * self._UNIT_BY_POS[pos]
        # optional
        opt = {
            4: "Q",
            5: "U",
            6: "V",
            7: "spec_index",
            8: "ref_freq",
            9: "rot_meas",
            10: "major_axis",
            11: "minor_axis",
            12: "pa",
        }
        for pos, canon in opt.items():
            if pos < len(self.cols_mapping) and self.cols_mapping[pos] > -1:
                col = colnames[self.cols_mapping[pos]]
                unit = self._UNIT_BY_POS.get(pos, u.dimensionless_unscaled)
                kwargs[canon] = float(row[col]) * unit
        return kwargs

    # Karabo loader (all columns have unit)
    def _load_karabo(self) -> "SkyModel":
        sky_model_module = require_karabo_module("karabo.simulation.sky_model")
        with fits.open(self.fpath) as hdul:
            hdu1 = hdul[1]
            unit_mapping: Dict[str, UnitBase] = {}
            for idx in self.cols_mapping:
                if idx < 0:
                    continue
                col = hdu1.columns[idx]
                if col.unit is not None:
                    mapped = mapping_unit(col.unit)
                    unit_mapping[col.unit] = (
                        u.Unit(mapped) if mapped else u.dimensionless_unscaled
                    )

            prefix_mapping = sky_model_module.SkyPrefixMapping(
                ra=hdu1.columns.names[self.cols_mapping[1]],
                dec=hdu1.columns.names[self.cols_mapping[2]],
                stokes_i=hdu1.columns.names[self.cols_mapping[3]],
                stokes_q=hdu1.columns.names[self.cols_mapping[4]]
                if self.cols_mapping[4] > -1
                else None,
                stokes_u=hdu1.columns.names[self.cols_mapping[5]]
                if self.cols_mapping[5] > -1
                else None,
                stokes_v=hdu1.columns.names[self.cols_mapping[6]]
                if self.cols_mapping[6] > -1
                else None,
                spectral_index=hdu1.columns.names[self.cols_mapping[7]]
                if self.cols_mapping[7] > -1
                else None,
                ref_freq=hdu1.columns.names[self.cols_mapping[8]]
                if self.cols_mapping[8] > -1
                else None,
                rm=hdu1.columns.names[self.cols_mapping[9]]
                if self.cols_mapping[9] > -1
                else None,
                major=hdu1.columns.names[self.cols_mapping[10]]
                if self.cols_mapping[10] > -1
                else None,
                minor=hdu1.columns.names[self.cols_mapping[11]]
                if self.cols_mapping[11] > -1
                else None,
                pa=hdu1.columns.names[self.cols_mapping[12]]
                if self.cols_mapping[12] > -1
                else None,
                id=hdu1.columns.names[self.cols_mapping[0]]
                if self.cols_mapping[0] > -1
                else None,
            )

        units_sources = sky_model_module.SkySourcesUnits(
            stokes_i=u.Jy / u.beam,
            stokes_q=u.Jy / u.beam,
            stokes_u=u.Jy / u.beam,
            stokes_v=u.Jy / u.beam,
            ref_freq=u.MHz,
            major=u.arcsec,
            minor=u.arcsec,
            pa=u.deg,
            rm=u.rad / u.m**2,
        )

        try:
            sky_model = SkyModel.get_sky_model_from_fits(
                fits_file=self.fpath,
                prefix_mapping=prefix_mapping,
                unit_mapping=unit_mapping,
                units_sources=units_sources,
                min_freq=None,
                max_freq=None,
                encoded_freq=None,
                memmap=False,
            )
            logger.info(f"Loaded FITS via Karabo: {self.fpath}")
            return sky_model
        except u.core.UnitConversionError as exc:
            logger.error(f"Beam-unit conversion failed ({exc}); retrying without beam.")

        units_sources = sky_model_module.SkySourcesUnits(
            stokes_i=u.Jy,
            stokes_q=u.Jy,
            stokes_u=u.Jy,
            stokes_v=u.Jy,
            ref_freq=u.MHz,
            major=u.arcsec,
            minor=u.arcsec,
            pa=u.deg,
            rm=u.rad / u.m**2,
        )
        sky_model = SkyModel.get_sky_model_from_fits(
            fits_file=self.fpath,
            prefix_mapping=prefix_mapping,
            unit_mapping=unit_mapping,
            units_sources=units_sources,
            min_freq=None,
            max_freq=None,
            encoded_freq=None,
            memmap=False,
        )
        logger.info(f"Loaded FITS via Karabo (no-beam retry): {self.fpath}")
        return sky_model
