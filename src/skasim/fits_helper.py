"""fits_helper.py — fuzzy column-name resolver for FITS catalogues.

Provides auto-detection of relevant hcolumns (RA, Dec, flux, etc.) from arbitrary
FITS table names by exact-, substring- and fuzzy-matching against known aliases.
"""

from __future__ import annotations

import difflib
from typing import Dict, List, Optional

import astropy.units as u
import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.units import UnitBase
from .sky import SkyModel, Source

from loguru import logger

from .runtime import require_karabo_module
from .utils import mapping_unit

''' TODO: Review & simplify.
# known aliases for each FITS field
COLUMN_ALIASES: Dict[str, List[str]] = {
    "id": [
        "id", "source_id", "objid", "name", "sourceid", "catalog_id",
    ],
    "ra": [
        "ra", "raj2000", "ra_deg", "right_ascension", "ra_equ", "ra_decimal",
        "alphaj2000", "alpha", "ra2000",
    ],
    "dec": [
        "dec", "dej2000", "dec_deg", "declination", "dec_equ", "dec_decimal",
        "deltaj2000", "delta", "dec2000",
    ],
    "stokes_i": [
        "i", "stokes_i", "flux", "peak_flux", "total_flux", "int_flux",
        "flux_peak", "flux_total", "integrated_flux", "stokesi", "si",
        "f", "f_int", "f_peak", "s",
    ],
    "stokes_q": [
        "q", "stokes_q", "stokesq", "sq",
    ],
    "stokes_u": [
        "u", "stokes_u", "stokesu", "su",
    ],
    "stokes_v": [
        "v", "stokes_v", "stokesv", "sv",
    ],
    "spectral_index": [
        "alpha", "spectral_index", "si", "spidx", "spindex", "spi",
        "spectralindex", "spec_index",
    ],
    "ref_freq": [
        "ref_freq", "frequency", "freq", "ref_frequency", "eff_freq",
        "effective_freq", "freq_peak", "centre_freq", "center_freq",
        "frequency_hz", "frequency_mhz", "freq_hz", "freq_mhz",
    ],
    "rm": [
        "rm", "rotation_measure", "rot_meas", "rotationmeasure",
    ],
    "major": [
        "major", "major_axis", "maj", "dc_maj", "deconvolved_major",
        "maj_axis", "majaxis",
    ],
    "minor": [
        "minor", "minor_axis", "min", "dc_min", "deconvolved_minor",
        "min_axis", "minaxis",
    ],
    "pa": [
        "pa", "position_angle", "pos_angle", "dc_pa", "deconvolved_pa",
        "positionangle",
    ],
}


# unit defaults when the FITS header does not provide them
DEFAULT_UNIT_BY_NAME: Dict[str, UnitBase] = {
    "ra": u.deg,
    "dec": u.deg,
    "stokes_i": u.Jy,
    "stokes_q": u.Jy,
    "stokes_u": u.Jy,
    "stokes_v": u.Jy,
    "spectral_index": u.dimensionless_unscaled,
    "ref_freq": u.MHz,
    "rm": u.rad / u.m**2,
    "major": u.arcsec,
    "minor": u.arcsec,
    "pa": u.deg,
}


def _normalise(name: str) -> str:
    """lower-case, strip whitespace, remove underscores and hyphens."""
    return name.lower().strip().replace("_", "").replace("-", "").replace(" ", "")


def _best_match(
    candidates: List[str],
    aliases: List[str],
    cutoff: float = 0.6,
) -> Optional[str]:
    """Return the candidate that best matches any alias.

    Scoring:
        1. exact normalised match (highest priority)
        2. substring containment (alias in candidate or vice versa)
        3. difflib fuzzy match with cutoff
    """
    norm_candidates = {_normalise(c): c for c in candidates}
    norm_aliases = [_normalise(a) for a in aliases]

    # 1. exact match
    for na, alias in zip(norm_aliases, aliases):
        if na in norm_candidates:
            return norm_candidates[na]

    # 2. substring
    scores: Dict[str, float] = {}
    for nc, orig in norm_candidates.items():
        for na in norm_aliases:
            if na in nc or nc in na:
                # prefer shorter candidate (less chance of confusion with E_* variants)
                scores[orig] = max(scores.get(orig, 0.0), 1.0 - len(orig) * 0.001)
        if orig not in scores:
            # 3. fuzzy
            matches = difflib.get_close_matches(nc, norm_aliases, n=1, cutoff=cutoff)
            if matches:
                scores[orig] = difflib.SequenceMatcher(None, nc, matches[0]).ratio()

    if not scores:
        return None

    # pick highest score; tie-break by shorter original name
    best = sorted(scores.items(), key=lambda x: (-x[1], len(x[0])))[0][0]
    return best


class FitsColumnResolver:
    """Resolve FITS table column names to semantic roles used by Karabo SkyModel."""

    def __init__(self, column_names: List[str], cutoff: float = 0.6) -> None:
        self.column_names = list(column_names)
        self.cutoff = cutoff
        self.resolved: Dict[str, Optional[str]] = {}
        self._resolve_all()

    # ------------------------------------------------------------------
    # internal
    # ------------------------------------------------------------------

    def _resolve_all(self) -> None:
        matched: set[str] = set()
        for field, aliases in COLUMN_ALIASES.items():
            best = _best_match(self.column_names, aliases, self.cutoff)
            if best is not None:
                matched.add(best)
            self.resolved[field] = best

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------

    def get(self, field: str) -> Optional[str]:
        return self.resolved.get(field)

    def get_unit(self, column_name: str) -> UnitBase:
        """Return a sensible astropy Unit for *column_name* based on its role."""
        # reverse lookup: which field maps to this column?
        for field, matched in self.resolved.items():
            if matched == column_name:
                return DEFAULT_UNIT_BY_NAME.get(field, u.dimensionless_unscaled)
        return u.dimensionless_unscaled

    @property
    def prefix_mapping(self) -> SkyPrefixMapping:
        return SkyPrefixMapping(
            id=self.get("id"),
            ra=self.resolved["ra"],            # required — raises downstream if None
            dec=self.resolved["dec"],          # required
            stokes_i=self.resolved["stokes_i"],# required
            stokes_q=self.get("stokes_q"),
            stokes_u=self.get("stokes_u"),
            stokes_v=self.get("stokes_v"),
            spectral_index=self.get("spectral_index"),
            ref_freq=self.get("ref_freq"),
            rm=self.get("rm"),
            major=self.get("major"),
            minor=self.get("minor"),
            pa=self.get("pa"),
        )

    @property
    def unit_mapping(self) -> Dict[str, UnitBase]:
        """Mapping of *column name* -> astropy Unit, for every column in the table."""
        mapping: Dict[str, UnitBase] = {}
        for name in self.column_names:
            mapping[name] = self.get_unit(name)
        return mapping

    @property
    def units_sources(self) -> SkySourcesUnits:
        """SkySourcesUnits for the Karabo loader."""
        # try with /beam first; caller retries without /beam on UnitConversionError
        return SkySourcesUnits(
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


def build_prefix_and_unit_mapping(
    column_names: List[str],
    cutoff: float = 0.6,
) -> tuple[SkyPrefixMapping, Dict[str, UnitBase], SkySourcesUnits]:
    """Convenience factory — returns (prefix_mapping, unit_mapping, units_sources)."""
    resolver = FitsColumnResolver(column_names, cutoff=cutoff)
    return resolver.prefix_mapping, resolver.unit_mapping, resolver.units_sources
'''

# ---------------------------------------------------------------------------
# FitsCatalogLoader — load SkyModel from a FITS table with explicit column mapping.
# ---------------------------------------------------------------------------

class FitsCatalogLoader:
    """Load a `SkyModel` from a FITS table using a `column_mapping` string.

    Automatically falls back to a custom loader (no Karabo) when mapped columns lack TUNIT in the FITS header.
    When all mapped columnsprovide TUNIT, the path via `SkyPrefixMapping` + `SkySourcesUnits` is used.
    """

    # default astropy units for each position in column_mapping
    _UNIT_BY_POS: Dict[int, UnitBase] = {
        1: u.deg,       # ra
        2: u.deg,       # dec
        3: u.Jy,        # I
        4: u.Jy,        # Q
        5: u.Jy,        # U
        6: u.Jy,        # V
        7: u.dimensionless_unscaled,  # spectral_index
        8: u.MHz,       # ref_freq
        9: u.rad / u.m**2,  # rot_meas
        10: u.arcsec,   # major_axis
        11: u.arcsec,   # minor_axis
        12: u.deg,      # pa
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
                    unit_mapping[col.unit] = u.Unit(mapped) if mapped else u.dimensionless_unscaled

            prefix_mapping = sky_model_module.SkyPrefixMapping(
                ra=hdu1.columns.names[self.cols_mapping[1]],
                dec=hdu1.columns.names[self.cols_mapping[2]],
                stokes_i=hdu1.columns.names[self.cols_mapping[3]],
                stokes_q=hdu1.columns.names[self.cols_mapping[4]]
                if self.cols_mapping[4] > -1 else None,
                stokes_u=hdu1.columns.names[self.cols_mapping[5]]
                if self.cols_mapping[5] > -1 else None,
                stokes_v=hdu1.columns.names[self.cols_mapping[6]]
                if self.cols_mapping[6] > -1 else None,
                spectral_index=hdu1.columns.names[self.cols_mapping[7]]
                if self.cols_mapping[7] > -1 else None,
                ref_freq=hdu1.columns.names[self.cols_mapping[8]]
                if self.cols_mapping[8] > -1 else None,
                rm=hdu1.columns.names[self.cols_mapping[9]]
                if self.cols_mapping[9] > -1 else None,
                major=hdu1.columns.names[self.cols_mapping[10]]
                if self.cols_mapping[10] > -1 else None,
                minor=hdu1.columns.names[self.cols_mapping[11]]
                if self.cols_mapping[11] > -1 else None,
                pa=hdu1.columns.names[self.cols_mapping[12]]
                if self.cols_mapping[12] > -1 else None,
                id=hdu1.columns.names[self.cols_mapping[0]]
                if self.cols_mapping[0] > -1 else None,
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
