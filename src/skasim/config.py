from __future__ import annotations

from datetime import datetime
from typing import List, Literal, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

# spectral-grid defaults (used when all three are omitted)
_DEFAULT_BW_MHZ = 100.0
_DEFAULT_NCH = 8
_DEFAULT_DF_MHZ = 12.5
CatalogueName = Literal["MIGHTEE", "GLEAM", "SKAMid"]
_CATALOGUE_NAMES = {
    "MIGHTEE": "MIGHTEE",
    "GLEAM": "GLEAM",
    "SKAMID": "SKAMid",
}
_CATALOGUE_MIGRATION_MESSAGE = (
    "Numeric catalogue IDs were removed in skasim 0.2; use named catalogues "
    "such as MIGHTEE, GLEAM, or SKAMid."
)


class ObsConfig(BaseModel):
    """observation parameters for Karabo."""

    freq_mhz: float = Field(700.0, gt=0)
    bandwidth_mhz: Optional[float] = Field(default=None)
    n_channels: Optional[int] = Field(default=None)
    delta_freq_mhz: Optional[float] = Field(default=None)
    seconds: int = Field(600, gt=0)
    phase_center_ra_deg: Optional[float] = None
    phase_center_dec_deg: Optional[float] = None
    start_time: Optional[datetime] = None

    @model_validator(mode="after")
    def _resolve_spectral_grid(self):
        bw = self.bandwidth_mhz
        nch = self.n_channels
        df = self.delta_freq_mhz
        defined = sum(v is not None for v in (bw, nch, df))

        if defined == 0:
            # all omitted — apply defaults
            self.bandwidth_mhz = _DEFAULT_BW_MHZ
            self.n_channels = _DEFAULT_NCH
            self.delta_freq_mhz = _DEFAULT_DF_MHZ
            return self

        if defined == 1:
            raise ValueError(
                "at least two of bandwidth_mhz, n_channels, delta_freq_mhz are required"
            )

        if defined == 3:
            # all provided — tolerate if mathematically consistent
            if abs(bw - nch * df) > 1e-6:
                raise ValueError(
                    f"inconsistent grid: {bw=} ≠ {nch} × {df} = {nch * df}"
                )
            return self

        # two provided, derive the third
        if bw is None:
            self.bandwidth_mhz = nch * df
        elif nch is None:
            self.n_channels = max(1, round(bw / df))
            self.bandwidth_mhz = self.n_channels * df
        else:  # df is None
            self.delta_freq_mhz = bw / nch

        return self


class ImgConfig(BaseModel):
    """imaging parameters passed to OSKAR / WSClean."""

    pixels: int = 512
    fov_deg: Optional[float] = None
    imaging_niter: int = 1000
    robust: float = 0.0
    imager: Literal["oskar-dirty", "wsclean"] = "oskar-dirty"
    wsclean_command: str = "wsclean"
    algorithm: Literal["oskar_dirty", "wsclean_clean"] = "oskar_dirty"

    @field_validator("pixels")
    @classmethod
    def _power_of_two(cls, v: int) -> int:
        if v < 64:
            raise ValueError("pixels must be >= 64")
        return v


class SimConfig(BaseModel):
    """simulation settings."""

    telescope: str = "SKA1MID"
    telescope_version: Optional[str] = None

    # sky input (pipeline resolves one explicit source, else generated sources)
    sky_file: Optional[str] = None
    sky_format: Literal["auto", "fits", "json", "pickle", "random"] = "auto"
    catalogue: Optional[CatalogueName] = None
    column_mapping: Optional[str] = "0,1,2,3,4,5,6,7,8,9,10,11,12"
    scale_I: float = 1.0

    # inline / random source generation
    source_names: Optional[List[str]] = None
    source_intensities: Optional[List[float]] = None
    I: List[float] = Field(default=[10.0])
    Q: Optional[float] = None
    U: Optional[float] = None
    V: Optional[float] = None
    ref_freq_hz: Optional[List[float]] = None

    # foreground json
    json_fg: Optional[str] = None

    # field center string
    center: Optional[str] = None

    # noise / rms
    rms: bool = False
    rms_value: float = 0.0
    rms_sigma: float = 3.0

    # wsclean iterations
    niter: int = 5000

    # nested configs
    observation: ObsConfig = ObsConfig()
    imaging: ImgConfig = ImgConfig()

    output_prefix: Optional[str] = None
    overwrite: bool = False
    cleaning: bool = False

    @field_validator("catalogue", mode="before")
    @classmethod
    def _normalise_catalogue(cls, value):
        if value in (None, "", 0):
            return None
        if isinstance(value, int) or (isinstance(value, str) and value.isdigit()):
            raise ValueError(_CATALOGUE_MIGRATION_MESSAGE)
        if isinstance(value, str):
            key = value.upper()
            if key in _CATALOGUE_NAMES:
                return _CATALOGUE_NAMES[key]
        return value

    @model_validator(mode="after")
    def _validate_one_sky_model_source(self):
        explicit_sources = [
            source
            for source in (self.sky_file, self.catalogue)
            if source is not None
        ]
        if len(explicit_sources) > 1:
            raise ValueError(
                "Provide one sky model source per run; choose a file-backed "
                "sky model or a named catalogue."
            )
        if self.source_intensities is not None and explicit_sources:
            raise ValueError(
                "Source intensity flags are only valid in generated source mode."
            )
        return self
