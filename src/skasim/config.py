from __future__ import annotations

from datetime import datetime
from typing import List, Literal, Optional

from pydantic import BaseModel, Field, field_validator, model_validator

# spectral-grid defaults (used when all three are omitted)
_DEFAULT_BW_MHZ = 100.0
_DEFAULT_NCH = 8
_DEFAULT_DF_MHZ = 12.5


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

    # sky input (pipeline resolves priority: sky_file > catalogue > random)
    sky_file: Optional[str] = None
    sky_format: Literal["auto", "fits", "json", "pickle", "random"] = "auto"
    catalogue: Literal[0, 1, 2, 3] = 0
    column_mapping: Optional[str] = "0,1,2,3,4,5,6,7,8,9,10,11,12"
    scale_I: float = 1.0

    # inline / random source generation
    source_names: Optional[List[str]] = None
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
