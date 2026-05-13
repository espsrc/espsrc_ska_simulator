from __future__ import annotations

from datetime import datetime
from typing import List, Literal, Optional

from pydantic import BaseModel, Field, field_validator


class ObsConfig(BaseModel):
    """observation parameters for Karabo."""

    # TODO: review defaults
    freq_mhz: float = 700.0
    bandwidth_mhz: float = 100.0
    n_channels: int = 8
    delta_freq_mhz: Optional[float] = None
    seconds: int = 600
    phase_center_ra_deg: Optional[float] = None
    phase_center_dec_deg: Optional[float] = None
    start_time: Optional[datetime] = None

    # just an example, TODO: add more validations
    @field_validator("seconds")
    @classmethod
    def _positive_time(cls, v: int) -> int:
        if v <= 0:
            raise ValueError("Observation time must be > 0 seconds")
        return v


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
