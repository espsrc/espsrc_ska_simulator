from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Annotated, List, Literal, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, field_validator, model_validator

# spectral-grid defaults (used when all three are omitted)
_DEFAULT_BW_MHZ = 100.0
_DEFAULT_NCH = 8
_DEFAULT_DF_MHZ = 12.5
CatalogName = Literal["MIGHTEE", "GLEAM", "SKAMid"]
_CATALOG_NAMES = {
    "MIGHTEE": "MIGHTEE",
    "GLEAM": "GLEAM",
    "SKAMID": "SKAMid",
}
_CATALOG_MIGRATION_MESSAGE = (
    "Numeric catalog IDs were removed in skasim 0.2; use named catalogs "
    "such as MIGHTEE, GLEAM, or SKAMid."
)


def _normalise_catalog_value(value):
    if value is None or value == "":
        return None
    if isinstance(value, int) or (isinstance(value, str) and value.isdigit()):
        raise ValueError(_CATALOG_MIGRATION_MESSAGE)
    if isinstance(value, str):
        key = value.upper()
        if key in _CATALOG_NAMES:
            return _CATALOG_NAMES[key]
    return value


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #


def _require_existing_path(value: str) -> str:
    path = Path(value).expanduser()
    if not path.exists():
        raise ValueError(f"model file does not exist: {value}")
    return value


# --------------------------------------------------------------------------- #
# model entry types (typed multi-model API)
# --------------------------------------------------------------------------- #


class ComponentSkyModelEntry(BaseModel):
    """Existing catalog/component sky-model entry."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["component_sky_model"]
    path: Optional[str] = None
    catalog: Optional[CatalogName] = None
    sky_format: Literal["auto", "fits", "json", "pickle", "random"] = "auto"
    column_mapping: Optional[str] = "0,1,2,3,4,5,6,7,8,9,10,11,12"
    flux_scale: float = 1.0

    @field_validator("path")
    @classmethod
    def _path_exists(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return value
        return _require_existing_path(value)

    @field_validator("catalog", mode="before")
    @classmethod
    def _normalise_catalog(cls, value):
        return _normalise_catalog_value(value)

    @model_validator(mode="after")
    def _validate_one_component_source(self):
        defined = sum(value is not None for value in (self.path, self.catalog))
        if defined != 1:
            raise ValueError(
                "component_sky_model requires exactly one of path or catalog."
            )
        return self


class ContinuumIAlphaModelEntry(BaseModel):
    """Continuum image model with a spatially varying spectral index."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["continuum_i_alpha"]
    stokes_i: str
    alpha: str
    reference_frequency_hz: float = Field(gt=0)

    @field_validator("stokes_i", "alpha")
    @classmethod
    def _path_exists(cls, value: str) -> str:
        return _require_existing_path(value)


class CasaTaylorTermsModelEntry(BaseModel):
    """Existing CASA Taylor-term image model set."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["casa_taylor_terms"]
    tt0: str
    tt1: Optional[str] = None
    reference_frequency_hz: float = Field(gt=0)

    @field_validator("tt0", "tt1")
    @classmethod
    def _path_exists(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return value
        return _require_existing_path(value)


class StaticStokesMapsModelEntry(BaseModel):
    """Static Stokes map set. Schema-ready for the phase-2 backend path."""

    model_config = ConfigDict(extra="forbid")

    type: Literal["static_stokes_maps"]
    stokes_i: Optional[str] = None
    stokes_q: Optional[str] = None
    stokes_u: Optional[str] = None
    stokes_v: Optional[str] = None

    @field_validator("stokes_i", "stokes_q", "stokes_u", "stokes_v")
    @classmethod
    def _path_exists(cls, value: Optional[str]) -> Optional[str]:
        if value is None:
            return value
        return _require_existing_path(value)

    @model_validator(mode="after")
    def _at_least_one_stokes_map(self):
        if not any((self.stokes_i, self.stokes_q, self.stokes_u, self.stokes_v)):
            raise ValueError("static_stokes_maps requires at least one Stokes map.")
        return self


ModelEntry = Annotated[
    Union[
        ComponentSkyModelEntry,
        ContinuumIAlphaModelEntry,
        CasaTaylorTermsModelEntry,
        StaticStokesMapsModelEntry,
    ],
    Field(discriminator="type"),
]


# --------------------------------------------------------------------------- #
# ObsConfig
# --------------------------------------------------------------------------- #


class ObsConfig(BaseModel):
    """observation parameters for Karabo."""

    model_config = ConfigDict(extra="forbid")

    frequency_mhz: float = Field(700.0, gt=0)
    bandwidth_mhz: Optional[float] = Field(default=None)
    n_channels: Optional[int] = Field(default=None)
    channel_width_mhz: Optional[float] = Field(default=None)
    observation_time_s: int = Field(600, gt=0)
    phase_center_ra_deg: Optional[float] = None
    phase_center_dec_deg: Optional[float] = None
    start_time: Optional[datetime] = None

    @model_validator(mode="after")
    def _resolve_spectral_grid(self):
        bw = self.bandwidth_mhz
        nch = self.n_channels
        df = self.channel_width_mhz
        defined = sum(v is not None for v in (bw, nch, df))

        if defined == 0:
            # all omitted — apply defaults
            self.bandwidth_mhz = _DEFAULT_BW_MHZ
            self.n_channels = _DEFAULT_NCH
            self.channel_width_mhz = _DEFAULT_DF_MHZ
            return self

        if defined == 1:
            raise ValueError(
                "at least two of bandwidth_mhz, n_channels, channel_width_mhz are required"
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
            self.channel_width_mhz = bw / nch

        return self


# --------------------------------------------------------------------------- #
# ImgConfig
# --------------------------------------------------------------------------- #


class ImgConfig(BaseModel):
    """imaging parameters passed to OSKAR / WSClean.

    WSClean-specific flags (mgain, multiscale, auto-threshold, etc.) are
    ignored when ``imager`` is ``oskar-dirty``; they only affect argv building
    for ``wsclean``.
    """

    model_config = ConfigDict(extra="forbid")

    tag: str = "default"
    pixels: int = 512
    fov_deg: Optional[float] = None
    robust: float = 0.0
    imager: Literal["oskar-dirty", "wsclean"] = "oskar-dirty"
    wsclean_command: str = "wsclean"
    clean_iterations: int = 5000

    # WSClean-only flags.  None means "use the skasim default".
    mgain: Optional[float] = None
    multiscale: Optional[bool] = None
    multiscale_scales: Optional[List[int]] = None
    auto_threshold: Optional[float] = None
    auto_mask: Optional[float] = None
    local_rms: Optional[bool] = None
    join_channels: Optional[bool] = None
    channels_out: Optional[int] = None
    padding: Optional[float] = None
    threads: Optional[int] = None

    @field_validator("pixels")
    @classmethod
    def _min_pixels(cls, v: int) -> int:
        if v < 64:
            raise ValueError("pixels must be >= 64")
        return v

    @field_validator("tag")
    @classmethod
    def _valid_tag(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("tag must not be empty")
        if any(c in v for c in [" ", "/", "\\", ":", "*", "?", '"', "<", ">", "|"]):
            raise ValueError("tag must not contain whitespace or path-special chars")
        return v


# --------------------------------------------------------------------------- #
# SimConfig
# --------------------------------------------------------------------------- #


class SimConfig(BaseModel):
    """simulation settings."""

    model_config = ConfigDict(extra="forbid")

    telescope: str = "SKA1MID"
    telescope_version: Optional[str] = None

    # sky input (pipeline resolves one explicit source, else generated sources)
    models: List[ModelEntry] = Field(default_factory=list)
    sky_file: Optional[str] = None
    sky_format: Literal["auto", "fits", "json", "pickle", "random"] = "auto"
    catalog: Optional[CatalogName] = None
    fits_image: Optional[str] = None
    column_mapping: Optional[str] = "0,1,2,3,4,5,6,7,8,9,10,11,12"
    flux_scale: float = 1.0

    # inline / random source generation
    source_flux_jy: Optional[List[float]] = Field(default=None)
    stokes_q_jy: Optional[List[float]] = None
    stokes_u_jy: Optional[List[float]] = None
    stokes_v_jy: Optional[List[float]] = None

    # field center string
    center: Optional[str] = None

    # noise / rms
    rms: bool = False
    rms_value: float = 0.0
    rms_sigma: float = 3.0
    noise_rms_start: Optional[float] = None
    noise_rms_end: Optional[float] = None

    # nested configs
    observation: ObsConfig = ObsConfig()
    imaging: List[ImgConfig] = Field(default_factory=lambda: [ImgConfig()])

    # UV coverage (shadeMS) identical for any imaging pass
    uv_coverage: bool = True
    shadems_command: str = "shadems"
    uv_coverage_canvas_size: int = 600

    output_dir: Optional[str] = None
    overwrite: bool = False

    # optional run metadata (config-file only; not a CLI argument)
    title: Optional[str] = None
    description: Optional[str] = None

    @model_validator(mode="before")
    @classmethod
    def _wrap_single_imaging(cls, data):
        """Backward compat: wrap a single imaging dict into a list."""
        if isinstance(data, dict) and "imaging" in data:
            img = data["imaging"]
            if isinstance(img, dict):
                if "tag" not in img:
                    img["tag"] = "default"
                data["imaging"] = [img]
        return data

    @model_validator(mode="before")
    @classmethod
    def _reject_generated_intensities_with_explicit_source(cls, data):
        if not isinstance(data, dict):
            return data
        # typed models mode: reject legacy sky-model fields
        if data.get("models") and any(
            field in data and data.get(field) not in (None, [])
            for field in (
                "sky_file",
                "catalog",
                "fits_image",
                "source_flux_jy",
                "stokes_q_jy",
                "stokes_u_jy",
                "stokes_v_jy",
            )
        ):
            raise ValueError(
                "Use typed models without legacy sky_file, catalog, fits_image, or "
                "generated source fields."
            )
        # legacy mode: same check as before
        if any(
            data.get(field) not in (None, [])
            for field in (
                "source_flux_jy",
                "stokes_q_jy",
                "stokes_u_jy",
                "stokes_v_jy",
            )
        ) and (
            data.get("sky_file") is not None
            or data.get("catalog") is not None
            or data.get("fits_image") is not None
        ):
            raise ValueError(
                "Generated source flux and polarization flags are only valid in "
                "generated source mode."
            )
        return data

    @field_validator("catalog", mode="before")
    @classmethod
    def _normalise_catalog(cls, value):
        return _normalise_catalog_value(value)

    @model_validator(mode="after")
    def _validate_one_sky_model_source(self):
        if self.models:
            component_count = sum(
                1
                for model in self.models
                if getattr(model, "type", None) == "component_sky_model"
            )
            if component_count > 1:
                raise ValueError("Provide at most one component_sky_model entry.")
            self.source_flux_jy = []
            self.stokes_q_jy = None
            self.stokes_u_jy = None
            self.stokes_v_jy = None
            return self

        explicit_sources = [
            source
            for source in (self.sky_file, self.catalog, self.fits_image)
            if source is not None
        ]
        if len(explicit_sources) > 1:
            raise ValueError(
                "Provide one sky model source per run; choose a file-backed "
                "sky model, named catalog, or FITS image."
            )
        if explicit_sources:
            self.source_flux_jy = []
            self.stokes_q_jy = None
            self.stokes_u_jy = None
            self.stokes_v_jy = None
        elif self.source_flux_jy is None:
            self.source_flux_jy = [10.0]
        elif not self.source_flux_jy:
            raise ValueError(
                "Generated source mode requires at least one flux density."
            )
        else:
            n_sources = len(self.source_flux_jy)
            for field_name in ("stokes_q_jy", "stokes_u_jy", "stokes_v_jy"):
                values = getattr(self, field_name)
                if values is not None and len(values) != n_sources:
                    raise ValueError(
                        f"{field_name} must contain {n_sources} values to match "
                        "source_flux_jy."
                    )
        return self

    @model_validator(mode="after")
    def _reject_duplicate_tags(self):
        tags = [img.tag for img in self.imaging]
        if len(tags) != len(set(tags)):
            raise ValueError("duplicate imaging tags are not allowed")
        return self

    @model_validator(mode="after")
    def _require_at_least_one_imaging(self):
        if not self.imaging:
            raise ValueError("at least one imaging block is required")
        return self

    @field_validator("uv_coverage_canvas_size")
    @classmethod
    def _positive_uv_coverage_canvas_size(cls, v: int) -> int:
        if v < 1:
            raise ValueError("uv_coverage_canvas_size must be >= 1")
        return v
