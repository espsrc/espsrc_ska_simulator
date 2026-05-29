"""tests/test_config.py"""

import pytest
from pydantic import ValidationError

from skasim.config import ImgConfig, ObsConfig, SimConfig

# ---------------------------------------------------------------------------
# ObsConfig
# ---------------------------------------------------------------------------


def test_obs_config_defaults():
    """all three spectral-grid parameters omitted -> defaults resolved."""
    obs = ObsConfig()
    assert obs.frequency_mhz == 700.0
    assert obs.bandwidth_mhz == 100.0  # default resolved
    assert obs.n_channels == 8  # default resolved
    assert obs.channel_width_mhz == 12.5
    assert obs.observation_time_s == 600
    assert obs.phase_center_ra_deg is None
    assert obs.phase_center_dec_deg is None
    assert obs.start_time is None


def test_obs_config_custom_values():
    """custom values are stored correctly."""
    obs = ObsConfig(
        frequency_mhz=1420.0,
        bandwidth_mhz=200.0,
        n_channels=16,
        observation_time_s=3600,
        phase_center_ra_deg=150.0,
        phase_center_dec_deg=2.5,
    )
    assert obs.frequency_mhz == pytest.approx(1420.0)
    assert obs.bandwidth_mhz == pytest.approx(200.0)
    assert obs.n_channels == 16
    assert obs.channel_width_mhz == pytest.approx(12.5)  # derived
    assert obs.observation_time_s == 3600
    assert obs.phase_center_ra_deg == pytest.approx(150.0)
    assert obs.phase_center_dec_deg == pytest.approx(2.5)


@pytest.mark.parametrize("bad_observation_time_s", [0, -1, -100])
def test_obs_config_non_positive_observation_time_s_raises(bad_observation_time_s):
    """observation_time_s must be strictly positive."""
    with pytest.raises(ValidationError):
        ObsConfig(observation_time_s=bad_observation_time_s)


# ---------------------------------------------------------------------------
# consistent_setup
# ---------------------------------------------------------------------------


def test_consistent_setup_bw_and_nch():
    """bandwidth + n_channels -> delta_freq derived."""
    obs = ObsConfig(bandwidth_mhz=100.0, n_channels=8)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.channel_width_mhz == pytest.approx(12.5)


def test_consistent_setup_bw_and_df():
    """bandwidth + channel width -> n_channels derived (rounded)."""
    obs = ObsConfig(bandwidth_mhz=100.0, channel_width_mhz=12.5)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.channel_width_mhz == pytest.approx(12.5)


def test_consistent_setup_nch_and_df():
    """no bandwidth -> bandwidth derived from n_channels * channel width."""
    obs = ObsConfig(n_channels=8, channel_width_mhz=12.5)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.channel_width_mhz == pytest.approx(12.5)


# ---------------------------------------------------------------------------
# spectral-grid validator — invalid combinations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_msg",
    [
        # all three given but inconsistent
        (
            dict(bandwidth_mhz=100.0, n_channels=8, channel_width_mhz=10.0),
            "inconsistent grid",
        ),
        # only one given — bandwidth
        (
            dict(bandwidth_mhz=100.0),
            "at least two of bandwidth_mhz, n_channels, channel_width_mhz",
        ),
        # only one given — n_channels
        (
            dict(n_channels=8),
            "at least two of bandwidth_mhz, n_channels, channel_width_mhz",
        ),
        # only one given — delta_freq
        (
            dict(channel_width_mhz=12.5),
            "at least two of bandwidth_mhz, n_channels, channel_width_mhz",
        ),
    ],
)
def test_obs_config_invalid_spectral_grid(kwargs, expected_msg):
    """model_validator rejects inconsistent or under-specified combos."""
    with pytest.raises(ValidationError, match=expected_msg):
        ObsConfig(**kwargs)


# ---------------------------------------------------------------------------
# ImgConfig
# ---------------------------------------------------------------------------


def test_img_config_defaults():
    """default imaging parameters."""
    img = ImgConfig()
    assert img.pixels == 512
    assert img.fov_deg is None
    assert img.robust == 0.0
    assert img.imager == "oskar-dirty"
    assert img.wsclean_command == "wsclean"


def test_img_config_wsclean_imager():
    """wsclean imager selection is explicit."""
    img = ImgConfig(imager="wsclean", pixels=1024)
    assert img.imager == "wsclean"
    assert img.pixels == 1024


def test_img_config_rejects_removed_algorithm_field():
    """algorithm was removed as stored config state."""
    with pytest.raises(ValidationError):
        ImgConfig(algorithm="wsclean_clean")


@pytest.mark.parametrize("bad_pixels", [32, 63, 0, -128])
def test_img_config_too_small_pixels_raises(bad_pixels):
    """image size valid size."""
    with pytest.raises(ValidationError, match="pixels must be >= 64"):
        ImgConfig(pixels=bad_pixels)


# ---------------------------------------------------------------------------
# SimConfig
# ---------------------------------------------------------------------------


def test_sim_config_defaults():
    """default settings"""
    cfg = SimConfig()
    assert cfg.telescope == "SKA1MID"
    assert cfg.telescope_version is None
    assert cfg.sky_file is None
    assert cfg.sky_format == "auto"
    assert cfg.catalog is None
    assert cfg.column_mapping == "0,1,2,3,4,5,6,7,8,9,10,11,12"
    assert cfg.flux_scale == 1.0
    assert cfg.source_flux_jy == [10.0]
    assert cfg.stokes_q_jy is None
    assert cfg.stokes_u_jy is None
    assert cfg.stokes_v_jy is None
    assert cfg.center is None
    assert cfg.rms is False
    assert cfg.rms_value == 0.0
    assert cfg.rms_sigma == 3.0
    assert cfg.output_dir is None
    assert cfg.overwrite is False
    # nested configurations (post refactor)
    assert isinstance(cfg.observation, ObsConfig)
    assert isinstance(cfg.imaging, list)
    assert cfg.imaging[0].pixels == 512
    assert cfg.imaging[0].clean_iterations == 5000
    assert cfg.imaging[0].tag == "default"
    assert cfg.observation.frequency_mhz == 700.0


def test_sim_config_explicit_telescope_version():
    """telescope with version."""
    cfg = SimConfig(telescope="SKA1LOW", telescope_version="AA0.5")
    assert cfg.telescope == "SKA1LOW"
    assert cfg.telescope_version == "AA0.5"


def test_sim_config_accepts_named_catalog():
    """built-in catalogs are selected by name."""
    cfg = SimConfig(catalog="MIGHTEE")
    assert cfg.catalog == "MIGHTEE"
    assert cfg.source_flux_jy == []
    assert cfg.stokes_q_jy is None
    assert cfg.stokes_u_jy is None
    assert cfg.stokes_v_jy is None


def test_sim_config_accepts_serialized_catalog_run_empty_source_flux():
    """Manifest round-trips keep explicit-source runs loadable."""
    cfg = SimConfig(catalog="MIGHTEE", source_flux_jy=[])

    assert cfg.catalog == "MIGHTEE"
    assert cfg.source_flux_jy == []


def test_sim_config_numeric_catalog_zero_has_migration_message():
    """catalog=0 triggers the numeric migration error, not silent None."""
    with pytest.raises(ValidationError, match="Numeric catalog IDs were removed"):
        SimConfig(catalog=0)


def test_sim_config_numeric_catalog_has_migration_message():
    """numeric catalog IDs are removed in 0.2."""
    with pytest.raises(ValidationError, match="Numeric catalog IDs were removed"):
        SimConfig(catalog=1)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"I": [1.0]},
        {"Q": [0.1]},
        {"U": [0.1]},
        {"V": [0.1]},
        {"cleaning": True},
        {"source_names": ["a"]},
        {"ref_freq_hz": [1.4e9]},
        {"json_fg": "fg.json"},
        {"output_prefix": "old"},
        {"niter": 10},
        {"scale_I": 2.0},
    ],
)
def test_sim_config_rejects_removed_fields(kwargs):
    """Deprecated 0.1 config fields are not accepted in the 0.2 access layer."""
    with pytest.raises(ValidationError):
        SimConfig(**kwargs)


def test_sim_config_rejects_model_and_catalog_together():
    """0.2 accepts one explicit sky model source per run."""
    with pytest.raises(ValidationError, match="one sky model source"):
        SimConfig(sky_file="sources.json", catalog="GLEAM")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"sky_file": "sources.json"},
        {"catalog": "GLEAM"},
    ],
)
def test_sim_config_rejects_source_flux_jy_with_explicit_source(kwargs):
    """Source flux flags are only valid in generated source mode."""
    with pytest.raises(ValidationError, match="generated source mode"):
        SimConfig(source_flux_jy=[1.0], **kwargs)


def test_sim_config_rejects_stokes_values_with_explicit_source():
    """Generated-source polarization flags cannot be combined with catalogs."""
    with pytest.raises(ValidationError, match="generated source mode"):
        SimConfig(catalog="MIGHTEE", stokes_q_jy=[0.1])


def test_sim_config_rejects_empty_generated_intensities():
    """Generated source mode needs at least one flux density."""
    with pytest.raises(ValidationError, match="at least one flux density"):
        SimConfig(source_flux_jy=[])


def test_sim_config_accepts_generated_source_polarization():
    """Generated-source Stokes Q/U/V values match source fluxes by position."""
    cfg = SimConfig(
        source_flux_jy=[1.0, 2.0],
        stokes_q_jy=[0.1, 0.2],
        stokes_u_jy=[0.0, 0.1],
        stokes_v_jy=[0.0, -0.1],
    )

    assert cfg.stokes_q_jy == [0.1, 0.2]
    assert cfg.stokes_u_jy == [0.0, 0.1]
    assert cfg.stokes_v_jy == [0.0, -0.1]


def test_sim_config_rejects_stokes_length_mismatch():
    """Generated-source Stokes vectors must align with the flux-density list."""
    with pytest.raises(ValidationError, match="stokes_q_jy must contain 2 values"):
        SimConfig(source_flux_jy=[1.0, 2.0], stokes_q_jy=[0.1])


def test_sim_config_nested_observation_override():
    """nested ObsConfig can be fully replaced."""
    cfg = SimConfig(observation=ObsConfig(observation_time_s=1200, frequency_mhz=900.0))
    assert cfg.observation.observation_time_s == 1200
    assert cfg.observation.frequency_mhz == pytest.approx(900.0)
    # other nested defaults remain untouched
    assert cfg.imaging[0].pixels == 512


def test_sim_config_nested_imaging_override():
    """nested ImgConfig can be fully replaced."""
    cfg = SimConfig(imaging=[ImgConfig(tag="main", pixels=2048, imager="wsclean")])
    assert cfg.imaging[0].pixels == 2048
    assert cfg.imaging[0].imager == "wsclean"
    # other nested defaults remain untouched
    assert cfg.observation.observation_time_s == 600


def test_sim_config_wraps_single_imaging_dict():
    """Backward compat: a single imaging dict becomes a list of one."""
    cfg = SimConfig(imaging={"tag": "legacy", "pixels": 256})
    assert isinstance(cfg.imaging, list)
    assert len(cfg.imaging) == 1
    assert cfg.imaging[0].tag == "legacy"
    assert cfg.imaging[0].pixels == 256


def test_sim_config_rejects_duplicate_tags():
    """Duplicate imaging tags are not allowed."""
    with pytest.raises(ValidationError, match="duplicate imaging tags"):
        SimConfig(
            imaging=[
                ImgConfig(tag="a", pixels=256),
                ImgConfig(tag="a", pixels=512),
            ]
        )


def test_sim_config_rejects_empty_imaging_list():
    """At least one imaging block is required."""
    with pytest.raises(ValidationError, match="at least one imaging block"):
        SimConfig(imaging=[])


def test_img_config_tag_validation():
    """Tag must not be empty or contain path-special chars."""
    with pytest.raises(ValidationError, match="tag must not be empty"):
        ImgConfig(tag="")
    with pytest.raises(ValidationError, match="tag must not contain whitespace"):
        ImgConfig(tag="foo bar")
    with pytest.raises(ValidationError, match="tag must not contain whitespace"):
        ImgConfig(tag="foo/bar")
