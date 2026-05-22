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
    assert obs.freq_mhz == 700.0
    assert obs.bandwidth_mhz == 100.0  # default resolved
    assert obs.n_channels == 8  # default resolved
    assert obs.delta_freq_mhz == 12.5
    assert obs.seconds == 600
    assert obs.phase_center_ra_deg is None
    assert obs.phase_center_dec_deg is None
    assert obs.start_time is None


def test_obs_config_custom_values():
    """custom values are stored correctly."""
    obs = ObsConfig(
        freq_mhz=1420.0,
        bandwidth_mhz=200.0,
        n_channels=16,
        seconds=3600,
        phase_center_ra_deg=150.0,
        phase_center_dec_deg=2.5,
    )
    assert obs.freq_mhz == pytest.approx(1420.0)
    assert obs.bandwidth_mhz == pytest.approx(200.0)
    assert obs.n_channels == 16
    assert obs.delta_freq_mhz == pytest.approx(12.5)  # derived
    assert obs.seconds == 3600
    assert obs.phase_center_ra_deg == pytest.approx(150.0)
    assert obs.phase_center_dec_deg == pytest.approx(2.5)


@pytest.mark.parametrize("bad_seconds", [0, -1, -100])
def test_obs_config_non_positive_seconds_raises(bad_seconds):
    """seconds must be strictly positive."""
    with pytest.raises(ValidationError):
        ObsConfig(seconds=bad_seconds)


# ---------------------------------------------------------------------------
# consistent_setup
# ---------------------------------------------------------------------------


def test_consistent_setup_bw_and_nch():
    """bandwidth + n_channels -> delta_freq derived."""
    obs = ObsConfig(bandwidth_mhz=100.0, n_channels=8)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.delta_freq_mhz == pytest.approx(12.5)


def test_consistent_setup_bw_and_df():
    """bandwidth + delta_freq -> n_channels derived (rounded)."""
    obs = ObsConfig(bandwidth_mhz=100.0, delta_freq_mhz=12.5)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.delta_freq_mhz == pytest.approx(12.5)


def test_consistent_setup_nch_and_df():
    """no bandwidth -> bandwidth derived from n_channels * delta_freq."""
    obs = ObsConfig(n_channels=8, delta_freq_mhz=12.5)
    assert obs.bandwidth_mhz == pytest.approx(100.0)
    assert obs.n_channels == 8
    assert obs.delta_freq_mhz == pytest.approx(12.5)


# ---------------------------------------------------------------------------
# spectral-grid validator — invalid combinations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kwargs, expected_msg",
    [
        # all three given but inconsistent
        (
            dict(bandwidth_mhz=100.0, n_channels=8, delta_freq_mhz=10.0),
            "inconsistent grid",
        ),
        # only one given — bandwidth
        (
            dict(bandwidth_mhz=100.0),
            "at least two of bandwidth_mhz, n_channels, delta_freq_mhz",
        ),
        # only one given — n_channels
        (
            dict(n_channels=8),
            "at least two of bandwidth_mhz, n_channels, delta_freq_mhz",
        ),
        # only one given — delta_freq
        (
            dict(delta_freq_mhz=12.5),
            "at least two of bandwidth_mhz, n_channels, delta_freq_mhz",
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
    assert img.imaging_niter == 1000
    assert img.robust == 0.0
    assert img.algorithm == "oskar_dirty"


def test_img_config_custom_algorithm():
    """imaging initialization."""
    img = ImgConfig(algorithm="wsclean_clean", pixels=1024)
    assert img.algorithm == "wsclean_clean"
    assert img.pixels == 1024


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
    assert cfg.catalogue is None
    assert cfg.column_mapping == "0,1,2,3,4,5,6,7,8,9,10,11,12"
    assert cfg.scale_I == 1.0
    assert cfg.I == [10.0]
    assert cfg.Q is None
    assert cfg.U is None
    assert cfg.V is None
    assert cfg.ref_freq_hz is None
    assert cfg.json_fg is None
    assert cfg.center is None
    assert cfg.rms is False
    assert cfg.rms_value == 0.0
    assert cfg.rms_sigma == 3.0
    assert cfg.niter == 5000
    assert cfg.output_prefix is None
    assert cfg.overwrite is False
    assert cfg.cleaning is False
    # nested configurations (post refactor)
    assert isinstance(cfg.observation, ObsConfig)
    assert isinstance(cfg.imaging, ImgConfig)
    assert cfg.observation.freq_mhz == 700.0
    assert cfg.imaging.pixels == 512


def test_sim_config_explicit_telescope_version():
    """telescope with version."""
    cfg = SimConfig(telescope="SKA1LOW", telescope_version="AA0.5")
    assert cfg.telescope == "SKA1LOW"
    assert cfg.telescope_version == "AA0.5"


def test_sim_config_accepts_named_catalogue():
    """built-in catalogues are selected by name."""
    cfg = SimConfig(catalogue="MIGHTEE")
    assert cfg.catalogue == "MIGHTEE"


def test_sim_config_scalar_I_rejected():
    """I must be a list; error otherwise"""
    # this error was documented in the original pipeline
    with pytest.raises(ValidationError):
        SimConfig(I=10.0)


def test_sim_config_numeric_catalogue_has_migration_message():
    """numeric catalogue IDs are removed in 0.2."""
    with pytest.raises(ValidationError, match="Numeric catalogue IDs were removed"):
        SimConfig(catalogue=1)


def test_sim_config_rejects_model_and_catalogue_together():
    """0.2 accepts one explicit sky model source per run."""
    with pytest.raises(ValidationError, match="one sky model source"):
        SimConfig(sky_file="sources.json", catalogue="GLEAM")


def test_sim_config_nested_observation_override():
    """nested ObsConfig can be fully replaced."""
    cfg = SimConfig(observation=ObsConfig(seconds=1200, freq_mhz=900.0))
    assert cfg.observation.seconds == 1200
    assert cfg.observation.freq_mhz == pytest.approx(900.0)
    # other nested defaults remain untouched
    assert cfg.imaging.pixels == 512


def test_sim_config_nested_imaging_override():
    """nested ImgConfig can be fully replaced."""
    cfg = SimConfig(imaging=ImgConfig(pixels=2048, algorithm="wsclean_clean"))
    assert cfg.imaging.pixels == 2048
    assert cfg.imaging.algorithm == "wsclean_clean"
    # other nested defaults remain untouched
    assert cfg.observation.seconds == 600
