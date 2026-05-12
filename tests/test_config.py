"""tests/test_config.py"""

import pytest
from pydantic import ValidationError

from skasim.config import ImgConfig, ObsConfig, SimConfig

# ---------------------------------------------------------------------------
# ObsConfig
# ---------------------------------------------------------------------------


def test_obs_config_defaults():
    """default values."""
    obs = ObsConfig()
    assert obs.freq_mhz == 700.0
    assert obs.bandwidth_mhz == 100.0
    assert obs.n_channels == 8
    assert obs.delta_freq_mhz == 0.1
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
        delta_freq_mhz=12.5,
        seconds=3600,
        phase_center_ra_deg=150.0,
        phase_center_dec_deg=2.5,
    )
    assert obs.freq_mhz == pytest.approx(1420.0)
    assert obs.bandwidth_mhz == pytest.approx(200.0)
    assert obs.n_channels == 16
    assert obs.delta_freq_mhz == pytest.approx(12.5)
    assert obs.seconds == 3600
    assert obs.phase_center_ra_deg == pytest.approx(150.0)
    assert obs.phase_center_dec_deg == pytest.approx(2.5)


@pytest.mark.parametrize("bad_seconds", [0, -1, -100])
def test_obs_config_non_positive_seconds_raises(bad_seconds):
    """seconds must be strictly positive."""
    with pytest.raises(ValidationError, match="Observation time must be > 0"):
        ObsConfig(seconds=bad_seconds)


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
    assert cfg.catalogue == 0
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


def test_sim_config_scalar_I_rejected():
    """I must be a list; error otherwise"""
    # this error was documented in the original pipeline
    with pytest.raises(ValidationError):
        SimConfig(I=10.0)


def test_sim_config_invalid_catalogue():
    """non-existing catalogue, id out of {0,1,2,3}."""
    # TODO: review implementation of ID=3
    with pytest.raises(ValidationError):
        SimConfig(catalogue=4)


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
