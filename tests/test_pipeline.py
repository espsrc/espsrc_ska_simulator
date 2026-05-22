"""tests/test_pipeline.py — unit tests for pipeline functions that don't invoke Karabo simulators."""

import json
import os
from pathlib import Path
from unittest.mock import MagicMock, mock_open, patch

import astropy.units as u
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from skasim.config import ImgConfig, ObsConfig, SimConfig
from skasim.manifest import create_run_context
from skasim.sky import SkyModel
from skasim.pipeline import (
    _load_sky_from_file,
    _load_sky_from_fits,
    build_sky_model,
    compute_fov,
    parse_center,
    source_ref_get_best_observation_time,
)


class _FakeTelescope:
    def plot_telescope(self, file):
        Path(file).write_text("plot", encoding="utf-8")

# --------------------------------------------------------------------------- #
# parse_center
# --------------------------------------------------------------------------- #


def test_parse_center_none_returns_fallback():
    """None center string returns the fallback coordinate."""
    fallback = SkyCoord(10 * u.deg, 20 * u.deg)
    result = parse_center(None, fallback)
    assert result.ra.value == pytest.approx(10.0)
    assert result.dec.value == pytest.approx(20.0)


def test_parse_center_valid_hmsdms_string():
    """Space-separated HMS/DMS string parses correctly."""
    fallback = SkyCoord(0 * u.deg, 0 * u.deg)
    result = parse_center("10h01m35.1s 2d41m41s", fallback)
    assert result.ra.to(u.deg).value == pytest.approx(150.3962, abs=1e-3)
    assert result.dec.to(u.deg).value == pytest.approx(2.6947, abs=1e-4)


def test_parse_center_with_colons():
    """Colon-separated format also parses correctly."""
    fallback = SkyCoord(0 * u.deg, 0 * u.deg)
    result = parse_center("10:01:35.1 02:41:41", fallback)
    assert result.ra.to(u.deg).value == pytest.approx(150.3962, abs=1e-3)


def test_parse_center_invalid_returns_fallback():
    """Unparseable string falls back gracefully."""
    fallback = SkyCoord(5 * u.deg, -10 * u.deg)
    result = parse_center("not-a-coordinate", fallback)
    assert result.ra.value == pytest.approx(5.0)
    assert result.dec.value == pytest.approx(-10.0)


# --------------------------------------------------------------------------- #
# compute_fov
# --------------------------------------------------------------------------- #


def test_compute_fov_uses_explicit_fov_deg():
    """When fov_deg is set, return that value converted to radians."""
    config = SimConfig(
        observation=ObsConfig(freq_mhz=700, seconds=1),
        imaging=ImgConfig(fov_deg=2.5),
    )
    freq = 700 * u.MHz
    fov = compute_fov(config, freq)
    assert fov.to(u.deg).value == pytest.approx(2.5, abs=1e-6)


@pytest.mark.parametrize(
    "telescope,expected_diameter_m",
    [
        ("SKA1MID", 15.0),
        ("SKA1LOW", 38.0),
    ],
)
def test_compute_fov_diffraction_limit(telescope, expected_diameter_m):
    """When fov_deg is None, compute diffraction-limited FoV."""
    config = SimConfig(
        telescope=telescope,
        observation=ObsConfig(freq_mhz=700, seconds=1),
        imaging=ImgConfig(fov_deg=None),
    )
    freq = 700 * u.MHz
    fov = compute_fov(config, freq)
    wavelength = freq.to(u.m, equivalencies=u.spectral()).value
    expected_rad = 1.25 * wavelength / expected_diameter_m
    assert fov.to(u.rad).value == pytest.approx(expected_rad, rel=1e-6)


def test_compute_fov_positive():
    """Default config produces a positive FoV."""
    config = SimConfig(observation=ObsConfig(seconds=1))
    freq = config.observation.freq_mhz * u.MHz
    fov = compute_fov(config, freq)
    assert fov.to(u.deg).value > 0


# --------------------------------------------------------------------------- #
# create_run_context
# --------------------------------------------------------------------------- #


def test_create_run_context_creates_directory_and_log_file(tmp_path):
    """create_run_context creates working directory, returns RunContext with paths."""
    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        config = SimConfig(output_prefix="test_run", telescope="SKA1MID")
        ctx = create_run_context(config)
        assert ctx.work_dir.is_absolute()
        assert ctx.work_dir.name.startswith("test_run_SKA1MID")
        assert ctx.work_dir.is_dir()
        assert ctx.log_path.exists()
        assert ctx.manifest_path.exists()
        assert ctx.manifest.run_id.startswith("test_run_SKA1MID")
        assert ctx.manifest.status == "running"
    finally:
        os.chdir(old_cwd)


def test_create_run_context_no_prefix_uses_timestamp(tmp_path):
    """When output_prefix is None, directory name includes a timestamp-like string."""
    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        config = SimConfig(telescope="SKA1MID")
        ctx = create_run_context(config)
        assert len(ctx.work_dir.name) >= 8
        assert ctx.work_dir.is_dir()
        assert ctx.manifest_path.exists()
    finally:
        os.chdir(old_cwd)


# --------------------------------------------------------------------------- #
# helpers for JSON fixture
# --------------------------------------------------------------------------- #

def _make_ctx(tmp_path, config):
    """Build a RunContext inside tmp_path for tests that need one."""
    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        ctx = create_run_context(config)
        return ctx
    finally:
        os.chdir(old_cwd)

def _single_source_json():
    return {
        "ra": 150.0,
        "dec": 2.5,
        "I": 1.0,
        "Q": 0.0,
        "U": 0.0,
        "V": 0.0,
        "ref_freq": 700e6,
        "spec_index": 0.0,
        "rot_meas": 0.0,
        "major_axis": 0.0,
        "minor_axis": 0.0,
        "pa": 0.0,
        "true_redshift": 0.0,
        "obs_redshift": 0.0,
    }


# --------------------------------------------------------------------------- #
# _load_sky_from_file — inline JSON via mock_open
# --------------------------------------------------------------------------- #


def test_load_sky_from_json_file():
    """JSON catalogue loads successfully and returns a SkyModel with one source."""
    data = [_single_source_json()]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky = _load_sky_from_file("catalogue.json")
    assert sky.sources is not None
    assert len(sky.sources) == 1


def test_load_sky_from_json_scales_intensity():
    """scale_I != 1 multiplies source intensities."""
    data = [_single_source_json()]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky_before = _load_sky_from_file("catalogue.json", scale_I=1.0)
    mopen2 = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen2):
        sky_after = _load_sky_from_file("catalogue.json", scale_I=3.0)
    assert len(sky_after.sources) == 1
    assert sky_after.sources[0, 2] == pytest.approx(sky_before.sources[0, 2] * 3.0)


def test_load_sky_from_json_assigns_ref_freq_when_zero():
    """If JSON source has ref_freq == 0, the loader assigns ref_freq_hz or frequency."""
    src = _single_source_json()
    src["ref_freq"] = 0
    data = [src, {**src, "ra": 150.1}, {**src, "ra": 149.9}]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky = _load_sky_from_file("catalogue.json", ref_freq_hz=1.42e9)
    assert sky.sources is not None
    # reduced_form only exports (ra, dec, I) so ref_freq is not inspectable
    # through sources array; smoke test that the path completes without crash
    assert len(sky.sources) == 3


def test_load_sky_from_file_unknown_extension():
    """Passing an unsupported extension raises ValueError."""
    with pytest.raises(ValueError, match="Unsupported sky-file extension"):
        _load_sky_from_file("foo.txt")


def test_load_sky_from_file_empty_json_raises():
    """An empty JSON array raises ValueError about no sources."""
    mopen = mock_open(read_data="[]")
    with patch("builtins.open", mopen), pytest.raises(ValueError, match="No sources found in JSON"):
        _load_sky_from_file("empty.json")


# --------------------------------------------------------------------------- #
# build_sky_model — random source generation
# --------------------------------------------------------------------------- #


def test_build_sky_model_random_source_count(tmp_path):
    """Random source generation produces exactly len(I) sources."""
    config = SimConfig(I=[1.0, 5.0, 10.0, 20.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)
    assert len(sky.sources) == 4
    assert center is not None


def test_build_sky_model_random_single_source(tmp_path):
    """A single intensity produces one source."""
    config = SimConfig(I=[42.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)
    assert len(sky.sources) == 1


def test_build_sky_model_uses_source_intensities(tmp_path):
    """Generated source intensities create one generated source per value."""
    config = SimConfig(source_intensities=[1.0, 5.0, 10.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)

    assert len(sky.sources) == 3
    assert center is not None


@pytest.mark.parametrize(
    ("catalogue", "loader_name"),
    [
        ("MIGHTEE", "get_MIGHTEE_Sky"),
        ("GLEAM", "get_GLEAM_Sky"),
    ],
)
def test_build_sky_model_named_catalogue(tmp_path, monkeypatch, catalogue, loader_name):
    """A named built-in catalogue selects the matching catalogue source."""
    fake_sky = SkyModel(np.array([[10.0, 20.0, 1.0]]))
    monkeypatch.setattr(
        SkyModel,
        loader_name,
        staticmethod(lambda: fake_sky),
        raising=False,
    )
    config = SimConfig(catalogue=catalogue)
    ctx = _make_ctx(tmp_path, config)

    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)

    assert sky is fake_sky
    assert center.ra.value == pytest.approx(10.0)
    assert ctx.manifest.milestones[-1].details["format"] == catalogue


# --------------------------------------------------------------------------- #
# build_sky_model — invalid configurations
# --------------------------------------------------------------------------- #


def test_build_sky_model_unsupported_catalogue_raises(tmp_path):
    """Unsupported catalogue names raise ValueError inside build_sky_model."""
    config = SimConfig(observation=ObsConfig(seconds=1), catalogue="MIGHTEE")
    # bypass pydantic validation; build_sky_model has its own ValueError guard
    object.__setattr__(config, "catalogue", "UNKNOWN")
    ctx = _make_ctx(tmp_path, config)
    with pytest.raises(ValueError, match="Catalogue UNKNOWN not available"):
        build_sky_model(ctx, fov=0.5 * u.deg)


# --------------------------------------------------------------------------- #
# source_ref_get_best_observation_time
# --------------------------------------------------------------------------- #


def test_source_ref_get_best_observation_time():
    """Returns an astropy Time around culmination for a mock telescope location."""
    center = SkyCoord(150 * u.deg, 2.5 * u.deg)
    telescope = MagicMock()
    telescope.centre_latitude = -30.0
    telescope.centre_longitude = 116.0
    telescope.centre_altitude = 300.0
    best_time = source_ref_get_best_observation_time(center, telescope)
    assert best_time is not None
    assert hasattr(best_time, "iso")  # astropy Time


def test_run_uses_resolved_wsclean_imager(tmp_path, monkeypatch):
    """run() selects imaging from config.imaging.imager and records it."""
    import skasim.pipeline as pipeline

    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        center = SkyCoord(10 * u.deg, 20 * u.deg)
        fake_sky = MagicMock()
        fake_sky.get_center.return_value = center
        called = []

        monkeypatch.setattr(pipeline, "build_telescope", lambda ctx: _FakeTelescope())
        monkeypatch.setattr(pipeline, "compute_fov", lambda config, freq: 0.2 * u.deg)
        monkeypatch.setattr(pipeline, "build_sky_model", lambda ctx, fov: (fake_sky, center))
        monkeypatch.setattr(
            pipeline,
            "build_observation",
            lambda ctx, center, telescope: (
                object(),
                700 * u.MHz,
                100 * u.MHz,
                8,
                12.5 * u.MHz,
                650 * u.MHz,
            ),
        )
        monkeypatch.setattr(
            pipeline,
            "run_simulation",
            lambda ctx, telescope, observation, sky_model: ctx.visibility_path,
        )
        monkeypatch.setattr(
            pipeline,
            "run_dirty_imaging",
            lambda *args, **kwargs: pytest.fail("dirty imaging should not run"),
        )
        monkeypatch.setattr(
            pipeline,
            "run_wsclean_imaging",
            lambda *args, **kwargs: called.append("wsclean"),
        )

        config = SimConfig(
            output_prefix="imager_run",
            imaging=ImgConfig(imager="wsclean", algorithm="wsclean_clean"),
            observation=ObsConfig(seconds=1),
        )
        pipeline.run(config)

        manifest = json.loads(
            (tmp_path / "imager_run_SKA1MID" / "run_manifest.json").read_text(
                encoding="utf-8"
            )
        )
    finally:
        os.chdir(old_cwd)

    assert called == ["wsclean"]
    assert manifest["config"]["imaging"]["imager"] == "wsclean"
    imaging_done = [
        item for item in manifest["milestones"] if item["name"] == "imaging_completed"
    ][0]
    assert imaging_done["details"]["Imager"] == "wsclean"
