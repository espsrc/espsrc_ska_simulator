"""tests/test_pipeline.py — unit tests for pipeline functions that don't invoke Karabo simulators."""

import json
import os
from enum import Enum
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
    resolve_telescope_version,
    run_simulation,
    source_ref_get_best_observation_time,
)


class _FakeTelescope:
    def plot_telescope(self, file):
        Path(file).write_text("plot", encoding="utf-8")


class _FakeVersions(Enum):
    SKA_OST_ARRAY_CONFIG_2_3_1 = "ska-ost-array-config-2.3.1"

# --------------------------------------------------------------------------- #
# telescope
# --------------------------------------------------------------------------- #


def test_resolve_telescope_version_accepts_enum_name():
    """CLI telescope version names are converted to Karabo enum members."""
    module = type(
        "Module",
        (),
        {"OSKAR_TELESCOPE_TO_VERSIONS": {"SKA-LOW-AAstar": _FakeVersions}},
    )

    version = resolve_telescope_version(
        module,
        "SKA-LOW-AAstar",
        "SKA_OST_ARRAY_CONFIG_2_3_1",
    )

    assert version is _FakeVersions.SKA_OST_ARRAY_CONFIG_2_3_1


def test_resolve_telescope_version_accepts_enum_value():
    """Karabo enum values are also accepted for telescope versions."""
    module = type(
        "Module",
        (),
        {"OSKAR_TELESCOPE_TO_VERSIONS": {"SKA-LOW-AAstar": _FakeVersions}},
    )

    version = resolve_telescope_version(
        module,
        "SKA-LOW-AAstar",
        "ska-ost-array-config-2.3.1",
    )

    assert version is _FakeVersions.SKA_OST_ARRAY_CONFIG_2_3_1


def test_resolve_telescope_version_rejects_unknown_version():
    """Bad telescope versions fail before Karabo raises opaque enum errors."""
    module = type(
        "Module",
        (),
        {"OSKAR_TELESCOPE_TO_VERSIONS": {"SKA-LOW-AAstar": _FakeVersions}},
    )

    with pytest.raises(ValueError, match="Accepted versions"):
        resolve_telescope_version(module, "SKA-LOW-AAstar", "bad-version")


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
        observation=ObsConfig(frequency_mhz=700, observation_time_s=1),
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
        observation=ObsConfig(frequency_mhz=700, observation_time_s=1),
        imaging=ImgConfig(fov_deg=None),
    )
    freq = 700 * u.MHz
    fov = compute_fov(config, freq)
    wavelength = freq.to(u.m, equivalencies=u.spectral()).value
    expected_rad = 1.25 * wavelength / expected_diameter_m
    assert fov.to(u.rad).value == pytest.approx(expected_rad, rel=1e-6)


def test_compute_fov_positive():
    """Default config produces a positive FoV."""
    config = SimConfig(observation=ObsConfig(observation_time_s=1))
    freq = config.observation.frequency_mhz * u.MHz
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
        config = SimConfig(output_dir="test_run", telescope="SKA1MID")
        ctx = create_run_context(config)
        assert ctx.work_dir.is_absolute()
        assert ctx.work_dir.name == "test_run"
        assert ctx.work_dir.is_dir()
        assert ctx.log_path.exists()
        assert ctx.manifest_path.exists()
        assert ctx.manifest.run_id == "test_run"
        assert ctx.manifest.status == "running"
    finally:
        os.chdir(old_cwd)


def test_create_run_context_no_prefix_uses_timestamp(tmp_path):
    """When output_dir is None, directory name includes a timestamp-like string."""
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
    """JSON catalog loads successfully and returns a SkyModel with one source."""
    data = [_single_source_json()]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky = _load_sky_from_file("catalog.json")
    assert sky.sources is not None
    assert len(sky.sources) == 1


def test_load_sky_from_json_preserves_shape_metadata():
    """JSON sky models retain source ellipse metadata for previews."""
    src = _single_source_json()
    src["major_axis"] = 120.0
    src["minor_axis"] = 30.0
    src["pa"] = 42.0
    src["spec_index"] = -0.63
    mopen = mock_open(read_data=json.dumps([src]))

    with patch("builtins.open", mopen):
        sky = _load_sky_from_file("catalog.json")

    rendered = sky.to_json()[0]
    assert rendered["major_axis"] == pytest.approx(120.0)
    assert rendered["minor_axis"] == pytest.approx(30.0)
    assert rendered["pa"] == pytest.approx(42.0)
    assert rendered["spec_index"] == pytest.approx(-0.63)


def test_load_sky_from_json_scales_intensity():
    """flux_scale != 1 multiplies source intensities."""
    data = [_single_source_json()]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky_before = _load_sky_from_file("catalog.json", flux_scale=1.0)
    mopen2 = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen2):
        sky_after = _load_sky_from_file("catalog.json", flux_scale=3.0)
    assert len(sky_after.sources) == 1
    assert sky_after.sources[0, 2] == pytest.approx(sky_before.sources[0, 2] * 3.0)


def test_load_sky_from_json_assigns_ref_freq_when_zero():
    """If JSON source has ref_freq == 0, the loader assigns observing frequency."""
    src = _single_source_json()
    src["ref_freq"] = 0
    data = [src, {**src, "ra": 150.1}, {**src, "ra": 149.9}]
    mopen = mock_open(read_data=json.dumps(data))
    with patch("builtins.open", mopen):
        sky = _load_sky_from_file("catalog.json", frequency=1420 * u.MHz)
    assert sky.sources is not None
    assert len(sky.sources) == 3
    assert sky.to_json()[0]["ref_freq"] == pytest.approx(1420e6)


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


def test_build_sky_model_file_records_sky_model_output(tmp_path):
    """File-backed sky model sources are recorded as sky_model outputs."""
    source_file = tmp_path / "sources.json"
    source_file.write_text(json.dumps([_single_source_json()]), encoding="utf-8")
    config = SimConfig(sky_file=str(source_file))
    ctx = _make_ctx(tmp_path, config)

    build_sky_model(ctx, fov=0.2 * u.deg)

    sky_outputs = [
        output for output in ctx.manifest.outputs if output.kind == "sky_model"
    ]
    assert sky_outputs[0].path == str(source_file.resolve())


def test_build_sky_model_random_source_count(tmp_path):
    """Random source generation produces one source per configured intensity."""
    config = SimConfig(source_flux_jy=[1.0, 5.0, 10.0, 20.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)
    assert len(sky.sources) == 4
    assert center is not None


def test_build_sky_model_random_single_source(tmp_path):
    """A single intensity produces one source."""
    config = SimConfig(source_flux_jy=[42.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)
    assert len(sky.sources) == 1


def test_build_sky_model_uses_source_flux_jy(tmp_path):
    """Generated source intensities create one generated source per value."""
    config = SimConfig(source_flux_jy=[1.0, 5.0, 10.0])
    ctx = _make_ctx(tmp_path, config)
    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)

    assert len(sky.sources) == 3
    assert center is not None


def test_build_sky_model_uses_generated_source_polarization(tmp_path):
    """Generated source Stokes Q/U/V values are passed to the sky model."""
    config = SimConfig(
        source_flux_jy=[1.0, 5.0],
        stokes_q_jy=[0.1, 0.2],
        stokes_u_jy=[0.0, 0.3],
        stokes_v_jy=[0.0, -0.1],
    )
    ctx = _make_ctx(tmp_path, config)
    sky, _ = build_sky_model(ctx, fov=0.2 * u.deg)

    assert len(sky.sources) == 2
    assert sky.sources[0, 3] == pytest.approx(0.1)
    assert sky.sources[1, 3] == pytest.approx(0.2)
    assert sky.sources[1, 4] == pytest.approx(0.3)
    assert sky.sources[1, 5] == pytest.approx(-0.1)


@pytest.mark.parametrize(
    ("catalog", "loader_name"),
    [
        ("MIGHTEE", "get_MIGHTEE_Sky"),
        ("GLEAM", "get_GLEAM_Sky"),
    ],
)
def test_build_sky_model_named_catalog(tmp_path, monkeypatch, catalog, loader_name):
    """A named built-in catalog selects the matching catalog source."""
    fake_sky = SkyModel(np.array([[10.0, 20.0, 1.0]]))
    monkeypatch.setattr(
        SkyModel,
        loader_name,
        staticmethod(lambda: fake_sky),
        raising=False,
    )
    config = SimConfig(catalog=catalog)
    ctx = _make_ctx(tmp_path, config)

    sky, center = build_sky_model(ctx, fov=0.2 * u.deg)

    assert sky is fake_sky
    assert center.ra.value == pytest.approx(10.0)
    assert ctx.manifest.milestones[-1].details["format"] == catalog


# --------------------------------------------------------------------------- #
# build_sky_model — invalid configurations
# --------------------------------------------------------------------------- #


def test_build_sky_model_unsupported_catalog_raises(tmp_path):
    """Unsupported catalog names raise ValueError inside build_sky_model."""
    config = SimConfig(observation=ObsConfig(observation_time_s=1), catalog="MIGHTEE")
    # bypass pydantic validation; build_sky_model has its own ValueError guard
    object.__setattr__(config, "catalog", "UNKNOWN")
    ctx = _make_ctx(tmp_path, config)
    with pytest.raises(ValueError, match="Catalog UNKNOWN not available"):
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
            output_dir="imager_run",
            imaging=ImgConfig(imager="wsclean"),
            observation=ObsConfig(observation_time_s=1),
        )
        pipeline.run(config)

        manifest = json.loads(
            (tmp_path / "imager_run" / "run_manifest.json").read_text(
                encoding="utf-8"
            )
        )
    finally:
        os.chdir(old_cwd)

    assert called == ["wsclean"]
    assert manifest["config"]["imaging"]["imager"] == "wsclean"
    assert (tmp_path / "imager_run" / "weblog.html").exists()
    assert any(output["kind"] == "weblog" for output in manifest["outputs"])
    imaging_done = [
        item for item in manifest["milestones"] if item["name"] == "imaging_completed"
    ][0]
    assert imaging_done["details"]["Imager"] == "wsclean"


def test_run_renders_weblog_on_failure(tmp_path, monkeypatch):
    """Failed runs save the manifest and render a weblog with the error."""
    import skasim.pipeline as pipeline

    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        monkeypatch.setattr(
            pipeline,
            "build_telescope",
            lambda ctx: (_ for _ in ()).throw(RuntimeError("boom")),
        )
        config = SimConfig(output_dir="failed_run")

        with pytest.raises(RuntimeError, match="boom"):
            pipeline.run(config)

        work_dir = tmp_path / "failed_run"
        manifest = json.loads((work_dir / "run_manifest.json").read_text())
        weblog = (work_dir / "weblog.html").read_text(encoding="utf-8")
    finally:
        os.chdir(old_cwd)

    assert manifest["status"] == "failed"
    assert "boom" in manifest["errors"][0]
    assert any(output["kind"] == "weblog" for output in manifest["outputs"])
    assert "failed" in weblog
    assert "boom" in weblog


def test_run_records_failure_milestone_details_as_dict(tmp_path, monkeypatch):
    """Failure milestones remain serializable when a phase raises."""
    import skasim.pipeline as pipeline

    old_cwd = os.getcwd()
    os.chdir(tmp_path)
    try:
        center = SkyCoord(10 * u.deg, 20 * u.deg)
        fake_sky = MagicMock()
        fake_sky.get_center.return_value = center
        monkeypatch.setattr(
            pipeline,
            "build_telescope",
            lambda ctx: _FakeTelescope(),
        )
        monkeypatch.setattr(pipeline, "compute_fov", lambda config, freq: 0.2 * u.deg)
        monkeypatch.setattr(
            pipeline,
            "build_sky_model",
            lambda ctx, fov: (fake_sky, center),
        )
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
            lambda ctx, telescope, observation, sky_model: (_ for _ in ()).throw(
                RuntimeError("phase failed")
            ),
        )
        config = SimConfig(output_dir="failed_milestone")

        with pytest.raises(RuntimeError, match="phase failed"):
            pipeline.run(config)

        manifest = json.loads(
            (tmp_path / "failed_milestone" / "run_manifest.json").read_text(
                encoding="utf-8"
            )
        )
    finally:
        os.chdir(old_cwd)

    failed = [
        item for item in manifest["milestones"] if item["name"] == "simulation_failed"
    ][0]
    assert failed["details"] == {"error": "phase failed"}


def test_run_simulation_does_not_rebuild_observation(tmp_path, monkeypatch):
    """Simulation uses the already-built observation instead of duplicating setup."""
    class FakeInterferometerSimulation:
        params = None
        run_args = None

        def __init__(self, **params):
            FakeInterferometerSimulation.params = params

        def run_simulation(self, **kwargs):
            FakeInterferometerSimulation.run_args = kwargs

    def fake_require(module_name):
        if module_name == "karabo.simulation.interferometer":
            return MagicMock(InterferometerSimulation=FakeInterferometerSimulation)
        if module_name == "karabo.simulator_backend":
            return MagicMock(SimulatorBackend=MagicMock(OSKAR="OSKAR"))
        raise AssertionError(module_name)

    monkeypatch.setattr("skasim.pipeline.require_karabo_module", fake_require)
    monkeypatch.setattr(
        "skasim.pipeline.build_observation",
        lambda *args, **kwargs: pytest.fail("build_observation should not run"),
    )
    config = SimConfig(output_dir=str(tmp_path / "sim"))
    ctx = create_run_context(config)
    sky_model = MagicMock()

    visibility_path = run_simulation(ctx, object(), object(), sky_model)

    assert visibility_path == ctx.visibility_path
    assert FakeInterferometerSimulation.params["channel_bandwidth_hz"] == pytest.approx(
        config.observation.channel_width_mhz * 1e6
    )
    assert FakeInterferometerSimulation.run_args["observation"] is not None
