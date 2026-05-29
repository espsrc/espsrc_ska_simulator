"""Image production behavior."""

import importlib.util
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from skasim.config import ImgConfig, ObsConfig, SimConfig
from skasim.imaging import (
    SKY_MODEL_CMAP,
    _compact_source_mask,
    _flux_marker_sizes,
    _padded_limits,
    _plot_sky_model_sources,
    _sky_model_ellipses,
    _sky_model_position_angle,
    build_shadems_uv_coverage_argv,
    build_wsclean_argv,
    collect_wsclean_outputs,
    shadems_uv_coverage_env,
    run_wsclean_command,
    write_uv_coverage_plot,
    wsclean_output_prefix,
    write_fits_preview,
    write_sky_model_previews,
)
from skasim.manifest import create_run_context


def test_build_wsclean_argv_uses_default_command():
    """The default WSClean command builds an argv list starting with wsclean."""
    config = SimConfig(
        imaging=ImgConfig(imager="wsclean", wsclean_command="wsclean", pixels=256)
    )

    argv = build_wsclean_argv(
        config=config,
        visibility_path=Path("visibilities.MS"),
        fov=0.2 * u.deg,
        output_prefix="run-clean",
    )

    assert argv[0] == "wsclean"
    assert "-name" in argv
    assert argv[argv.index("-name") + 1] == "run-clean"
    assert argv[-1] == "visibilities.MS"


def test_build_wsclean_argv_parses_singularity_command():
    """Containerized WSClean commands are parsed into argv without shell text."""
    command = "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"
    config = SimConfig(
        imaging=ImgConfig(imager="wsclean", wsclean_command=command, pixels=256)
    )

    argv = build_wsclean_argv(
        config=config,
        visibility_path=Path("visibilities.MS"),
        fov=0.2 * u.deg,
        output_prefix="run-clean",
    )

    assert argv[:4] == [
        "singularity",
        "exec",
        "/mnt/software/containers/wsclean-3.10-dysco.sif",
        "wsclean",
    ]
    assert "-name" in argv


def test_build_wsclean_argv_caps_channels_out_to_available_channels():
    """Small smoke runs should not request more WSClean outputs than channels."""
    config = SimConfig(
        observation=ObsConfig(bandwidth_mhz=25.0, n_channels=2),
        imaging=ImgConfig(imager="wsclean", pixels=256),
    )

    argv = build_wsclean_argv(
        config=config,
        visibility_path=Path("visibilities.MS"),
        fov=0.2 * u.deg,
        output_prefix="run-clean",
    )

    assert argv[argv.index("-channels-out") + 1] == "2"


def test_run_wsclean_command_uses_argv_and_working_directory(tmp_path, monkeypatch):
    """WSClean execution avoids shell execution and uses an explicit cwd."""
    calls = []

    def fake_run(argv, **kwargs):
        calls.append((argv, kwargs))

        class Result:
            stdout = "ok"

        return Result()

    monkeypatch.setattr("skasim.imaging.subprocess.run", fake_run)

    result = run_wsclean_command(["wsclean", "-name", "run-clean"], tmp_path)

    assert result.stdout == "ok"
    assert calls[0][0] == ["wsclean", "-name", "run-clean"]
    assert calls[0][1]["cwd"] == str(tmp_path)
    assert calls[0][1]["shell"] is False
    assert calls[0][1]["check"] is True


def test_build_shadems_uv_coverage_argv_uses_verified_configuration(tmp_path):
    """shadeMS UV coverage uses U/V axes and a square canvas."""
    config = SimConfig(
        imaging=ImgConfig(
            shadems_command="python -m shade_ms",
            uv_coverage_canvas_size=600,
        )
    )

    argv = build_shadems_uv_coverage_argv(
        config=config,
        visibility_path=tmp_path / "visibilities.MS",
        output_dir=tmp_path,
        png_name="run_uvcoverage.png",
        title="run uv coverage",
    )

    assert argv[:3] == ["python", "-m", "shade_ms"]
    assert "--xaxis" in argv
    assert argv[argv.index("--xaxis") + 1] == "u"
    assert "--yaxis" in argv
    assert argv[argv.index("--yaxis") + 1] == "v"
    assert argv[argv.index("--xcanvas") + 1] == "600"
    assert argv[argv.index("--ycanvas") + 1] == "600"
    assert argv[argv.index("--spread-pix") + 1] == "2"
    assert "--no-lim-save" in argv


def test_shadems_uv_coverage_env_uses_writable_cache_dirs(tmp_path):
    """shadeMS receives writable Matplotlib and Numba cache directories."""
    env = shadems_uv_coverage_env(tmp_path)

    assert Path(env["MPLCONFIGDIR"]).is_dir()
    assert Path(env["NUMBA_CACHE_DIR"]).is_dir()
    assert Path(env["MPLCONFIGDIR"]).parent.parent == tmp_path
    assert Path(env["NUMBA_CACHE_DIR"]).parent.parent == tmp_path


def test_write_uv_coverage_plot_records_manifest_outputs(tmp_path, monkeypatch):
    """shadeMS plot generation records both PNG and command log outputs."""
    config = SimConfig(output_dir=str(tmp_path / "run"))
    ctx = create_run_context(config)
    visibility_path = ctx.work_dir / "visibilities.MS"
    visibility_path.mkdir()

    def fake_run(argv, work_dir):
        (work_dir / "run_uvcoverage.png").write_bytes(
            b"\x89PNG\r\n\x1a\n\x00\x00\x00\rIHDR"
            b"\x00\x00\x00\x01\x00\x00\x00\x01\x08\x02"
            b"\x00\x00\x00\x90wS\xde\x00\x00\x00\x00IEND\xaeB`\x82"
        )

        return SimpleNamespace(stdout="shadeMS ok\n", stderr="")

    monkeypatch.setattr("skasim.imaging.run_shadems_command", fake_run)

    png_path = write_uv_coverage_plot(ctx, visibility_path)

    assert png_path == ctx.work_dir / "run_uvcoverage.png"
    uv_outputs = [output for output in ctx.manifest.outputs if output.role == "uv_coverage"]
    assert [output.kind for output in uv_outputs] == ["plot", "log"]
    assert uv_outputs[0].path == "run_uvcoverage.png"
    assert uv_outputs[0].metadata["tool"] == "shadems"
    assert (ctx.work_dir / "run_uvcoverage_shadems.log").read_text(encoding="utf-8") == "shadeMS ok\n"


def test_collect_wsclean_outputs_matches_only_configured_prefix(tmp_path):
    """WSClean output discovery ignores files from other runs."""
    expected = [
        tmp_path / "run-clean-MFS-dirty.fits",
        tmp_path / "run-clean-MFS-image.fits",
    ]
    unrelated = [
        tmp_path / "old-run-MFS-image.fits",
        tmp_path / "wsclean-0000-image.fits",
    ]
    for path in expected + unrelated:
        path.write_text("fits", encoding="utf-8")

    outputs = collect_wsclean_outputs(tmp_path, "run-clean")

    assert outputs == expected


def test_wsclean_output_prefix_is_run_scoped(tmp_path):
    """WSClean -name uses a stable prefix scoped to the current run."""
    config = SimConfig(output_dir=str(tmp_path / "example"))
    ctx = create_run_context(config)

    prefix = wsclean_output_prefix(ctx)

    assert prefix == "example_wsclean"


def test_write_fits_preview_uses_aplpy_cmasher_renderer(tmp_path):
    """WSClean FITS previews are rendered as WCS-aware APLpy PNGs."""
    pytest.importorskip("aplpy")
    pytest.importorskip("cmasher")

    fits_path = tmp_path / "image.fits"
    png_path = tmp_path / "image.png"
    data = np.array(
        [
            [0.0, 1.0e-6, 2.0e-6],
            [3.0e-6, 4.0e-6, 5.0e-6],
            [6.0e-6, 7.0e-6, 8.0e-6],
        ]
    )
    header = fits.Header()
    header["CTYPE1"] = "RA---SIN"
    header["CTYPE2"] = "DEC--SIN"
    header["CRVAL1"] = 150.0
    header["CRVAL2"] = 2.0
    header["CRPIX1"] = 2.0
    header["CRPIX2"] = 2.0
    header["CDELT1"] = -1.0 / 3600.0
    header["CDELT2"] = 1.0 / 3600.0
    header["BUNIT"] = "Jy/beam"
    header["BMAJ"] = 2.0 / 3600.0
    header["BMIN"] = 1.5 / 3600.0
    header["BPA"] = 35.0
    fits.writeto(fits_path, data, header, overwrite=True)

    write_fits_preview(fits_path, png_path, "Ignored title")

    assert png_path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_write_sky_model_previews_writes_full_and_fov_pngs(tmp_path):
    """Sky model previews include full-catalog and FoV-matched views."""
    from astropy.coordinates import SkyCoord

    from skasim.sky import SkyModel

    sky_model = SkyModel(
        np.array(
            [
                [150.0, 2.0, 1.0, 0.0, 0.0, 0.0, 700e6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [150.1, 2.1, 0.5, 0.0, 0.0, 0.0, 700e6, 0.0, 0.0, 4.0, 2.0, 0.0, 0.0, 0.0],
            ]
        )
    )

    outputs = write_sky_model_previews(
        sky_model,
        SkyCoord(150.0 * u.deg, 2.0 * u.deg),
        1.0 * u.deg,
        tmp_path,
        "example",
    )

    assert outputs == [
        ("example_sky_model.png", "sky_model"),
        ("example_sky_model_fov.png", "sky_model_fov"),
    ]
    for path, _ in outputs:
        assert (tmp_path / path).read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


def test_sky_model_previews_use_equal_axis_scale(tmp_path, monkeypatch):
    """Sky-model FoV circles stay circular by using equal RA/Dec data scaling."""
    from matplotlib.axes import Axes

    aspect_calls = []
    box_calls = []
    original_set_aspect = Axes.set_aspect
    original_set_box_aspect = Axes.set_box_aspect

    def record_set_aspect(self, aspect, *args, **kwargs):
        aspect_calls.append((aspect, kwargs.get("adjustable")))
        return original_set_aspect(self, aspect, *args, **kwargs)

    def record_set_box_aspect(self, aspect, *args, **kwargs):
        box_calls.append(aspect)
        return original_set_box_aspect(self, aspect, *args, **kwargs)

    monkeypatch.setattr(Axes, "set_aspect", record_set_aspect)
    monkeypatch.setattr(Axes, "set_box_aspect", record_set_box_aspect)

    _plot_sky_model_sources(
        tmp_path / "sky.png",
        np.array([150.0]),
        np.array([2.0]),
        np.array([1.0]),
        np.array([0.0]),
        np.array([0.0]),
        np.array([0.0]),
        norm=None,
        title="sky",
        xlim=(150.5, 149.5),
        ylim=(1.5, 2.5),
        fov_circle=(150.0, 2.0, 0.5),
    )

    assert ("equal", "datalim") in aspect_calls
    assert 1 in box_calls


def test_sky_model_ellipses_preserve_source_shape_metadata():
    """Sky model source previews use major/minor axes and position angle."""
    ellipses = _sky_model_ellipses(
        np.array([150.0]),
        np.array([2.0]),
        np.array([18.0]),
        np.array([6.0]),
        np.array([35.0]),
    )

    assert len(ellipses) == 1
    ellipse = ellipses[0]
    assert ellipse.width == pytest.approx(18.0 / 3600.0)
    assert ellipse.height == pytest.approx(6.0 / 3600.0)
    assert ellipse.angle == pytest.approx(55.0)


def test_sky_model_position_angle_uses_astronomical_convention():
    """PA 0 is north-south and PA 90 is east-west in Matplotlib coordinates."""
    assert _sky_model_position_angle(0.0) == pytest.approx(90.0)
    assert _sky_model_position_angle(90.0) == pytest.approx(0.0)


def test_compact_source_mask_depends_on_plot_width():
    """Small sources become crosses when they are unresolved in the plotted FoV."""
    mask = _compact_source_mask(
        np.array([3.0, 100.0]),
        plot_width_deg=1.0,
    )

    assert mask.tolist() == [True, False]


def test_flux_marker_sizes_scale_with_flux_density():
    """Compact-source cross size encodes flux density."""
    sizes = _flux_marker_sizes(np.array([0.01, 1.0]))

    assert sizes[1] > sizes[0]


def test_reference_catalog_generator_writes_ds9_regions(tmp_path):
    """Reference JSON catalog generation has a matching DS9 region writer."""
    script = Path(__file__).resolve().parents[1] / "scripts" / "generate_gaussian_catalog.py"
    spec = importlib.util.spec_from_file_location("generate_gaussian_catalog", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    region_path = tmp_path / "catalog.reg"
    module.write_ds9_regions(
        [
            {
                "ra": 150.0,
                "dec": 2.0,
                "I": 1.2,
                "spec_index": -0.5,
                "major_axis": 10.0,
                "minor_axis": 4.0,
                "pa": 35.0,
            }
        ],
        region_path,
    )

    text = region_path.read_text(encoding="utf-8")
    assert "fk5" in text
    assert 'ellipse(150.0000000000,2.0000000000,10.000000",4.000000",55.000000)' in text
    assert "point(150.0000000000,2.0000000000) # point=cross 8" in text


def test_reference_catalog_generator_uses_broad_demo_distributions(tmp_path, monkeypatch):
    """Reference catalog fluxes and sizes exercise compact, faint, and extended sources."""
    script = Path(__file__).resolve().parents[1] / "scripts" / "generate_gaussian_catalog.py"
    spec = importlib.util.spec_from_file_location("generate_gaussian_catalog", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.chdir(tmp_path)

    module.generate_catalog(seed=42, n_sources=100)

    import json

    sources = json.loads((tmp_path / "demo_output" / "reference_gaussian_catalog.json").read_text())
    fluxes = np.array([source["I"] for source in sources])
    major_axes = np.array([source["major_axis"] for source in sources])
    axis_ratios = np.array([source["minor_axis"] / source["major_axis"] for source in sources])
    position_angles = np.array([source["pa"] for source in sources])
    spectral_indices = np.array([source["spec_index"] for source in sources])

    assert len(sources) == 100
    assert np.count_nonzero(fluxes >= 0.06) == 2
    assert np.count_nonzero((fluxes >= 1.0e-5) & (fluxes <= 2.0e-5)) >= 10
    assert np.median(fluxes) == pytest.approx(1.0e-3, abs=5.0e-4)
    assert np.count_nonzero(major_axes < 1.0) >= 5
    assert np.count_nonzero((major_axes >= 1.0) & (major_axes < 15.0)) >= 40
    assert np.count_nonzero((major_axes >= 15.0) & (major_axes < 60.0)) >= 15
    assert np.count_nonzero(major_axes >= 60.0) >= 5
    assert axis_ratios.min() < 0.35
    assert axis_ratios.max() > 0.9
    assert position_angles.min() < 10.0
    assert position_angles.max() > 170.0
    assert spectral_indices.mean() == pytest.approx(-0.5, abs=0.1)
    assert (tmp_path / "demo_output" / "reference_gaussian_catalog.reg").exists()


def test_sky_model_previews_use_reversed_colormap():
    """Sky model previews render brighter sources darker than faint sources."""
    from matplotlib.collections import PatchCollection

    ellipses = _sky_model_ellipses(
        np.array([150.0]),
        np.array([2.0]),
        np.array([18.0]),
        np.array([6.0]),
        np.array([35.0]),
    )
    collection = PatchCollection(ellipses, cmap=SKY_MODEL_CMAP)

    assert collection.cmap.name == "viridis_r"


def test_padded_limits_use_source_coordinates():
    """Full sky-model previews set axes from source coordinates, not 0..1 defaults."""
    lower, upper = _padded_limits(np.array([149.5, 150.5]))

    assert lower < 149.5
    assert upper > 150.5


def test_run_wsclean_imaging_uses_run_prefix_and_stable_outputs(
    tmp_path, monkeypatch
):
    """WSClean imaging records only stable outputs for the current prefix."""
    from skasim.imaging import run_wsclean_imaging

    class FakeImage:
        def __init__(self, path):
            self.path = path

        def plot(self, filename, **kwargs):
            Path(filename).write_text("png", encoding="utf-8")

    def fake_require(module_name):
        if module_name == "karabo.imaging.image":
            return SimpleNamespace(Image=FakeImage)
        if module_name == "karabo.imaging.imager_wsclean":
            return SimpleNamespace(TMP_PREFIX_CUSTOM="tmp", TMP_PURPOSE_CUSTOM="test")
        if module_name == "karabo.util.file_handler":
            return SimpleNamespace(
                FileHandler=lambda: SimpleNamespace(get_tmp_dir=lambda **kwargs: tmp_path)
            )
        raise AssertionError(module_name)

    captured_argv = []

    def fake_run(argv, work_dir):
        captured_argv.append(argv)
        prefix = argv[argv.index("-name") + 1]
        (work_dir / f"{prefix}-MFS-image.fits").write_text("fits", encoding="utf-8")
        (work_dir / "old-run-MFS-image.fits").write_text("old", encoding="utf-8")

        class Result:
            stdout = "ok"

        return Result()

    monkeypatch.setattr("skasim.imaging.require_karabo_module", fake_require)
    monkeypatch.setattr("skasim.imaging.run_wsclean_command", fake_run)
    monkeypatch.setattr(
        "skasim.imaging.write_fits_preview",
        lambda img_path, png_path, title: Path(png_path).write_text(
            "png", encoding="utf-8"
        ),
    )

    config = SimConfig(
        output_dir=str(tmp_path / "example"),
        imaging=ImgConfig(imager="wsclean"),
    )
    ctx = create_run_context(config)

    run_wsclean_imaging(ctx, ctx.visibility_path, 0.2 * u.deg)

    prefix = wsclean_output_prefix(ctx)
    output_paths = [output.path for output in ctx.manifest.outputs]
    assert captured_argv[0][captured_argv[0].index("-name") + 1] == prefix
    assert f"{prefix}-MFS-image.fits" in output_paths
    assert f"{prefix}-MFS-image.png" in output_paths
    assert all("old-run" not in output for output in output_paths)
    assert all("_bw" not in output for output in output_paths)


def test_run_wsclean_imaging_cleanup_keeps_run_scoped_outputs(tmp_path, monkeypatch):
    """Cleanup does not delete outputs when the run prefix starts like a temp file."""
    from skasim.imaging import run_wsclean_imaging

    class FakeImage:
        def __init__(self, path):
            self.path = path

        def plot(self, filename, **kwargs):
            Path(filename).write_text("png", encoding="utf-8")

    def fake_require(module_name):
        if module_name == "karabo.imaging.image":
            return SimpleNamespace(Image=FakeImage)
        if module_name == "karabo.imaging.imager_wsclean":
            return SimpleNamespace(TMP_PREFIX_CUSTOM="tmp", TMP_PURPOSE_CUSTOM="test")
        if module_name == "karabo.util.file_handler":
            return SimpleNamespace(
                FileHandler=lambda: SimpleNamespace(get_tmp_dir=lambda **kwargs: tmp_path)
            )
        raise AssertionError(module_name)

    def fake_run(argv, work_dir):
        prefix = argv[argv.index("-name") + 1]
        (work_dir / f"{prefix}-MFS-image.fits").write_text("fits", encoding="utf-8")
        (work_dir / "wsclean-0000-temp.fits").write_text("temp", encoding="utf-8")

        class Result:
            stdout = "ok"

        return Result()

    monkeypatch.setattr("skasim.imaging.require_karabo_module", fake_require)
    monkeypatch.setattr("skasim.imaging.run_wsclean_command", fake_run)
    monkeypatch.setattr(
        "skasim.imaging.write_fits_preview",
        lambda img_path, png_path, title: Path(png_path).write_text(
            "png", encoding="utf-8"
        ),
    )

    config = SimConfig(
        output_dir=str(tmp_path / "wsclean-00case"),
        imaging=ImgConfig(imager="wsclean"),
    )
    ctx = create_run_context(config)

    run_wsclean_imaging(ctx, ctx.visibility_path, 0.2 * u.deg)

    prefix = wsclean_output_prefix(ctx)
    output_paths = [output.path for output in ctx.manifest.outputs]
    assert f"{prefix}-MFS-image.fits" in output_paths
    assert not (ctx.work_dir / "wsclean-0000-temp.fits").exists()
