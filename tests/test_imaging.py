"""Image production behavior."""

from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pytest
from astropy.io import fits

from skasim.config import ImgConfig, ObsConfig, SimConfig
from skasim.imaging import (
    build_wsclean_argv,
    collect_wsclean_outputs,
    run_wsclean_command,
    wsclean_output_prefix,
    write_fits_preview,
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
