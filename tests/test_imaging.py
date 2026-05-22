"""Image production behavior."""

from pathlib import Path

import astropy.units as u

from skasim.config import ImgConfig, SimConfig
from skasim.imaging import build_wsclean_argv, run_wsclean_command


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
