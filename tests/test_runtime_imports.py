"""Runtime import behavior for lightweight skasim use."""

import os
import subprocess
import sys
from pathlib import Path

import pytest

from skasim.config import ObsConfig, SimConfig


def test_cli_help_does_not_require_karabo():
    """CLI help is available before the full Karabo runtime is installed."""
    repo_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(repo_root / "src"), env.get("PYTHONPATH", "")]
    )

    result = subprocess.run(
        [sys.executable, "-m", "skasim.cli", "--help"],
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert "--model" in result.stdout
    assert "--catalogue" in result.stdout


def test_simulation_without_karabo_fails_with_installation_message(tmp_path):
    """Full simulation execution reports the missing Karabo runtime clearly."""
    from skasim.pipeline import run

    config = SimConfig(
        output_prefix=str(tmp_path / "missing_karabo"),
        observation=ObsConfig(seconds=1),
    )

    with pytest.raises(RuntimeError, match="Karabo.*conda.*skasim"):
        run(config)
