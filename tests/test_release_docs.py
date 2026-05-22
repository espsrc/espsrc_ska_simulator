"""Release documentation checks."""

from pathlib import Path

import yaml


def test_environment_file_defines_skasim_conda_runtime():
    """The full simulation environment is named skasim and includes Karabo."""
    data = yaml.safe_load(Path("environment.yml").read_text(encoding="utf-8"))

    assert data["name"] == "skasim"
    assert "conda-forge" in data["channels"]
    assert "i4ds" in data["channels"]
    assert "karabo-pipeline" in data["dependencies"]


def test_installation_docs_describe_pip_and_conda_paths():
    """Installation docs separate lightweight pip use from full simulations."""
    text = Path("docs/installation.rst").read_text(encoding="utf-8")

    assert "pip-only" in text
    assert "conda env create -f environment.yml" in text
    assert "Karabo extra" not in text
    assert "--wsclean-command" in text
    assert "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean" in text
    assert "example" in text.lower()
