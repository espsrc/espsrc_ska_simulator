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


def test_release_examples_use_0_2_cli_language():
    """Release-facing examples use 0.2 CLI flags."""
    combined = "\n".join(
        [
            Path("README.md").read_text(encoding="utf-8"),
            Path("docs/guide.rst").read_text(encoding="utf-8"),
            Path("docs/examples.rst").read_text(encoding="utf-8"),
        ]
    )

    assert "--flux-density" in combined
    assert "--stokes-i" in combined
    assert "--imager" in combined
    assert "--wsclean-command" in combined
    assert "--catalogue MIGHTEE" in combined
    assert "--catalogue 1" not in combined
    assert "--cleaning" not in combined
    assert "--I" not in combined


def test_smoke_docs_include_named_catalogue_and_fits_shapes():
    """Smoke docs include named-catalogue and file-backed FITS run shapes."""
    text = Path("docs/examples.rst").read_text(encoding="utf-8")

    assert "MeerKAT" in text
    assert "--catalogue MIGHTEE" in text
    assert "--imager wsclean" in text
    assert "--wsclean-command" in text
    assert "--model" in text
    assert "FITS" in text
