"""Release documentation checks."""

from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_environment_file_defines_skasim_conda_runtime():
    """The full simulation environment is named skasim and includes Karabo."""
    data = yaml.safe_load((REPO_ROOT / "environment.yml").read_text(encoding="utf-8"))

    assert data["name"] == "skasim"
    assert "conda-forge" in data["channels"]
    assert "i4ds" in data["channels"]
    assert "karabo-pipeline" in data["dependencies"]


def test_installation_docs_describe_pip_and_conda_paths():
    """Installation docs separate lightweight pip use from full simulations."""
    text = (REPO_ROOT / "docs/installation.rst").read_text(encoding="utf-8")

    assert "pip-only" in text
    assert "conda create -y -n skasim python=3.10" in text
    assert "karabo-pipeline \"cuda-version=11.7\"" in text
    assert "conda run -n skasim pip install -e ." in text
    assert "import karabo, oskar" in text
    assert "oskar_sim_interferometer" in text
    assert "Karabo extra" not in text
    assert "--wsclean-command" in text
    assert "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean" in text
    assert "smoke_mightee_wsclean_quick" in text
    assert "example" in text.lower()


def test_release_examples_use_0_2_cli_language():
    """Release-facing examples use 0.2 CLI flags."""
    examples = "\n".join(
        [
            (REPO_ROOT / "README.md").read_text(encoding="utf-8"),
            (REPO_ROOT / "docs/guide.rst").read_text(encoding="utf-8"),
            (REPO_ROOT / "docs/examples.rst").read_text(encoding="utf-8"),
        ]
    )

    assert "--flux-density" in examples
    assert "--stokes-q" in examples
    assert "--stokes-u" in examples
    assert "--stokes-v" in examples
    assert "--imager" in examples
    assert "--wsclean-command" in examples
    assert "--catalog MIGHTEE" in examples
    assert "--catalog 1" not in examples
    assert "--frequency-mhz" in examples
    assert "--observation-time" in examples


def test_release_docs_document_removed_0_1_surface():
    """Release docs state that deprecated 0.1 CLI and config fields are removed."""
    text = "\n".join(
        [
            (REPO_ROOT / "README.md").read_text(encoding="utf-8"),
            (REPO_ROOT / "docs/guide.rst").read_text(encoding="utf-8"),
        ]
    )

    assert "--I" in text
    assert "--cleaning" in text
    assert "--stokes-i" in text
    assert "--json" in text
    assert "--json-fg" in text
    assert "--freq" in text
    assert "--seconds" in text
    assert "--prefix" in text
    assert "--niter" in text
    assert "SimConfig.I" in text
    assert "SimConfig.Q" in text
    assert "SimConfig.U" in text
    assert "SimConfig.V" in text
    assert "SimConfig.cleaning" in text
    assert "SimConfig.output_prefix" in text
    assert "SimConfig.niter" in text
    assert "SimConfig.scale_I" in text
    assert "ImgConfig.algorithm" in text


def test_smoke_docs_include_named_catalog_and_fits_shapes():
    """Smoke docs include named-catalog and file-backed FITS run shapes."""
    text = (REPO_ROOT / "docs/examples.rst").read_text(encoding="utf-8")

    assert "MeerKAT" in text
    assert "--catalog MIGHTEE" in text
    assert "--imager wsclean" in text
    assert "--wsclean-command" in text
    assert "--output-dir" in text
    assert "--clean-iterations" in text
    assert "--model" in text
    assert "FITS" in text
