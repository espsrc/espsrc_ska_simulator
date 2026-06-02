"""Release documentation verification."""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_environment_file_defines_skasim_conda_runtime():
    """environment.yml declares the skasim name and verified channels."""
    text = (REPO_ROOT / "environment.yml").read_text(encoding="utf-8")
    assert "name: skasim" in text
    assert "- nvidia/label/cuda-11.7.0" in text
    assert "- i4ds" in text
    assert "- conda-forge" in text


def test_installation_docs_describe_pip_and_conda_paths():
    """Installation docs separate lightweight pip use from full simulations."""
    text = (REPO_ROOT / "docs/installation.rst").read_text(encoding="utf-8")

    assert "pip-only" in text
    assert "conda create -y -n skasim python=3.10" in text
    assert 'karabo-pipeline "cuda-version=11.7"' in text
    assert "conda run -n skasim pip install -e ." in text
    assert "import karabo, oskar" in text
    assert "oskar_sim_interferometer" in text
    assert "Karabo extra" not in text
    assert "--wsclean-command" in text
    assert (
        "singularity exec /mnt/software/containers/wsclean-3.10-dysco.sif wsclean"
        in text
    )
    assert "shadeMS" in text
    assert "shadems" in text.lower()
    assert "smoke_mightee_wsclean_quick" in text
    assert "example" in text.lower()


def test_release_examples_use_0_2_cli_language():
    """Examples in docs use the unified skasim CLI instead of internal helpers."""
    text = (REPO_ROOT / "docs/installation.rst").read_text(encoding="utf-8")
    assert "skasim \\\n      --output-dir" in text
    assert "skasim --help" in text
    assert "skasim-pipeline" not in text
    assert "python -m skasim" not in text


def test_release_docs_document_removed_0_1_surface():
    """Documentation reflects that numeric catalog names are removed."""
    text = (REPO_ROOT / "docs/installation.rst").read_text(encoding="utf-8")
    assert "named mightee catalog" in text.lower()
    assert "--catalog MIGHTEE" in text


def test_smoke_docs_include_named_catalog_and_fits_shapes():
    """The smoke run verified for 0.2 is documented."""
    text = (REPO_ROOT / "docs/installation.rst").read_text(encoding="utf-8")
    assert "smoke_mightee_wsclean_quick" in text
    assert "visibilities.MS" in text
    assert "shadeMS UV-coverage" in text
    assert "WSClean FITS products" in text
    assert "PNG previews" in text
    assert "run_manifest.json" in text
    assert "weblog.html" in text
