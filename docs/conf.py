"""Sphinx configuration for skasim."""

import os
import sys

# path to src/ so autodoc can import skasim
sys.path.insert(0, os.path.abspath("../src"))

project = "skasim"
copyright = "2026, Spanish SRC Team"
author = "Spanish SRC Team"
release = "0.1.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
]

# mock heavy dependencies not available doc build env
autodoc_mock_imports = [
    "karabo",
    "karabo_interferometer_simulation",
    "karabo.util.file_handler",
    "rascil",
    "ska_sdp_datamodels",
]

# napoleon settings for Google-style docstrings
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False

# autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]
