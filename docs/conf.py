"""Sphinx configuration for the ThermoScreening documentation."""

import os
import sys

sys.path.insert(0, os.path.abspath(".."))

project = "ThermoScreening"
author = "Stefanie Kröll, Josef M. Gallmetzer"
copyright = "2026, the ThermoScreening authors"

try:
    from ThermoScreening.version import __version__ as release
except Exception:  # pragma: no cover - version file may be absent in a bare checkout
    release = ""
version = release

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False

# autodoc imports the package; the optional calculation backends (tblite, xtb,
# dftb+) are imported lazily, so only mock what might be absent on a minimal
# docs builder.
autodoc_mock_imports = ["tblite"]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "alabaster"
html_static_path = []
