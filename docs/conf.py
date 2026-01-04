# Configuration file for the Sphinx documentation builder
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
from datetime import datetime

# -- Path setup --------------------------------------------------------------

# Add the Trefftz package to sys.path so autodoc can import it
sys.path.insert(0, os.path.abspath("../trefftz"))

# -- Project information -----------------------------------------------------

project = "Trefftz"
author = "Manuel Pena"
copyright = f"{datetime.now().year}, {author}"
version = "0.1.0"
release = version

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",      # Auto-generate documentation from docstrings
    "sphinx.ext.napoleon",     # Support for Google / NumPy style docstrings
    "sphinx.ext.viewcode",     # Add links to source code
    "sphinx.ext.intersphinx",  # Link to other projects (e.g., NumPy, SciPy)
    "sphinx.ext.todo",         # Optional: include TODOs
]

# Type hints in documentation
autodoc_typehints = "description"

# Templates
templates_path = ["_templates"]

# Exclude build / cache folders
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- HTML output -------------------------------------------------------------

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]

# Custom CSS (optional)
# html_css_files = ["custom.css"]

# -- Autodoc options ---------------------------------------------------------

autodoc_member_order = "bysource"
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}

# -- Intersphinx mapping -----------------------------------------------------

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "matplotlib": ("https://matplotlib.org/stable/", None),
}

# -- Todo options ------------------------------------------------------------

todo_include_todos = True

