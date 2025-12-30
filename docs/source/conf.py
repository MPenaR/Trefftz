# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Trefftz'
copyright = '2025, Manuel Pena'
author = 'Manuel Pena'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = []

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']



extensions = [
    "sphinx.ext.autodoc",           # needed for automatic docstring extraction
    "sphinx.ext.napoleon",          # NumPy / Google style docstrings
    "sphinx.ext.viewcode",          # optional, links to source code
    "sphinx_autodoc_typehints",     # optional, shows type hints
]
