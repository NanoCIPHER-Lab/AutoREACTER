# docs/source/conf.py

import os
import sys

# Ensure Sphinx can find your AutoREACTER package folder
sys.path.insert(0, os.path.abspath('../..'))
from AutoREACTER import __version__
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AutoREACTER'
version = __version__
release = __version__
copyright = '2026, Janitha Mahanthe, Jacob Gissinger'
author = 'Janitha Mahanthe, Jacob Gissinger'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',      # Automatically pull docstrings from Python
    'sphinx.ext.napoleon',     # Support for Google/NumPy-style docstrings
    'sphinx.ext.viewcode',     # Add links to highlighted source code
    'myst_parser',             # Allow Markdown (.md) files!
    
]

templates_path = ['_templates']
exclude_patterns = []

source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}
# -- Options for HTML output -------------------------------------------------

# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-out
# ut

# html_theme = 'furo'
# html_static_path = ['_static']
# html_logo = "_static/logo.png"make html


html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_logo = "_static/logo.png"
