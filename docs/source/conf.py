# docs/source/conf.py

import os
import sys

# Ensure Sphinx can find your AutoREACTER package folder
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[2]
INIT_FILE = ROOT / "AutoREACTER" / "__init__.py"

init_text = INIT_FILE.read_text(encoding="utf-8")
match = re.search(r'^__version__\s*=\s*[\'"]([^\'"]+)[\'"]', init_text, re.M)

if not match:
    raise RuntimeError("Could not find __version__ in AutoREACTER/__init__.py")

release = match.group(1)
version = release

try:
    myst_enable_extensions
except NameError:
    myst_enable_extensions = []

if "substitution" not in myst_enable_extensions:
    myst_enable_extensions.append("substitution")

myst_substitutions = {
    "autoreacter_version": release,
}
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'AutoREACTER'
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
