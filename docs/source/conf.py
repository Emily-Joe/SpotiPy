# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('../..'))

# -- Project information
project = 'SpotiPy'
copyright = '2026, Saqib Sumra'
author = 'Saqib Sumra'
release = 'v1.0.0'

# -- General configuration
extensions = [
    'sphinx.ext.napoleon',      # supports Google/NumPy style docstrings
    'sphinx.ext.viewcode',      # adds [source] links to API docs
    'sphinx.ext.autodoc',       # pulls docstrings from modules
    'autoapi.extension',        # auto-generates full API docs
    'myst_parser',              # write docs in Markdown
]

# AutoAPI settings - point to your SpotiPy source folder
autoapi_dirs = ['../../src']
autoapi_type = 'python'
autoapi_add_toctree_entry = True

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- HTML output
html_theme = 'furo'
html_static_path = ['_static']
