# docs/conf.py
import os
import sys
from datetime import date

# -- Project info -----------------------------------------------------
project = "toxindex"
author = "Your Name"
copyright = f"{date.today().year}, {author}"
release = "0.1.0"

# -- Path setup: make project importable ------------------------------
# Add the repo root to sys.path (so 'workflows', 'webserver', etc. can be imported)
sys.path.insert(0, os.path.abspath(".."))

# -- General config ---------------------------------------------------
extensions = [
    "sphinx.ext.autodoc",   # pull in docstrings
    "sphinx.ext.viewcode",  # link to highlighted source
    "sphinx.ext.napoleon",  # Google/NumPy-style docstrings (built-in)
    "myst_parser",  # Markdown support
]

# Autodoc behavior tweaks
autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "show-inheritance": True,
}
autoclass_content = "both"  # include class docstring + __init__ docstring
autodoc_typehints = "description"

# -- HTML output ------------------------------------------------------
# html_theme = "furo"  # falls back to default if not installed