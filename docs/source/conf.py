# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
import os
import sys

sys.path.insert(0, os.path.abspath('../../src'))

project = 'FastGEMF'
copyright = '2024, Mohammad Hossein Samaei'
author = 'Mohammad Hossein Samaei'
release = 'v0.1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',        # Automatically generate docs from docstrings
    'sphinx.ext.napoleon',       # Supports Google and NumPy style docstrings
    'sphinx.ext.viewcode',       # Adds links to highlighted source code
    'sphinx.ext.todo',           # Supports todo comments in docs
    #'sphinx.ext.mathjax ' ,        
]

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
