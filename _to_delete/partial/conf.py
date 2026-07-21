# Sphinx configuration for the RefET documentation
import refet

project = 'RefET'
author = 'Charles Morton'
release = refet.__version__
version = release

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'myst_parser',
]

# The class docstrings live on __init__ in this codebase
autoclass_content = 'init'
autodoc_member_order = 'bysource'
autosummary_generate = True

napoleon_google_docstring = True
napoleon_numpy_docstring = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
}

myst_enable_extensions = ['colon_fence']

exclude_patterns = ['_build']

html_theme = 'furo'
html_title = f'RefET {release}'
