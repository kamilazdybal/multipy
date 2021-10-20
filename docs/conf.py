# -- Path setup --------------------------------------------------------------

import os
import sys
sys.path.insert(0, os.path.abspath('../multipy/'))

# -- Project information -----------------------------------------------------

project = 'multipy'
copyright = '2021, James C. Sutherland, Kamila Zdybał'
author = 'James C. Sutherland, Kamila Zdybał'
release = '1.0.0'

extensions = [
    "sphinx.ext.autodoc",
]

autosectionlabel_prefix_document = True

templates_path = ['_templates']

source_suffix = '.rst'

master_doc = 'index'

language = 'English'

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

pygments_style = 'sphinx'

html_static_path = []

import jupyter_sphinx_theme
html_theme = "jupyter"
html_sidebars = {'**': ['sidebartoc.html']}
html_theme_path = jupyter_sphinx_theme.get_html_theme_path()
