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

# html_theme = "furo"
#
# html_logo = "images/multipy-logo.svg"
#
# html_theme_options = {
#     "dark_css_variables": {
#         "color-problematic": "#b30000",
#         "color-foreground-primary": "black",
#         "color-foreground-secondary": "#5a5c63",
#         "color-foreground-muted": "#72747e",
#         "color-foreground-border": "#878787",
#         "color-background-primary": "white",
#         "color-background-secondary": "#f8f9fb",
#         "color-background-hover": "#efeff4ff",
#         "color-background-hover--transparent": "#efeff400",
#         "color-background-border": "#eeebee",
#         "color-inline-code-background": "#f2f2f2",
#
#         # Announcements
#         "color-announcement-background": "#000000dd",
#         "color-announcement-text": "#eeebee",
#
#         # Brand colors
#         "color-brand-primary": "#2962ff",
#         "color-brand-content": "#2a5adf",
#
#         # Highlighted text (search)
#         "color-highlighted-background": "#ddeeff",
#
#         # GUI Labels
#         "color-guilabel-background": "#ddeeff80",
#         "color-guilabel-border": "#bedaf580",
#
#         # API documentation
#         "color-api-highlight-on-target": "#ffffcc",
#
#         # Admonitions
#         "color-admonition-background": "transparent",
#     },
# }
