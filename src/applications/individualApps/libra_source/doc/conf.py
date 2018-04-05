# -*- coding: utf-8 -*-
#
# Default documentation build configuration file for Sphinx.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.

import sys
import os
import sphinx_rtd_theme

# General configuration
# ---------------------

sys.path.insert(0, '/sbiasfw/lab/basis/2.1.2/centos6/lib/python/basis/sphinx/ext/')
extensions = ['sphinx.ext.todo']  # Use 'rst2pdf.pdfbuilder' to enable pdf build, if rst2pdf is installed.

authors = ['Brad M. Keller', 'Michael M.K. Hsieh']
output_name = 'LIBRA_Software_Manual'.replace(' ', '_') # otherwise generated Makefile's do not work

today     = ''
today_fmt = '%B %d, %Y'

exclude_patterns = ['**/.svn', '**/.git']
exclude_patterns.extend([])

templates_path = []

source_suffix = '.rst'
master_doc = 'index'

project   = 'LIBRA'
version   = '1.0'
release   = '1.0.4'
copyright = '2014-2016 University of Pennsylvania'

highlight_language = 'matlab'

# Options for HTML output
# -----------------------

html_static_path = ['static']

#html_theme = "sphinx_rtd_theme"
#html_theme_path = [sphinx_rtd_theme.get_html_theme_path()]

html_theme = 'cbica'
html_theme_path = ['', './sphinx-themes']
if html_theme_path[0] == '': html_theme_path.pop(0)

html_theme_options = {'rellinks': ['about', 'download', 'manual', 'faq', 'publications', 'people']} # 'faq'

html_sidebars = {'**': [ 'localtoc.html', 'sourcelink.html', 'searchbox.html']}

html_style = ''
if html_style == '': html_style = None

html_logo = ''
if html_logo == '': html_logo = None

html_title = ''
if html_title == '': html_title = 'LIBRA'

html_short_title = ''
if html_short_title == '': html_short_title = html_title

htmlhelp_basename = 'libra_doc'


# Options for LaTeX output
# ------------------------

# This value determines how to group the document tree into LaTeX source files. 
# It must be a list of tuples (startdocname, targetname, title, author, documentclass, toctree_only)
latex_documents     = [(master_doc, output_name + '.tex', 'LIBRA Software Manual', ' \\and '.join(authors), 'howto', False)]
latex_use_parts     = False
latex_use_modindex  = False
latex_logo          = ''
if latex_logo == '': latex_logo = None
latex_show_urls     = 'no'
latex_show_pagerefs = False

PREAMBLE = """
\setcounter{page}{1}
\pagenumbering{arabic}
"""

latex_elements = {'printindex': '', 'fncychap': '', 'figure_align': 'H', 'preamble': PREAMBLE}


# Options for MAN output
# ----------------------

man_pages = [(master_doc, output_name, 'LIBRA package for documentation.', authors, 7)]


# Options for Texinfo output
# --------------------------

texinfo_documents = [(master_doc, output_name, '', '@*'.join(authors), 'LIBRA', 'LIBRA package for documentation.', 'SBIA', False)]


## Options for PDF output 
## ----------------------

## Grouping the document tree into PDF files. List of tuples
## (source start file, target name, title, author, options).
##
## If there is more than one author, separate them with \\.
## For example: r'Guido van Rossum\\Fred L. Drake, Jr., editor'
##
## The options element is a dictionary that lets you override
## this config per-document.
## For example,
## ('index', u'MyProject', u'My Project', u'Author Name',
##  dict(pdf_compressed = True))
## would mean that specific document would be compressed
## regardless of the global pdf_compressed setting.


#pdf_documents = [(master_doc, output_name + '.pdf', 'LIBRA Software Manual', ' \\and '.join(authors), 'howto', False)]
##[('index', u'EMSDocumentation', u'EMS Documentation', u'Thyag Sundaramoorthy'),]


## A comma-separated list of custom stylesheets. Example:
#pdf_stylesheets = ['sphinx','kerning','a4']


## A list of folders to search for stylesheets. Example:
#pdf_style_path = ['.', '_styles']


## Create a compressed PDF
## Use True/False or 1/0
## Example: compressed=True
##pdf_compressed = False


## A colon-separated list of folders to search for fonts. Example:
## pdf_font_path = ['/usr/share/fonts', '/usr/share/texmf-dist/fonts/']


## Language to be used for hyphenation support
##pdf_language = "en_US"


## Mode for literal blocks wider than the frame. Can be
## overflow, shrink or truncate
##pdf_fit_mode = "shrink"


## Section level that forces a break page.
## For example: 1 means top-level sections start in a new page
## 0 means disabled
##pdf_break_level = 0


## When a section starts in a new page, force it to be 'even', 'odd',
## or just use 'any'
##pdf_breakside = 'any'


## Insert footnotes where they are defined instead of
## at the end.
##pdf_inline_footnotes = True


## verbosity level. 0 1 or 2
##pdf_verbosity = 0


## If false, no index is generated.
##pdf_use_index = True


## If false, no modindex is generated.
##pdf_use_modindex = True


## If false, no coverpage is generated.
##pdf_use_coverpage = True


## Name of the cover page template to use
##pdf_cover_template = 'sphinxcover.tmpl'


## Documents to append as an appendix to all manuals.
##pdf_appendices = []


## Enable experimental feature to split table cells. Use it
## if you get "DelayedTable too big" errors
##pdf_splittables = False


## Set the default DPI for images
##pdf_default_dpi = 72


## Enable rst2pdf extension modules (default is only vectorpdf)
## you need vectorpdf if you want to use sphinx's graphviz support
##pdf_extensions = ['vectorpdf']


## Page template name for "regular" pages
##pdf_page_template = 'cutePage'


## Show Table Of Contents at the beginning?
##pdf_use_toc = True


## How many levels deep should the table of contents be?
#pdf_toc_depth = 9999


## Add section number to section references
#pdf_use_numbered_links = False


## Background images fitting mode
#pdf_fit_background_mode = 'scale'


# Options for extensions
# ----------------------

doxylink                = {}
breathe_projects        = {}
breathe_default_project = ''
