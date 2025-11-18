# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
# Check if venv sphinx-build exists, otherwise use system sphinx-build
SPHINXBUILD   ?= $(shell test -f ../venv/bin/sphinx-build && echo ../venv/bin/sphinx-build || echo sphinx-build)
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile html_en gettext_update

# Build English version
html_en:
	@$(SPHINXBUILD) -b html -D language='en' "$(SOURCEDIR)" "$(BUILDDIR)/html_en" $(SPHINXOPTS) $(O)

# Update translation files
gettext_update:
	@$(SPHINXBUILD) -b gettext "$(SOURCEDIR)" "$(BUILDDIR)/gettext" $(SPHINXOPTS) $(O)
	../venv/bin/sphinx-intl update -p "$(BUILDDIR)/gettext" -l en

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
