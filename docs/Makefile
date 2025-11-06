# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SPHINXPROJ    = mineye-terranigma
SOURCEDIR     = source
BUILDDIR      = build
VENV          = ~/.venv/2025

# Check if we're in the correct virtual environment
CHECK_VENV := $(shell echo $$VIRTUAL_ENV)
EXPECTED_VENV := $(shell readlink -f $(VENV))

# Target to check and activate venv
check-venv:
ifndef SKIP_VENV_CHECK
	@if [ "$(CHECK_VENV)" != "$(EXPECTED_VENV)" ]; then \
		echo ""; \
		echo "⚠️  Wrong virtual environment!"; \
		echo ""; \
		echo "Current: $(CHECK_VENV)"; \
		echo "Expected: $(EXPECTED_VENV)"; \
		echo ""; \
		echo "Please activate the correct environment:"; \
		echo "  source $(VENV)/bin/activate"; \
		echo ""; \
		echo "Or skip this check with: SKIP_VENV_CHECK=1 make html"; \
		echo ""; \
		exit 1; \
	fi
endif

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile check-venv

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile check-venv
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

html-noplot: export SPHX_GLR_THEME = $(theme)
html-noplot: check-venv
	$(SPHINXBUILD) -D plot_gallery=0 -b html $(ALLSPHINXOPTS) $(SOURCEDIR) $(BUILDDIR)/html
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."

