
# Minimal makefile for Sphinx documentation

SHELL        := /usr/bin/env bash

# CONFIG
VENV         := $(HOME)/.venv/2025
PYTHON       := $(VENV)/bin/python
SPHINXBUILD  := $(PYTHON) -m sphinx       # ← module form guarantees venv python
SPHINXOPTS   :=
SPHINXPROJ   := mineye-terranigma
SOURCEDIR    := source
BUILDDIR     := build

# Sanity check: print which Python will be used (helpful in logs)
.PHONY: which-python
which-python:
	@echo "Using Python: $$($(PYTHON) -c 'import sys; print(sys.executable)')"

# Optional venv check (won't block build; we force the venv via $(PYTHON) anyway)
.PHONY: check-venv
check-venv:
ifndef SKIP_VENV_CHECK
	@cur="$${VIRTUAL_ENV:-}"; exp="$$(readlink -f $(VENV))"; \
	if [[ "$$cur" != "$$exp" ]]; then \
		echo ""; \
		echo "ℹ️  Not in expected venv. That's okay — this Makefile will still use:"; \
		echo "    $(VENV)/bin/python"; \
		echo "Current: $$cur"; \
		echo "Expected: $$exp"; \
		echo ""; \
	fi
endif

# Put it first so that "make" without argument is like "make help".
.PHONY: help
help: which-python
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: Makefile
%: Makefile check-venv which-python
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

# Build without running the gallery plots
.PHONY: html-noplot
html-noplot: check-venv which-python
	@$(SPHINXBUILD) -D plot_gallery=0 -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"
	@echo
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
