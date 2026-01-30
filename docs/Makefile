
# Minimal makefile for Sphinx documentation

SHELL        := /usr/bin/env bash

# CONFIG
# If VENV is not set, try to use the active virtual environment or just 'python'
VENV         ?= $(shell which python | sed 's/\/bin\/python//')
PYTHON       := $(shell if [ -x "$(VENV)/bin/python" ]; then echo "$(VENV)/bin/python"; else echo "python"; fi)
SPHINXBUILD  := $(PYTHON) -m sphinx       # ← module form guarantees venv python
SPHINXOPTS   :=
SPHINXPROJ   := mineye-terranigma
SOURCEDIR    := source
BUILDDIR     := build

# Sanity check: print which Python will be used (helpful in logs)
.PHONY: which-python
which-python:
	@echo "Using Python: $(PYTHON)"
	@$(PYTHON) --version

# Optional venv check
.PHONY: check-venv
check-venv:
ifndef SKIP_VENV_CHECK
	@echo "ℹ️  Using Python from: $(PYTHON)"
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
