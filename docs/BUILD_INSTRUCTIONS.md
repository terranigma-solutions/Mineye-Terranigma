# Building Sphinx Documentation

## Prerequisites

The Makefile is configured to check that you're using the correct virtual environment before building.

**Required virtual environment:** `~/.venv/2025`

## Setup Instructions

### 1. Activate the Virtual Environment

```bash
source ~/.venv/2025/bin/activate
```

### 2. Install Documentation Dependencies

```bash
cd /home/leguark/PycharmProjects/Mineye-Terranigma
pip install sphinx sphinx-gallery myst-parser matplotlib
```

**Note:** If you encounter Python 3.13 compatibility issues with Sphinx, you may need to:
```bash
pip install --upgrade sphinx
```

Or use a Python 3.10-3.12 environment instead.

### 3. Install the Mineye Package

To enable version info in documentation:

```bash
pip install -e .
```

## Building the Documentation

### Standard Build (with venv check)

```bash
cd docs
make html
```

The Makefile will automatically check if you're in the correct virtual environment.
If not, you'll see:

```
⚠️  Wrong virtual environment!

Current: <your current venv>
Expected: /home/leguark/.venv/2025

Please activate the correct environment:
  source ~/.venv/2025/bin/activate
```

### Build Without Gallery Execution

If you want to build without running example scripts:

```bash
cd docs
make html-noplot
```

### Clean Build

To start fresh:

```bash
cd docs
make clean
make html
```

## Viewing the Documentation

After successful build:

```bash
cd docs/build/html
python -m http.server 8000
```

Then open http://localhost:8000 in your browser.

## Common Issues

### Issue: Wrong Virtual Environment

**Error:**
```
⚠️  Wrong virtual environment!
```

**Solution:**
```bash
source ~/.venv/2025/bin/activate
```

### Issue: Missing Sphinx

**Error:**
```
sphinx-build: command not found
```

**Solution:**
```bash
pip install sphinx sphinx-gallery myst-parser matplotlib
```

### Issue: Python 3.13 Compatibility

**Error:**
```
ImportError: cannot import name 'Union' from 'types'
```

**Solution:**
Use Python 3.10-3.12, or upgrade sphinx:
```bash
pip install --upgrade sphinx>=8.0
```

### Issue: Missing Mineye Package

**Warning:**
```
Warning: mineye package not found
```

**Solution:**
```bash
cd /home/leguark/PycharmProjects/Mineye-Terranigma
pip install -e .
```

This enables version info but is not required for building docs.

## Makefile Targets

- `make help` - Show available targets
- `make html` - Build HTML documentation (checks venv)
- `make html-noplot` - Build without running example galleries (checks venv)
- `make clean` - Remove build directory
- `make check-venv` - Check if correct venv is activated

## Configuration

To change the required virtual environment, edit `docs/Makefile`:

```makefile
VENV = ~/.venv/2025  # Change this line
```

## Troubleshooting Build Failures

1. Check you're in the right venv: `echo $VIRTUAL_ENV`
2. Check sphinx is installed: `sphinx-build --version`
3. Check Python version: `python --version` (should be 3.10-3.12)
4. Try clean build: `make clean && make html`
5. Check full error log in `/tmp/sphinx-err-*.log`

## Automated Builds

For CI/CD or automated builds, you can bypass the venv check by setting:

```bash
export SKIP_VENV_CHECK=1
make html
```

(Note: This feature would need to be added to the Makefile if needed)
