# Sphinx Documentation Adaptation Guide for Mineye-Terranigma

This guide details all steps needed to adapt the GemPy Sphinx documentation configuration for the Mineye-Terranigma project.

## Overview

The current docs folder contains Sphinx configuration files copied from a GemPy project. These need to be customized to generate documentation for Mineye-Terranigma.

---

## 1. Critical Configuration Changes (`docs/source/conf.py`)

### 1.1 Project Information (Lines 104-116)

**Current:**
```python
project = 'GemPy'
copyright = u'2023-{}, Gempy Probability Developers'.format(year)
with open(os.path.join(os.path.dirname(__file__), '../../AUTHORS.rst'), 'r') as f:
    author = f.read()
version = gempy.__version__
release = gempy.__version__
```

**Action Required:**
- [ ] Change `project = 'Mineye-Terranigma'`
- [ ] Update copyright to your organization/name
- [ ] Create `AUTHORS.rst` file (see Section 2.1) OR remove the author file reading logic
- [ ] Update version/release to use your package version OR hardcode it

**Questions for you:**
- What should the copyright holder name be?
Terranigma Solutions GmbH
- Do you want to maintain an AUTHORS.rst file?
Maybe but then it would be the company
- Does `mineye` package have a `__version__` attribute?
not yet but it should. Can you make a setup.py based on the gempy one:

import os
from os import path

from setuptools import setup, find_packages


def read_requirements(file_name, base_path=""):
# Construct the full path to the requirements file
full_path = os.path.join(base_path, file_name)
requirements = []
with open(full_path, "r", encoding="utf-8") as f:
for line in f:
# Strip whitespace and ignore comments
line = line.strip()
if line.startswith("#") or not line:
continue

            # Handle -r directive
            if line.startswith("-r "):
                referenced_file = line.split()[1]  # Extract the file name
                # Recursively read the referenced file, making sure to include the base path
                requirements.extend(read_requirements(referenced_file, base_path=base_path))
            else:
                requirements.append(line)

    return requirements




with open("README.md", "r") as fh:
long_description = fh.read()

setup(
name='gempy',
packages=find_packages(exclude=('test', 'docs', 'examples')),
install_requires=read_requirements("requirements.txt", "requirements"),
extras_require={
"opt": read_requirements("optional-requirements.txt", "requirements"),
"base": read_requirements("base-requirements.txt", "requirements"),
},
url='https://github.com/cgre-aachen/gempy',
license='EUPL-1.2',
author='Miguel de la Varga, Alexander Zimmerman, Elisa Heim, Alexander Schaaf, Fabian Stamm, Florian Wellmann, Jan Niederau, Andrew Annex',
author_email='gempy@terranigma-solutions.com',
description='An Open-source, Python-based 3-D structural geological modeling software.',
long_description=long_description,
long_description_content_type='text/markdown',
keywords=['geology', '3-D modeling', 'structural geology', 'uncertainty'],
classifiers=[
"Development Status :: 5 - Production/Stable",
"Intended Audience :: Developers",
"License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)",
"Operating System :: OS Independent",
"Programming Language :: Python :: 3",
"Programming Language :: Python :: 3.10",
"Programming Language :: Python :: 3.11",
"Programming Language :: Python :: 3.12",
],
setup_requires=['setuptools_scm'],
use_scm_version={
"root"            : ".",
"relative_to"     : __file__,
"write_to"        : path.join("gempy", "_version.py"),
"fallback_version": "3.0.0"
},
)



### 1.2 Import Dependencies (Lines 20-29)

**Current:**
```python
import gempy
import pyvista
```

**Action Required:**
- [ ] Change `import gempy` to `import mineye` (or remove if not needed)
- [ ] Check if you need PyVista configuration (lines 31-38) - likely NOT needed if you're not generating 3D visualization plots in docs
- [ ] Remove or comment out PyVista section if not used

### 1.3 Sphinx Extensions (Lines 51-64)

**Current extensions are mostly fine, but verify:**
- [ ] `myst_parser` - for Markdown support (you have `devlog/ProbabilisticModeling.md`)
- [ ] `sphinx_gallery.gen_gallery` - only needed if you want auto-generated example galleries. Yes we do

**Action Required:**
- [ ] If you don't have example scripts to showcase, remove `'sphinx_gallery.gen_gallery'` from extensions
- [ ] Remove or comment out entire `sphinx_gallery_conf` section (lines 135-168) if not using galleries

### 1.4 Sphinx Gallery Configuration (Lines 135-168)

**Current:** Points to non-existent example directories

**Action Required - CRITICAL:**
- [ ] **Option A (Recommended):** Remove `'sphinx_gallery.gen_gallery'` from extensions AND delete the entire `sphinx_gallery_conf` dictionary
- [ ] **Option B:** Create example directories and scripts (see Section 5)

### 1.5 HTML Theme Configuration (Lines 178-204)

**Current:** Uses 'alabaster' theme with GemPy branding

**Action Required:**
- [ ] Update GitHub links:
  ```python
  'github_user': 'YOUR_GITHUB_USERNAME',
  'github_repo': 'Mineye-Terranigma',
  ```
- [ ] Update logo path or set to empty if no logo yet:
  ```python
  'logo': 'logos/mineye_logo.png',  # or remove this line
  ```
- [ ] Update `html_favicon` (line 204) or remove if no favicon yet

### 1.6 Makefile (`docs/Makefile`)

**Current:** Line 7: `SPHINXPROJ = gempy`

**Action Required:**
- [ ] Change to `SPHINXPROJ = mineye-terranigma` or `mineye`

---

## 2. Missing Files to Create

### 2.1 AUTHORS.rst (Referenced in conf.py line 109)

**Location:** `docs/AUTHORS.rst` (project root, two levels up from conf.py)

**Action Required:**
- [ ] Create file with content like:
  ```rst
  Miguel de la Varga
  [Add other contributors]
  ```
- [ ] OR modify conf.py to remove the file reading and hardcode author name

### 2.2 Example Directories (Referenced in conf.py lines 137-151)

**Current:** Points to non-existent directories:
- `../../examples/tutorials/0-intro`
- `../../examples/tutorials/1-first_example_of_inference`
- etc.

**Action Required:**
- [ ] **Recommended:** Delete the `sphinx_gallery_conf` section entirely if you don't need example galleries
- [ ] **Alternative:** Create these directories and add Python example scripts (see Section 5)

### 2.3 Logo and Favicon Files

**Current references:**
- `_static/logos/gempy.png` (line 183)
- `_static/logos/favicon.ico` (line 204)

**Existing in your project:**
- `_static/logos/` directory exists
- `_static/logos/logo_CGRE.png` (referenced in index.rst)
- `_static/logos/Terranigma.png` (referenced in index.rst)

**Action Required:**
- [ ] Choose which logo to use as the main sidebar logo
- [ ] Create a favicon.ico file OR remove the `html_favicon` line from conf.py
- [ ] Verify logo file names match what's in your `_static/logos/` directory

---

## 3. Content Files to Update

### 3.1 Main Index (`docs/source/index.rst`)

**Current:** Heavy GemPy branding and content

**Action Required:**
- [ ] Line 6: Update title from "GemPy Probability" to "Mineye-Terranigma"
- [ ] Lines 9-16: Rewrite overview section to describe your project
- [ ] Lines 25-39: Update or remove the toctree entries for example galleries
- [ ] Lines 49-92: Rewrite or remove the "Stochastic geological modeling" section
- [ ] Lines 94-98: Update references section
- [ ] Lines 106-110: Verify logo paths match your files

**Questions for you:**
- What is the main purpose/overview of Mineye-Terranigma?
- Should the docs include both work packages (Bayesian Segmentation + GemPy modeling)?

### 3.2 Installation Guide (`docs/source/installation.rst`)

**Current:** GemPy installation instructions

**Action Required:**
- [ ] Line 9-16: Update package name and PyPI instructions
- [ ] Lines 20-38: Remove or update Windows/MacOS specific sections for your dependencies
- [ ] Lines 45-75: Update developer installation section with your repo URL
- [ ] Lines 79-102: Update dependencies list to match your `requirements.txt`
- [ ] Lines 108-137: Remove or update conflicting packages section based on your stack

**Your current dependencies (from requirements.txt):**
- arviz, torch, gempy_probability, gempy_viewer, gempy, matplotlib, numpy, pandas, scipy, rasterio, scikit-image, scikit-learn, geopandas

**Action Required:**
- [ ] Simplify installation to just: `pip install -r requirements.txt`
- [ ] Document any special installation steps (e.g., PyTorch with CUDA)

### 3.3 API Reference (`docs/source/api_reference.rst`)

**Current:** Documents GemPy API

**Action Required - COMPLETE REWRITE:**
- [ ] Replace with your package structure. Example:
  ```rst
  Code
  ====

  .. toctree::
     :maxdepth: 3

  Mineye Package
  --------------
  .. currentmodule:: mineye
  .. autosummary::
      :toctree: Mineye API
      :template: base.rst

      [list your main public functions here]

  Submodules
  ----------
  .. currentmodule:: mineye.GeoModel
  .. autosummary::
      :toctree: GeoModel
      :template: base.rst

      [list submodule functions]
  ```

**Questions for you:**
- Which modules/functions should be publicly documented?
- Should it cover `mineye.GeoModel`, `mineye.GisHelpers`, `mineye.config`?

---

## 4. Optional Configuration

### 4.1 CNAME File (`docs/CNAME`)

**Current:** `gempy.rocks`

**Action Required:**
- [ ] **If hosting on GitHub Pages with custom domain:** Update to your domain
- [ ] **If using default GitHub Pages:** Delete this file
- [ ] **If not hosting online:** Delete this file

### 4.2 Remove Built Documentation

**Current:** Several built HTML files exist in `docs/` root

**Action Required:**
- [ ] Delete `docs/index.html`, `docs/.buildinfo`, `docs/main_docstrings` if present
- [ ] These will be regenerated when you build the docs

### 4.3 MyST Parser Configuration

**Current:** Enabled for Markdown support (you have `devlog/ProbabilisticModeling.md`)

**Action Required:**
- [ ] Verify line 99 includes `.md` files: `source_suffix = ['.rst', '.md']`
- [ ] If you want to include your devlog, add to index.rst toctree

---

## 5. Example Gallery Setup (Optional - Only if Desired)

If you want auto-generated example galleries like GemPy docs:

**Requirements:**
1. Create example script directories matching `examples_dirs` in conf.py
2. Write Python scripts (`.py` files) with special docstring format
3. Each script runs standalone and generates plots

**Directory structure needed:**
```
Mineye-Terranigma/
├── examples/
│   ├── basic_examples/
│   │   ├── plot_example_1.py
│   │   └── plot_example_2.py
│   └── advanced_examples/
│       └── plot_advanced_1.py
└── docs/
    └── source/
        └── ... (sphinx-gallery generates galleries here)
```

**Action Required:**
- [ ] Decide if you want example galleries
- [ ] If NO: Remove `sphinx_gallery.gen_gallery` extension and entire config
- [ ] If YES: Create example directories and scripts (this is a lot of work!)

---

## 6. Build and Test

### 6.1 Install Documentation Dependencies

**Action Required:**
```bash
pip install sphinx sphinx-gallery myst-parser
```

### 6.2 Build Documentation

**Action Required:**
```bash
cd docs
make html
```

**Expected output:** Built HTML in `docs/build/html/`

### 6.3 Common Build Errors and Fixes

**Error:** `FileNotFoundError: AUTHORS.rst`
- **Fix:** Create the file OR remove the file reading code from conf.py

**Error:** `ModuleNotFoundError: No module named 'gempy'`
- **Fix:** Remove `import gempy` and related code from conf.py

**Error:** `sphinx_gallery` errors about missing example directories
- **Fix:** Remove sphinx_gallery extension and configuration

**Error:** Missing logos/images
- **Fix:** Update paths in conf.py and index.rst to match your actual files

### 6.4 View Documentation Locally

**Action Required:**
```bash
# After successful build
cd build/html
python -m http.server 8000
# Open browser to http://localhost:8000
```

---

## 7. Information Needed from You

To complete the adaptation, please provide:

### Project Information
- [ ] **Project description:** What should the overview say about Mineye-Terranigma?
- [ ] **Copyright holder:** Who owns the copyright? (You, your organization, etc.)
- [ ] **Main features:** What are the key capabilities users should know about?
- [ ] **Target audience:** Who will use this documentation?

### Technical Details
- [ ] **Package version:** Should docs show version from `mineye.__version__` or hardcoded?
- [ ] **Public API:** Which modules/functions should be documented in API reference?
- [ ] **Example galleries:** Do you want auto-generated example galleries? (Recommended: NO for now)

### Hosting and Distribution
- [ ] **Documentation hosting:** GitHub Pages? ReadTheDocs? Local only?
- [ ] **Custom domain:** Do you have a custom domain for docs?
- [ ] **Target URL:** Where will the docs be accessible?

### Branding
- [ ] **Logo:** Do you have a logo file for the sidebar?
- [ ] **Favicon:** Do you have a favicon.ico file?
- [ ] **Color scheme:** Do you want to customize the Alabaster theme colors?

---

## 8. Quick Start Checklist (Minimum Viable Docs)

For a basic working documentation site, do these steps in order:

1. **conf.py modifications:**
   - [ ] Line 27: Comment out or remove `import gempy`
   - [ ] Line 28-38: Comment out entire PyVista section
   - [ ] Line 63: Remove `'sphinx_gallery.gen_gallery'` from extensions
   - [ ] Line 105: Change `project = 'Mineye-Terranigma'`
   - [ ] Line 107: Update copyright
   - [ ] Line 109-110: Comment out or remove AUTHORS.rst reading
   - [ ] Line 110: Replace with `author = 'Your Name'`
   - [ ] Line 112-113: Replace with hardcoded version or use `mineye.__version__` if it exists
   - [ ] Line 135-168: Delete entire `sphinx_gallery_conf` dictionary
   - [ ] Line 180-181: Update GitHub user/repo
   - [ ] Line 183: Remove logo line or update path

2. **Makefile:**
   - [ ] Line 7: Change `SPHINXPROJ = mineye-terranigma`

3. **index.rst:**
   - [ ] Update title and overview
   - [ ] Remove example gallery toctrees (lines 31-39)
   - [ ] Simplify to just show installation and API reference

4. **installation.rst:**
   - [ ] Update to show: `pip install -r requirements.txt`

5. **api_reference.rst:**
   - [ ] Simplify to just document main `mineye` module

6. **Build test:**
   - [ ] Run `make html` and fix any errors

---

## 9. Advanced Configuration (Future Enhancements)

Once basic docs are working, consider:

- [ ] Switch to modern theme (e.g., `sphinx_rtd_theme`, `furo`, `sphinx-book-theme`)
- [ ] Add more content pages (tutorials, theory, examples)
- [ ] Set up automatic building on GitHub Actions
- [ ] Configure ReadTheDocs hosting
- [ ] Add version switching for documentation
- [ ] Create Jupyter notebook examples (with `nbsphinx`)
- [ ] Add changelog/release notes page

---

## 10. Additional Resources

- **Sphinx Documentation:** https://www.sphinx-doc.org/
- **Alabaster Theme:** https://alabaster.readthedocs.io/
- **Sphinx Gallery:** https://sphinx-gallery.github.io/
- **MyST Parser (Markdown):** https://myst-parser.readthedocs.io/

---

## Summary

The main adaptation work involves:
1. Removing GemPy-specific imports and configuration
2. Updating project metadata (name, copyright, author)
3. Rewriting content files (index, installation, API reference)
4. Deciding whether to use Sphinx Gallery for examples
5. Testing the build process

**Estimated time:** 2-4 hours for basic setup, more if adding extensive content and examples.

**Recommended approach:** Start with Quick Start Checklist (Section 8) to get a minimal working version, then iterate to add more features.
