"""
Gravity Forward Modeling and Data Comparison
=============================================

This example demonstrates gravity forward modeling using a GemPy geological model.

The workflow includes:

* Loading a pre-built geological model
* Reading observed gravity measurements  
* Computing forward gravity response
* Comparing modeled vs observed data

.. note::
   This example requires the mineye package and data files.
   Install with: ``pip install -e .`` from project root
"""

# %%
# Workflow Overview
# -----------------

import sys
print("This is a placeholder example.")
print("Full probabilistic modeling examples require:")
print("  - mineye package installed (pip install -e .)")
print("  - Geological model data")
print("  - Gravity measurements data")
print("\nSee examples/geomodels/SimpleGravimetricResponse.py for the full script")

# %%
# The full example demonstrates:
#
# * Loading geological models and geophysical data
# * Configuring centered grids for gravity computation
# * Computing forward gravity response
# * Comparing modeled and observed data
# * Visualizing residuals for model refinement
#
# The residuals can be used in probabilistic inversion frameworks
# for uncertainty quantification.

