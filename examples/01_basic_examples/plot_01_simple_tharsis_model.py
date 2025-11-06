"""
Simple 3D Geological Model - Tharsis Region
============================================

This example creates a 3D geological model of the Tharsis region in Spain using GemPy.
The model includes plutonitic intrusions and sedimentary host rocks.

**Key features:**

* Loading structural data (orientations and contact points)
* Defining model extent and resolution
* Creating and computing a GemPy geological model
* Visualizing the model in 2D cross-sections
* Adding topography from a DEM file

.. note::
   This example requires the mineye package and data files to be present.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(1234)

# Configure matplotlib for non-interactive environments
matplotlib.use('Agg')

try:
    import gempy as gp
    import gempy_viewer as gpv
    GEMPY_AVAILABLE = True
except ImportError:
    GEMPY_AVAILABLE = False
    print("⚠ GemPy not installed or missing dependencies (torch)")

try:
    from mineye.GeoModel import Plotter
    from mineye.config import paths
    MINEYE_AVAILABLE = True
except ImportError:
    MINEYE_AVAILABLE = False
    print("⚠ Mineye package not installed")
    print("Install with: pip install -e . from project root")

# %%
# Configuration and Paths
# -----------------------
# Set up paths to input data files

if MINEYE_AVAILABLE and GEMPY_AVAILABLE:
    try:
        topography_path = paths.get_topography_path()
        mod_or_path = paths.get_orientations_path()
        mod_pts_path = paths.get_points_path()
        DATA_AVAILABLE = True
    except Exception as e:
        DATA_AVAILABLE = False
        print(f"⚠ Data files not accessible: {e}")
else:
    DATA_AVAILABLE = False

# %%
# Define Model Extent and Resolution
# -----------------------------------
# The model covers the Tharsis mining district

min_x = -709521
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
extent = [min_x, max_x, min_y, max_y, model_depth, max_z]

# Model resolution: 64x64x64 voxels
nx = ny = nz = 64
resolution = [nx, ny, nz]

print(f"Model extent: {extent}")
print(f"Model resolution: {resolution}")

# %%
# Build GemPy Geological Model
# -----------------------------
# Create the model with imported structural data

if DATA_AVAILABLE:
    simple_geo_model = gp.create_geomodel(
        project_name='simple_model',
        extent=extent,
        refinement=4,
        resolution=resolution,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=mod_or_path,
            path_to_surface_points=mod_pts_path,
        )
    )

    # Map geological units to surfaces
    gp.map_stack_to_surfaces(
        gempy_model=simple_geo_model,
        mapping_object={
            "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
        }
    )

    # Add topography from DEM
    gp.set_topography_from_file(
        grid=simple_geo_model.grid,
        filepath=topography_path
    )

    print("✓ Model structure created")
else:
    print("⚠ Skipping model creation - data not available")

# %%
# Compute the Geological Model
# -----------------------------
# Run the interpolation algorithm

if DATA_AVAILABLE:
    gp.compute_model(simple_geo_model)
    print("✓ Model computed successfully")

# %%
# Visualization: 2D Geological Map
# ---------------------------------
# Create a top-down geological map with topography

if DATA_AVAILABLE:
    p = gpv.plot_2d(
        simple_geo_model,
        section_names=['topography'],
        show_topography=True,
        show_lith=True,
        show_boundaries=True,
        show_data=True,
        legend=False
    )
    print("✓ 2D geological map created")
    plt.close('all')  # Clean up matplotlib figures

# %%
# Visualization: Cross-Section
# -----------------------------
# Create a vertical cross-section through the model

if DATA_AVAILABLE:
    Plotter.create_cross_section(simple_geo_model, cross_section=5)
    print("✓ Cross-section created")
    plt.close('all')  # Clean up matplotlib figures

# %%
# Summary
# -------
#
# This example demonstrated:
#
# * Creating a 3D geological model from structural data
# * Defining model extent and resolution
# * Adding topographic constraints
# * Visualizing the model in 2D
#
# The model can be used for:
#
# * Forward modeling of geophysical responses
# * Probabilistic inversion and uncertainty quantification
# * Resource estimation
#
# See the full script at: ``examples/geomodels/Simple_Model_Tharsis.py``

if not DATA_AVAILABLE:
    print("\n" + "="*60)
    print("To run this example with real data:")
    print("1. Install mineye: pip install -e .")
    print("2. Ensure data files are in examples/Data/")
    print("3. Run the full script: examples/geomodels/Simple_Model_Tharsis.py")
    print("="*60)
