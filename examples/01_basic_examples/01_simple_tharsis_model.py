"""
Simple Geological Model - Tharsis Region
=========================================

This example creates a simple 3D geological model of the Tharsis region using GemPy.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import gempy as gp

# Set random seed for reproducibility
np.random.seed(1234)

# %%
# Import paths configuration
from mineye.config import paths

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

# Model resolution: use octree with refinement level 5
resolution = None  # Using octrees instead of regular grid
refinement = 5

# %%
# Get Data Paths
# --------------
# Load paths to structural data

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()

print(f"Orientations: {mod_or_path}")
print(f"Points: {mod_pts_path}")

# %%
# Create GemPy Geological Model
# ------------------------------
# Build the model structure with imported structural data

simple_geo_model = gp.create_geomodel(
    project_name='simple_model',
    extent=extent,
    refinement=refinement,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

# %%
# Map Geological Units
# --------------------
# Define the stratigraphic stack

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }
)

# %%
# Compute the Model
# -----------------
# Run the interpolation algorithm

gp.compute_model(simple_geo_model)
print("✓ Model computed successfully")

# %%
# Model with Topography
# ---------------------
# Add topographic surface from DEM

topography_path = paths.get_topography_path()
gp.set_topography_from_file(
    grid=simple_geo_model.grid,
    filepath=topography_path,
    crop_to_extent=[-695558, simple_geo_model.grid.extent[2],
                    simple_geo_model.grid.extent[1], simple_geo_model.grid.extent[3]]
)

# Recompute with topography
gp.compute_model(simple_geo_model)
print("✓ Model with topography computed successfully")
