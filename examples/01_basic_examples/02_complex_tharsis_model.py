"""
Complex Geological Model - Tharsis Region
=========================================

This example creates a complex 3D geological model of the Tharsis region using GemPy.
It features an erosive Tournaisian Plutonite unit overlying a conformable Devonian sequence.
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

# Use octree refinement
resolution = None
refinement = 5

# %%
# Get Data Paths
# --------------
# Load paths to structural data

mod_or_path = paths.get_orientations_path_complex()
mod_pts_path = paths.get_points_path_complex()

print(f"Orientations: {mod_or_path}")
print(f"Points: {mod_pts_path}")

# %%
# Create Stratigraphic Stack Model
# ------------------------------
from mineye.GeoModel import helper_plotter

stratigraphic_geo_model = gp.create_geomodel(
    project_name='stratigraphic_stack_model',
    extent=extent,
    refinement=refinement,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=os.path.join(BASE_DIR, "points_orientations_incl_topo_all_df.csv"),
        path_to_surface_points=os.path.join(BASE_DIR, "points_incl_topo_all_df.csv"),
    )
)

gp.map_stack_to_surfaces(
    gempy_model=stratigraphic_geo_model,
    mapping_object={
        "Strat_Series1": ("Upper Carboniferous Volcanics","Mid Carboniferous Shales", "Visean Shales", "Upper Devonian Siliciclastics")
    }
)

topography_path = paths.get_topography_path()
gp.set_topography_from_file(
    grid=stratigraphic_geo_model.grid,
    filepath=topography_path,
    crop_to_extent=[-695558, stratigraphic_geo_model.grid.extent[2],
                    stratigraphic_geo_model.grid.extent[1], stratigraphic_geo_model.grid.extent[3]]
)

gempy_model = gp.compute_model(stratigraphic_geo_model)
helper_plotter.create_cross_section(stratigraphic_geo_model, cross_section=5, vertical_exaggeration=10)


# %%
# Create GemPy Geological Model
# ------------------------------

geo_model = gp.create_geomodel(
    project_name='complex_tharsis_model',
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

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Tournaisian_Series": ["Tournaisian Plutonites"],
        "Devonian_Series": [
            "Famennian Siliciclastics"
        ]
    }
)

# Verify stack order
print(geo_model.structural_frame)

# %%
# Model with Topography
# ---------------------
# Add topographic surface from DEM

topography_path = paths.get_topography_path()
gp.set_topography_from_file(
    grid=geo_model.grid,
    filepath=topography_path,
    crop_to_extent=[-695558, geo_model.grid.extent[2],
                    geo_model.grid.extent[1], geo_model.grid.extent[3]]
)

# %%
# Compute the Model
# -----------------

gp.compute_model(geo_model)
print("✓ Model computed successfully")

# %%
# Visualize the Model
# -------------------

import gempy_viewer as gpv
gpv.plot_3d(geo_model, ve=5, image=False)

