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
import matplotlib.pyplot as plt

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

mod_or_sed_path = paths.get_orientations_path_sed_complex()
mod_pts_sed_path = paths.get_points_path_sed_complex()

mod_or_plut_path = paths.get_orientations_path_magmatic_complex()
mod_pts_plut_path = paths.get_points_path_magmatic_complex()

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
    resolution=[64, 64, 64],  # Need regular grid for cross-sections and merging
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_sed_path,
        path_to_surface_points=mod_pts_sed_path,
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
# Adding Plutonite Unit
# ------------------------------
# Create a separate model for the plutonite body, then merge with stratigraphic model

plutonite_id = 1

# Need resolution set to get a regular grid for merging
plutonite_geo_model = gp.create_geomodel(
    project_name='plutonite_model',
    extent=extent,
    refinement=refinement,
    resolution=[64, 64, 64],  # Need regular grid for merging
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_plut_path,
        path_to_surface_points=mod_pts_plut_path,
    )
)

# Map the plutonite surface
gp.map_stack_to_surfaces(
    gempy_model=plutonite_geo_model,
    mapping_object={
        "Plutonite_Series": ["Tournaisian Plutonites"]
    }
)

# Compute the plutonite model
gp.compute_model(plutonite_geo_model)
print("Plutonite model computed")

# Get plutonite lithology block
plut_lith_block = plutonite_geo_model.solutions.raw_arrays.lith_block
plut_lith_block_reshaped = plut_lith_block.reshape(64, 64, 64)
plutonite_mask = plut_lith_block_reshaped == plutonite_id

# The interpolated grid of formation IDs is stored here
lith_block = stratigraphic_geo_model.solutions.raw_arrays.lith_block
lith_block_reshaped = lith_block.reshape(64, 64, 64)

# Get the coordinates of these voxels
voxel_coords = stratigraphic_geo_model.grid.regular_grid.values
voxel_coords_reshaped = voxel_coords.reshape(64, 64, 64, 3)

# Flip Y coordinate
y_min, y_max = np.min(voxel_coords[:, 1]), np.max(voxel_coords[:, 1])
voxel_coords_flipped = voxel_coords.copy()
voxel_coords_flipped[:, 1] = y_max - (voxel_coords[:, 1] - y_min)

# Insert the igneous rock locations into the stratigraphic model
lith_block_reshaped[plutonite_mask] = 6
lith_block_modified = lith_block_reshaped.flatten()

# Plot the combined model
sample_step = 10  # Take every 10th point

# Create interactive 3D plot
fig = plt.figure(figsize=(12, 9))
ax = fig.add_subplot(111, projection='3d')

# Create the scatter plot with color mapping based on Z values
scatter = ax.scatter(voxel_coords_flipped[:,0][::sample_step],
                    voxel_coords_flipped[:,1][::sample_step],
                    voxel_coords_flipped[:,2][::sample_step],
                    c=lith_block_modified[::sample_step],
                    cmap='viridis',
                    s=10,
                    alpha=0.6)

# Add colorbar
plt.colorbar(scatter, ax=ax, label='Rock Type', shrink=0.5)

# Set labels and title
ax.set_xlabel('X Coordinate')
ax.set_ylabel('Y Coordinate')
ax.set_zlabel('Z Elevation')

# Improve the view
ax.view_init(elev=90, azim=270)
plt.tight_layout()
plt.show()