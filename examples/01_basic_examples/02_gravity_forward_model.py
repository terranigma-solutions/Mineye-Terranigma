"""
Forward Gravity Modeling - Tharsis Region
==========================================

This example demonstrates forward gravity modeling using a 3D geological model.
We compute the gravity response of the geological structure at real measurement locations.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import gempy as gp
import geopandas as gpd

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
# Load paths to structural and geophysical data

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
gravity_data_path = paths.get_gravity_data_path()

print(f"Orientations: {mod_or_path}")
print(f"Points: {mod_pts_path}")
print(f"Gravity data: {gravity_data_path}")

# %%
# Create GemPy Geological Model
# ------------------------------
# Build the model structure with imported structural data

simple_geo_model = gp.create_geomodel(
    project_name='gravity_model',
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
# Load Gravity Observations
# --------------------------
# Read actual gravity measurements from the field

gravity_data = gpd.read_file(gravity_data_path)
observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

print(f"Number of gravity observations: {len(observed_gravity)}")
print(f"Gravity range: {observed_gravity.min():.2f} to {observed_gravity.max():.2f} mGal")

# %%
# Set Up Gravity Measurement Locations
# -------------------------------------
# Create measurement grid at actual device locations

xy_ravel = np.column_stack([
    np.array(gravity_data.geometry.x.values),
    np.array(gravity_data.geometry.y.values),
    np.full(len(gravity_data), 0)  # Set Z to surface elevation
])

print(f"Using {len(xy_ravel)} actual measurement points")

# %%
# Configure Centered Grid for Gravity Computation
# ------------------------------------------------
# Set up centered grid around each measurement point

gp.set_centered_grid(
    grid=simple_geo_model.grid,
    centers=xy_ravel,
    resolution=np.array([10, 10, 15]),
    radius=np.array([5000, 5000, 5000])
)

# %%
# Calculate Gravity Gradient
# ---------------------------
# Compute the gravity gradient (tz component) for the centered grid

gravity_gradient = gp.calculate_gravity_gradient(simple_geo_model.grid.centered_grid)
print(f"Gravity gradient tensor shape: {gravity_gradient.shape}")

# %%
# Set Density Values
# ------------------
# Define densities for geological units

density_plutonites = 2.9  # kg/m³
density_sedimentary_host = 2.3  # kg/m³

simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
    tz=gravity_gradient,
    densities=np.array([density_plutonites, density_sedimentary_host])
)

print(f"Plutonite density: {density_plutonites} g/cm³")
print(f"Sedimentary host density: {density_sedimentary_host} g/cm³")

# %%
# Compute Forward Gravity Model
# ------------------------------
# Run the interpolation and gravity computation

simple_geo_model.interpolation_options.mesh_extraction = False
sol = gp.compute_model(simple_geo_model)

print("✓ Forward gravity model computed successfully")

# %%
# Extract Gravity Results
# -----------------------
# Get the computed gravity response

grav = sol.gravity
print(f"\nComputed gravity values:")
print(f"  Shape: {grav.shape}")
print(f"  Range: {grav.min():.2f} to {grav.max():.2f} mGal")
print(f"  Mean: {grav.mean():.2f} mGal")

# %%
# Compare with Observations
# -------------------------
# Calculate residuals between observed and computed gravity

residuals = observed_gravity - grav
print(f"\nGravity residuals:")
print(f"  Mean: {residuals.mean():.2f} mGal")
print(f"  Std: {residuals.std():.2f} mGal")
print(f"  RMS: {np.sqrt(np.mean(residuals**2)):.2f} mGal")
