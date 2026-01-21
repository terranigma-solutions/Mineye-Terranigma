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
import torch
from matplotlib import pyplot as plt

import gempy as gp
import geopandas as gpd

from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize

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

observed_gravity_ugal= observed_gravity * 1000
norm_params = normalize(
    sol.gravity,
    observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)
grav = align_forward_to_observed(sol.gravity, norm_params)

print(f"\nComputed gravity values:")
print(f"  Shape: {grav.shape}")
print(f"  Range: {grav.min():.2f} to {grav.max():.2f} mGal")
print(f"  Mean: {grav.mean():.2f} mGal")

# %%
# Compare with Observations
# -------------------------
# Calculate residuals between observed and computed gravity

if isinstance(grav, torch.Tensor):
    grav = grav.detach().numpy()
    
residuals = observed_gravity_ugal - grav

print(f"\nGravity residuals:")
print(f"  Mean: {residuals.mean():.2f} mGal")
print(f"  Std: {residuals.std():.2f} mGal")
print(f"  RMS: {np.sqrt(np.mean(residuals**2)):.2f} mGal")

# %%
# Visualize the Model
# -------------------
# Create 2D and 3D visualizations of the geological model

import gempy_viewer as gpv

gpv.plot_2d(simple_geo_model)
gpv.plot_3d(simple_geo_model, ve=5, image=True)

# %%
# Visualize Gravity Results
# --------------------------
# Plot forward gravity and comparison with observations

from mineye.GeoModel.plotting.probabilistic_analysis import plot_fw_geophysics, plot_comparison

plot_fw_geophysics(grav, gravity_data, xy_ravel)
plot_comparison(observed_gravity, grav, xy_ravel)

# %%
# Visualize Forward Model Results
# --------------------------------
grav_forward = grav
print(f"✓ Forward model computed")
print(f"  Gravity range: [{grav_forward.min():.2f}, {grav_forward.max():.2f}] µGal")

fig, ax = plt.subplots(figsize=(10, 8))

scatter = ax.scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=grav_forward, s=50,
    cmap='viridis_r', alpha=0.8,
    edgecolors='black', linewidth=0.5
)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label(r'Forward Model Gravity (µGal)', fontsize=12)

ax.set_title('Forward Gravity Model', fontsize=14, fontweight='bold')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %%
# Compare with Observed Gravity
# ------------------------------

# Convert observed from mGal to µGal
observed_ugal = observed_gravity_ugal

# Create comparison plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Observed
sc1 = axes[0].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=observed_ugal, s=50, cmap='viridis_r',
    edgecolors='black', linewidth=0.5
)
axes[0].set_title('Observed Gravity')
axes[0].set_xlabel('X (m)')
axes[0].set_ylabel('Y (m)')
plt.colorbar(sc1, ax=axes[0], label='µGal')

# Forward model
sc2 = axes[1].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=grav_forward, s=50, cmap='viridis_r',
    edgecolors='black', linewidth=0.5
)
axes[1].set_title('Forward Model')
axes[1].set_xlabel('X (m)')
plt.colorbar(sc2, ax=axes[1], label='µGal')

# Residuals
sc3 = axes[2].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=residuals, s=50, cmap='RdBu_r',
    edgecolors='black', linewidth=0.5
)
axes[2].set_title('Residuals (Obs - Model)')
axes[2].set_xlabel('X (m)')
plt.colorbar(sc3, ax=axes[2], label='µGal')

plt.tight_layout()
plt.show()

# %%
# Correlation Analysis
# --------------------

fig, ax = plt.subplots(figsize=(8, 8))

ax.scatter(observed_ugal, grav_forward, alpha=0.6, s=60,
           edgecolors='black', linewidth=0.5)

# 1:1 line
lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
]
ax.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')

# Calculate statistics
correlation = np.corrcoef(observed_ugal, grav_forward)[0, 1]
rmse = np.sqrt(np.mean(residuals**2))

ax.set_xlabel('Observed (µGal)', fontsize=12)
ax.set_ylabel('Forward Model (µGal)', fontsize=12)
ax.set_title('Observed vs Modeled Gravity', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Add statistics text box
textstr = f'R = {correlation:.3f}\nR² = {correlation**2:.3f}\nRMSE = {rmse:.2f} µGal'
ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
        fontsize=11)

plt.tight_layout()
plt.show()

# %%
# Print Summary Statistics
# -------------------------

print("\n" + "="*50)
print("GRAVITY COMPARISON STATISTICS")
print("="*50)
print(f"Number of measurements: {len(observed_ugal)}")
print(f"\nObserved gravity:")
print(f"  Mean: {observed_ugal.mean():.2f} µGal")
print(f"  Std:  {observed_ugal.std():.2f} µGal")
print(f"  Range: [{observed_ugal.min():.2f}, {observed_ugal.max():.2f}] µGal")
print(f"\nForward model:")
print(f"  Mean: {grav_forward.mean():.2f} µGal")
print(f"  Std:  {grav_forward.std():.2f} µGal")
print(f"  Range: [{grav_forward.min():.2f}, {grav_forward.max():.2f}] µGal")
print(f"\nResiduals:")
print(f"  Mean: {residuals.mean():.2f} µGal")
print(f"  Std:  {residuals.std():.2f} µGal")
print(f"  RMSE: {rmse:.2f} µGal")
print(f"  MAE:  {np.abs(residuals).mean():.2f} µGal")
print(f"\nCorrelation: {correlation:.4f} (R² = {correlation**2:.4f})")
print("="*50)

# %%
# Summary
# -------
#
# This example demonstrated:
#
# * Loading geological models and gravity data
# * Setting up centered grids for gravity computation
# * Computing forward gravity response
# * Comparing modeled with observed data
# * Statistical analysis of model fit
#
# **Next steps:**
#
# * Use residuals for probabilistic inversion
# * Uncertainty quantification with Bayesian methods
# * Joint inversion with multiple data types
# sphinx_gallery_thumbnail_number = 4