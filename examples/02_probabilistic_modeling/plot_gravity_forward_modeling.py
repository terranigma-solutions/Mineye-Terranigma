"""
Gravity Forward Modeling and Data Comparison
=============================================

This example demonstrates gravity forward modeling using a GemPy geological model.
We compute the gravity response of the model and compare it with observed gravity data.

The workflow includes:

* Loading a pre-built geological model
* Reading observed gravity measurements
* Computing forward gravity response
* Comparing modeled vs observed data
* Visualizing the results
"""

# %%
# Setup and Configuration
# -----------------------

import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

import numpy as np
import gempy as gp
import geopandas as gpd
import matplotlib.pyplot as plt
import pickle
from mineye.config import paths

# Configure paths
BASE_DIR = paths.get_base_dir()
geophysical_dir = paths.get_geophysical_dir(BASE_DIR)
pickle_model_path = paths.get_pickle_model_path(BASE_DIR)

# Model parameters
density_plutonites = 2.9  # g/cm³
density_sedimentary_host = 2.3  # g/cm³
gravity_resolution = 20

# %%
# Load Geological Model and Gravity Data
# ---------------------------------------

try:
    with open(pickle_model_path, 'rb') as f:
        geo_model = pickle.load(f)

    gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
    observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

    print(f"✓ Loaded geological model: {geo_model.meta.project_name}")
    print(f"✓ Loaded {len(observed_gravity)} gravity observations")

    has_data = True
except FileNotFoundError as e:
    print(f"⚠ Data files not found: {e}")
    print("This example requires pre-computed model and gravity data.")
    has_data = False

# %%
# Extract Measurement Locations
# ------------------------------

if has_data:
    xy_coords = gravity_data[['geometry']].apply(
        lambda row: (row['geometry'].x, row['geometry'].y), axis=1
    )
    xy_ravel = np.array([[x, y] for x, y in xy_coords])

    # Add custom grid for gravity measurements
    geo_model.grid.create_custom_grid(xy_ravel=xy_ravel)

    print(f"✓ Created measurement grid with {len(xy_ravel)} points")

# %%
# Configure Geophysical Forward Modeling
# ---------------------------------------

if has_data:
    gp.set_custom_grid(geo_model.grid, xyz_coord=xy_ravel)
    gp.set_centered_grid(
        geo_model.grid,
        centers=xy_ravel,
        resolution=np.array([10, 10, 15]),
        radius=np.array([5000, 5000, 5000])
    )

# %%
# Compute Forward Gravity Response
# ---------------------------------

if has_data:
    # Enable gravity computation
    geo_model.geophysics_input.tz = geo_model.grid.centered_grid

    # Compute model
    sol = gp.compute_model(geo_model)

    # Extract gravity response
    simulated_gravity = sol.gravity

    print(f"✓ Computed gravity forward model")
    print(f"  Gravity range: [{simulated_gravity.min():.2f}, {simulated_gravity.max():.2f}] mGal")

# %%
# Compare Observed vs Modeled Gravity
# ------------------------------------

if has_data:
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # Observed gravity
    scatter1 = axes[0].scatter(
        xy_ravel[:, 0], xy_ravel[:, 1],
        c=observed_gravity, cmap='viridis', s=50
    )
    axes[0].set_title('Observed Gravity')
    axes[0].set_xlabel('X (m)')
    axes[0].set_ylabel('Y (m)')
    plt.colorbar(scatter1, ax=axes[0], label='mGal')

    # Modeled gravity
    scatter2 = axes[1].scatter(
        xy_ravel[:, 0], xy_ravel[:, 1],
        c=-simulated_gravity, cmap='viridis', s=50
    )
    axes[1].set_title('Modeled Gravity')
    axes[1].set_xlabel('X (m)')
    plt.colorbar(scatter2, ax=axes[1], label='mGal')

    # Residuals
    residuals = observed_gravity - (-simulated_gravity)
    scatter3 = axes[2].scatter(
        xy_ravel[:, 0], xy_ravel[:, 1],
        c=residuals, cmap='RdBu', s=50
    )
    axes[2].set_title('Residuals (Obs - Model)')
    axes[2].set_xlabel('X (m)')
    plt.colorbar(scatter3, ax=axes[2], label='mGal')

    plt.tight_layout()
    plt.show()

    # Print statistics
    print(f"\nGravity Statistics:")
    print(f"  Observed mean: {observed_gravity.mean():.2f} mGal")
    print(f"  Modeled mean: {-simulated_gravity.mean():.2f} mGal")
    print(f"  RMSE: {np.sqrt(np.mean(residuals**2)):.2f} mGal")

# %%
# Summary
# -------
#
# This example demonstrated:
#
# * Loading geological models and geophysical data
# * Configuring centered grids for gravity computation
# * Computing forward gravity response
# * Comparing modeled and observed data
#
# The residuals can be used to refine the geological model in
# probabilistic inversion frameworks.

if has_data:
    print("\n✓ Gravity forward modeling complete!")
else:
    print("\n⚠ Example requires data files - see installation guide for setup")
