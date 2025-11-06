"""
Gravity Forward Modeling with GemPy
====================================

This example demonstrates forward modeling of gravity data using a 3D geological model.
We compute the gravity response at measurement locations and compare with observed data.

**Workflow:**

1. Load a pre-computed geological model
2. Read observed gravity measurements
3. Set up centered grid for gravity computation
4. Compute forward gravity response
5. Compare modeled vs observed gravity
6. Analyze residuals and statistics

.. note::
   This example requires pre-computed model data and gravity measurements.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import matplotlib.pyplot as plt
import pickle

try:
    import gempy as gp
    import geopandas as gpd
    GEMPY_AVAILABLE = True
except ImportError:
    GEMPY_AVAILABLE = False
    print("⚠ GemPy/geopandas not installed or missing dependencies")

try:
    from mineye.config import paths
    MINEYE_AVAILABLE = True
except ImportError:
    MINEYE_AVAILABLE = False
    print("⚠ Mineye package not installed")

# %%
# Load Geological Model and Gravity Data
# ---------------------------------------

if MINEYE_AVAILABLE and GEMPY_AVAILABLE:
    try:
        BASE_DIR = paths.get_base_dir()
        geophysical_dir = paths.get_geophysical_dir(BASE_DIR)
        pickle_model_path = paths.get_pickle_model_path(BASE_DIR)

        # Load pre-computed geological model
        with open(pickle_model_path, 'rb') as f:
            geo_model = pickle.load(f)

        # Load observed gravity measurements
        gravity_data = gpd.read_file(
            f"{geophysical_dir}/cleaned_gravity_data.geojson"
        )
        observed_gravity = gravity_data['VALU_BOU267'].values  # mGal

        DATA_AVAILABLE = True
        print(f"✓ Loaded model: {geo_model.meta.project_name}")
        print(f"✓ Loaded {len(observed_gravity)} gravity observations")

    except Exception as e:
        DATA_AVAILABLE = False
        print(f"⚠ Data not accessible: {e}")
else:
    DATA_AVAILABLE = False

# %%
# Extract Measurement Locations
# ------------------------------
# Get XYZ coordinates of gravity measurement devices

if DATA_AVAILABLE:
    extent = geo_model.grid.regular_grid.extent
    min_x, max_x, min_y, max_y, min_z, max_z = extent

    # Extract measurement point coordinates
    xy_ravel = np.column_stack([
        np.array(gravity_data.geometry.x.values),
        np.array(gravity_data.geometry.y.values),
        np.full(len(gravity_data), max_z)  # Surface elevation
    ])

    print(f"✓ Using {len(xy_ravel)} measurement locations")
    print(f"  X range: [{xy_ravel[:, 0].min():.0f}, {xy_ravel[:, 0].max():.0f}]")
    print(f"  Y range: [{xy_ravel[:, 1].min():.0f}, {xy_ravel[:, 1].max():.0f}]")

# %%
# Configure Centered Grid for Gravity
# ------------------------------------
# Set up grid around each measurement point

if DATA_AVAILABLE:
    gp.set_centered_grid(
        grid=geo_model.grid,
        centers=xy_ravel,
        resolution=np.array([10, 10, 15]),
        radius=np.array([5000, 5000, 5000])
    )

    print("✓ Centered grid configured")
    print(f"  Resolution: 10x10x15 voxels per station")
    print(f"  Radius: 5km around each measurement")

# %%
# Compute Forward Gravity Model
# ------------------------------

if DATA_AVAILABLE:
    # Calculate gravity gradient (Tz component)
    gravity_gradient = gp.calculate_gravity_gradient(
        geo_model.grid.centered_grid
    )

    # Set densities: [sedimentary, plutonites]
    density_sedimentary = 2.3  # g/cm³
    density_plutonites = 2.9  # g/cm³

    geo_model.geophysics_input = gp.data.GeophysicsInput(
        tz=gravity_gradient,
        densities=np.array([density_sedimentary, density_plutonites])
    )

    # Compute forward model
    geo_model.interpolation_options.mesh_extraction = False
    sol = gp.compute_model(
        gempy_model=geo_model,
        engine_config=gp.data.GemPyEngineConfig(
            backend=gp.data.AvailableBackends.numpy,
            dtype='float32'
        )
    )

    grav_forward = sol.gravity
    print(f"✓ Forward model computed")
    print(f"  Gravity range: [{grav_forward.min():.2f}, {grav_forward.max():.2f}] µGal")

# %%
# Visualize Forward Model Results
# --------------------------------

if DATA_AVAILABLE:
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

if DATA_AVAILABLE:
    # Convert observed from mGal to µGal
    observed_ugal = observed_gravity * 1000

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
    residuals = observed_ugal - grav_forward
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

if DATA_AVAILABLE:
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

if DATA_AVAILABLE:
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
#
# See the full script at: ``examples/geomodels/SimpleGravimetricResponse.py``

if not DATA_AVAILABLE:
    print("\n" + "="*60)
    print("To run this example with real data:")
    print("1. Install mineye: pip install -e .")
    print("2. Run Simple_Model_Tharsis.py to create the model")
    print("3. Ensure gravity data is in examples/Data/geophysical/")
    print("="*60)
