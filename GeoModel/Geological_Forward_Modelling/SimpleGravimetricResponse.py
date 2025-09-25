import numpy as np
import gempy as gp
import gempy_viewer as gpv
import geopandas as gpd
import matplotlib.pyplot as plt
import os
import pickle

# ========== CONFIG ==========
BASE_DIR = os.getcwd()
data_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Input_Data'))
geomodel_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Output_Data'))
geophysical_dir = os.path.join(data_dir, 'Geophysical_Cleaned_Data')
forward_model_folder = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'GeoModel', 'Geological_Forward_Modelling'))

temp_inputs_dir = os.path.join(geomodel_dir, 'Simple-Models', 'temp_inputs')
pickle_model_path = os.path.join(temp_inputs_dir, 'simple_geo_model.pkl')
gravity_data_path = os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson')
contact_points_path = os.path.join(geomodel_dir, 'Simple-Models', 'contact_points.csv')

# Model parameters
density_plutonites = 2.9  # kg/m³
density_sedimentary_host = 2.3  # kg/m³
gravity_resolution = 20  # Number of gravity measurement points per axis

# Measurement grid options
use_actual_measurement_locations = False  # True: use gravity device locations, False: use regular grid

# Normalization options
normalize_data = True  # Enable/disable normalization
normalization_method = 'minmax'  # Options: 'zscore', 'minmax', 'mean_center', 'relative'

# Triggers
trigger_load_gravity_data = True
trigger_forward_modeling = True
trigger_comparison_plots = True
trigger_save_results = True

# ========== LOAD GEOLOGICAL MODEL AND GRAVITY DATA ==========
with open(pickle_model_path, 'rb') as f:
    geo_model = pickle.load(f)

# Load gravity data
gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
observed_gravity = gravity_data['VALU_BOU267'].values # in mGal

# p = gpv.plot_2d(
#     geo_model,
#     section_names=['topography'],   # this triggers the top-down geological map
#     show_topography=True,
#     show_lith=True,
#     show_boundaries=True,
#     show_data=True,
#     legend=False
# )

# ========== CREATE MEASUREMENT GRID ==========
extent = geo_model.grid.regular_grid.extent
min_x, max_x, min_y, max_y, min_z, max_z = extent

if use_actual_measurement_locations:
    # Use actual gravity measurement device locations
    print("Using actual gravity measurement locations...")
    xy_ravel = np.column_stack([
        np.array(gravity_data.geometry.x.values),
        np.array(gravity_data.geometry.y.values),
        np.full(len(gravity_data), max_z)  # Set Z to surface elevation
    ])
    print(f"Using {len(xy_ravel)} actual measurement points")
else:
    # Create regular measurement grid following official example pattern
    print("Using regular grid for measurements...")
    grav_res = gravity_resolution
    X = np.linspace(min_x, max_x, gravity_resolution)
    Y = np.linspace(min_y, max_y, gravity_resolution)
    Z = max_z  # Use topographic surface elevation as measurement height
    xyz = np.meshgrid(X, Y, Z)
    xy_ravel = np.vstack(list(map(np.ravel, xyz))).T
    print(f"Created regular grid: {gravity_resolution}x{gravity_resolution} = {len(xy_ravel)} points")

# ========== FORWARD MODELING ==========
if trigger_forward_modeling:
    print("Computing forward gravity model...")

    # Step 1: Set centered grid
    print("Setting up centered grid...")
    gp.set_centered_grid(
        grid=geo_model.grid,
        centers=xy_ravel,
        resolution=np.array([10, 10, 15]),
        radius=np.array([5000, 5000, 5000])
    )

    # Step 2: Calculate gravity gradient (tz component)
    print("Calculating gravity gradient...")
    gravity_gradient = gp.calculate_gravity_gradient(geo_model.grid.centered_grid)

    # Step 3: Configure geophysics input
    print("Configuring geophysics input...")
    geo_model.geophysics_input = gp.data.GeophysicsInput(
        tz=gravity_gradient,
        densities=np.array([density_sedimentary_host, density_plutonites])  # kg/m³ for different formations,
    )

    # Step 4: Compute forward model
    print("Computing forward model...")
    geo_model.interpolation_options.mesh_extraction = False
    sol = gp.compute_model(
        gempy_model=geo_model,
        engine_config=gp.data.GemPyEngineConfig(
            backend=gp.data.AvailableBackends.numpy,
            dtype='float32'
        )
    )

    grav = sol.gravity

    print("Forward modeling complete!")

print("\n=== Forward Modeling Complete ===")

# ========== PLOTTING ==========
if trigger_comparison_plots:
    if use_actual_measurement_locations:
        # For actual measurement locations, show as scatter plot with color-coded gravity values
        scatter = plt.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=grav, s=30,
                            cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)

        # Add colorbar for scatter plot
        cbar = plt.colorbar(scatter)
        cbar.set_label(r'Forward Model Gravity ($\mu$gal)')

        print(f"Plotting {len(xy_ravel)} actual measurement locations")

    else:
        # For regular grid, use imshow overlay as before
        grav_res = gravity_resolution

        # Add measurement grid points as small dots
        plt.scatter(xy_ravel[:, 0], xy_ravel[:, 1], s=1, c='white', alpha=0.5)

        # Overlay gravity field as image
        plt.imshow(grav.reshape(grav_res, grav_res),
                   extent=(xy_ravel[:, 0].min() + (xy_ravel[0, 0] - xy_ravel[1, 0]) / 2,
                           xy_ravel[:, 0].max() - (xy_ravel[0, 0] - xy_ravel[1, 0]) / 2,
                           xy_ravel[:, 1].min() + (xy_ravel[0, 1] - xy_ravel[grav_res, 1]) / 2,
                           xy_ravel[:, 1].max() - (xy_ravel[0, 1] - xy_ravel[grav_res, 1]) / 2),
                   cmap='viridis_r', origin='lower', alpha=.8)

        # Add colorbar for image
        cbar = plt.colorbar()
        cbar.set_label(r'Forward Model Gravity ($\mu$gal)')

        print(f"Plotting regular grid: {grav_res}x{grav_res} points")

    # Always show actual measurement point locations regardless of grid type
    actual_measurement_points = np.column_stack([
        np.array(gravity_data.geometry.x.values),
        np.array(gravity_data.geometry.y.values)
    ])

    plt.scatter(actual_measurement_points[:, 0], actual_measurement_points[:, 1],
               s=15, c='red', marker='x', alpha=0.8, linewidth=1.5,
               label=f'Actual Measurement Points (n={len(actual_measurement_points)})')

    plt.legend(loc='upper right', framealpha=0.9)
    plt.title('Forward Gravity Model Results')
    plt.show()

# ========== COMPARISON WITH OBSERVED DATA ==========
if trigger_comparison_plots and use_actual_measurement_locations:
    # Only do comparison if using actual measurement locations
    print("\n=== Observed vs Predicted Comparison ===")

    # Convert units: observed is in mGal, predicted in μGal
    observed_ugal = observed_gravity * 1000  # Convert mGal to μGal
    forward_model = grav.copy()

    # Apply normalization if enabled
    if normalize_data:
        print(f"Applying {normalization_method} normalization...")

        if normalization_method == 'zscore':
            # Z-score normalization (mean=0, std=1)
            obs_mean, obs_std = np.mean(observed_ugal), np.std(observed_ugal)
            fwd_mean, fwd_std = np.mean(forward_model), np.std(forward_model)

            observed_norm = (observed_ugal - obs_mean) / obs_std
            forward_norm = (forward_model - fwd_mean) / fwd_std

            unit_label = 'Z-score'

        elif normalization_method == 'minmax':
            # Min-max normalization (0 to 1)
            obs_min, obs_max = np.min(observed_ugal), np.max(observed_ugal)
            fwd_min, fwd_max = np.min(forward_model), np.max(forward_model)

            observed_norm = (observed_ugal - obs_min) / (obs_max - obs_min)
            forward_norm = (forward_model - fwd_min) / (fwd_max - fwd_min)

            unit_label = 'Normalized [0-1]'

        elif normalization_method == 'mean_center':
            # Mean centering (subtract mean)
            observed_norm = observed_ugal - np.mean(observed_ugal)
            forward_norm = forward_model - np.mean(forward_model)

            unit_label = 'Mean-centered (μGal)'

        elif normalization_method == 'relative':
            # Relative to range
            obs_range = np.max(observed_ugal) - np.min(observed_ugal)
            fwd_range = np.max(forward_model) - np.min(forward_model)

            observed_norm = observed_ugal / obs_range
            forward_norm = forward_model / fwd_range

            unit_label = 'Relative to range'

        residuals_norm = observed_norm - forward_norm

        print(f"  Observed stats (normalized): mean={np.mean(observed_norm):.3f}, std={np.std(observed_norm):.3f}")
        print(f"  Forward stats (normalized):  mean={np.mean(forward_norm):.3f}, std={np.std(forward_norm):.3f}")

    else:
        # No normalization - use original units with conversion
        observed_norm = observed_ugal
        forward_norm = forward_model
        residuals_norm = observed_ugal - forward_model
        unit_label = 'μGal'

    # Create comparison plot with 4 subplots (add correlation plot)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Observed gravity
    scatter1 = ax1.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=observed_norm,
                          s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax1.set_title(f'Observed Gravity{"" if not normalize_data else f" ({normalization_method})"}')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label(f'Observed ({unit_label})')

    # Plot 2: Forward model gravity
    scatter2 = ax2.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=forward_norm,
                          s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax2.set_title(f'Forward Model Gravity{"" if not normalize_data else f" ({normalization_method})"}')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label(f'Forward Model ({unit_label})')

    # Plot 3: Residuals (observed - predicted)
    scatter3 = ax3.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=residuals_norm,
                          s=30, cmap='RdBu_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax3.set_title(f'Residuals (Observed - Forward Model){"" if not normalize_data else f" ({normalization_method})"}')
    ax3.set_xlabel('X (m)')
    ax3.set_ylabel('Y (m)')
    cbar3 = plt.colorbar(scatter3, ax=ax3)
    cbar3.set_label(f'Residuals ({unit_label})')

    # Plot 4: Correlation plot (observed vs forward model)
    ax4.scatter(observed_norm, forward_norm, alpha=0.7, s=40, edgecolors='black', linewidth=0.5)
    ax4.set_xlabel(f'Observed ({unit_label})')
    ax4.set_ylabel(f'Forward Model ({unit_label})')
    ax4.set_title('Observed vs Forward Model Correlation')

    # Add 1:1 line
    lims = [min(ax4.get_xlim()[0], ax4.get_ylim()[0]), max(ax4.get_xlim()[1], ax4.get_ylim()[1])]
    ax4.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # Calculate and display correlation coefficient
    correlation = np.corrcoef(observed_norm, forward_norm)[0, 1]
    ax4.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=12)

    plt.tight_layout()
    plt.show()

    # Print enhanced statistics
    rmse = np.sqrt(np.mean(residuals_norm**2))
    mean_residual = np.mean(residuals_norm)
    std_residual = np.std(residuals_norm)
    mae = np.mean(np.abs(residuals_norm))  # Mean Absolute Error

    print(f"Comparison Statistics ({unit_label}):")
    print(f"  RMSE: {rmse:.4f}")
    print(f"  MAE:  {mae:.4f}")
    print(f"  Mean Residual: {mean_residual:.4f}")
    print(f"  Std Residual:  {std_residual:.4f}")
    print(f"  Correlation (R): {correlation:.4f}")
    print(f"  R²: {correlation**2:.4f}")
    print(f"  Number of points: {len(residuals_norm)}")

    if normalize_data:
        print(f"\nNormalization method: {normalization_method}")
        print(f"  Original observed range: {np.min(observed_ugal):.2f} to {np.max(observed_ugal):.2f} μGal")
        print(f"  Original forward range:  {np.min(forward_model):.2f} to {np.max(forward_model):.2f} μGal")
        print(f"  Normalized observed range: {np.min(observed_norm):.4f} to {np.max(observed_norm):.4f}")
        print(f"  Normalized forward range:  {np.min(forward_norm):.4f} to {np.max(forward_norm):.4f}")
