import os
import time
from typing import Any

from pandas import DataFrame

import gempy as gp
import gempy_viewer as gpv

import numpy as np
import geopandas as gpd


def test_simple_model_gravity(simple_geo_model, geophysical_dir):
    """Test reading and computing a geological model."""

    # Model parameters
    # Use actual gravity measurement device locations
    gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
    observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

    xy_ravel = np.column_stack([
            np.array(gravity_data.geometry.x.values),
            np.array(gravity_data.geometry.y.values),
            np.full(len(gravity_data), 0)  # Set Z to surface elevation
    ])
    
    _gravity_precomputations(
        density_plutonites=2.9, # kg/m³
        density_sedimentary_host=2.3, # kg/m³
        xy_ravel=xy_ravel,
        simple_geo_model=simple_geo_model
    )

    simple_geo_model.interpolation_options.mesh_extraction = False
    start_time = time.time()
    sol: gp.data.Solutions = gp.compute_model(simple_geo_model)
    elapsed_time = time.time() - start_time

    print(f"\n⏱️  Model computation time: {elapsed_time:.2f} seconds")

    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None

    grav = sol.gravity

    gpv.plot_3d(simple_geo_model, ve=5, image=True)


    if PLOT:=True:
        # For actual measurement locations, show as scatter plot with color-coded gravity values
        _plot_fw_gravity(grav, gravity_data, xy_ravel)
        _plot_comparison(observed_gravity, grav, xy_ravel)


def _plot_comparison(observed_gravity, grav, xy_ravel, 
                     normalization_method='zscore'):
    import matplotlib.pyplot as plt
    print("\n=== Observed vs Predicted Comparison ===")

    # Convert units: observed is in mGal, predicted in μGal
    observed_ugal = observed_gravity * 1000  # Convert mGal to μGal
    forward_model = grav.numpy().copy()

    # Apply normalization if enabled
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


        print(f"  Observed stats (normalized): mean={np.mean(observed_norm):.3f}, std={np.std(observed_norm):.3f}")
        print(f"  Forward stats (normalized):  mean={np.mean(forward_norm):.3f}, std={np.std(forward_norm):.3f}")
    else:
        raise ValueError(f"Invalid normalization method: {normalization_method}")

    residuals_norm = observed_norm - forward_norm
    
    # Create comparison plot with 4 subplots (add correlation plot)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Observed gravity
    scatter1 = ax1.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=observed_norm,
                           s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax1.set_title(f'Observed Gravity{"" if not True else f" ({normalization_method})"}')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label(f'Observed ({unit_label})')

    # Plot 2: Forward model gravity
    scatter2 = ax2.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=forward_norm,
                           s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax2.set_title(f'Forward Model Gravity{"" if not True else f" ({normalization_method})"}')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label(f'Forward Model ({unit_label})')

    # Plot 3: Residuals (observed - predicted)
    scatter3 = ax3.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=residuals_norm,
                           s=30, cmap='RdBu_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax3.set_title(f'Residuals (Observed - Forward Model){"" if not True else f" ({normalization_method})"}')
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


def _plot_fw_gravity(grav, gravity_data: DataFrame, xy_ravel: np.ndarray[tuple[Any, ...], np.dtype]):
    import matplotlib.pyplot as plt
    scatter = plt.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=grav, s=30,
                          cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)

    # Add colorbar for scatter plot
    cbar = plt.colorbar(scatter)
    cbar.set_label(r'Forward Model Gravity ($\mu$gal)')

    print(f"Plotting {len(xy_ravel)} actual measurement locations")

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


def _gravity_precomputations(density_plutonites: float, density_sedimentary_host: float, 
                             xy_ravel: np.ndarray, simple_geo_model: gp.data.GeoModel):
    print("Using actual gravity measurement locations...")
    print(f"Using {len(xy_ravel)} actual measurement points")

    print("Computing forward gravity model...")

    # Step 1: Set centered grid
    print("Setting up centered grid...")
    gp.set_centered_grid(
        grid=simple_geo_model.grid,
        centers=xy_ravel,
        resolution=np.array([10, 10, 15]),
        radius=np.array([5000, 5000, 5000])
    )

    # Step 2: Calculate gravity gradient (tz component)
    print("Calculating gravity gradient...")
    gravity_gradient = gp.calculate_gravity_gradient(simple_geo_model.grid.centered_grid)

    # Step 3: Configure geophysics input
    print("Configuring geophysics input...")
    simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
        tz=gravity_gradient,
        densities=np.array([density_sedimentary_host, density_plutonites])  # kg/m³ for different formations,
    )
