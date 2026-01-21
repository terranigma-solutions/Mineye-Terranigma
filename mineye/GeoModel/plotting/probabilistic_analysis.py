from typing import Any

import numpy as np
import torch
from matplotlib import pyplot as plt
from pandas import DataFrame

from mineye.GeoModel.geophysics import compute_alignment_params, align_forward_to_observed


def plot_comparison(observed, fw_values, xy_ravel, unit_label,
                    normalization_method='align_to_reference'):
    import matplotlib.pyplot as plt
    print("\n=== Observed vs Predicted Comparison ===")

    # Convert units: observed is in mGal, predicted in μGal
    observed_ugal = observed * 1000  # Convert mGal to μGal
    forward_model = -fw_values

    params = compute_alignment_params(
        observed=observed_ugal,
        baseline_forward=forward_model,
        method=normalization_method
    )

    forward_norm = align_forward_to_observed(
        forward=forward_model,
        params=params,
    )

    if isinstance(forward_norm, torch.Tensor):
        forward_norm = forward_norm.numpy()

    plot_geophysics_comparison(forward_norm, normalization_method, observed_ugal, xy_ravel, unit_label=unit_label)


def plot_geophysics_comparison(forward_norm, normalization_method, observed_ugal,
                               xy_ravel, show=True, unit_label=r'$\mu$Gal'):
    observed_norm = observed_ugal

    residuals_norm = observed_norm - forward_norm

    # Compute SHARED colorbar limits for observed and forward plots
    vmin_shared = min(np.min(observed_norm), np.min(forward_norm))
    vmax_shared = max(np.max(observed_norm), np.max(forward_norm))

    # Create comparison plot with 4 subplots (add correlation plot)
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12), layout='constrained')

    # Plot 1: Observed gravity (with shared colorbar limits)
    scatter1 = ax1.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=observed_norm,
                           s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5,
                           vmin=vmin_shared, vmax=vmax_shared)
    ax1.set_title(f'Observed {"" if not True else f" ({normalization_method})"}')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label(f'Observed ({unit_label})')

    # Plot 2: Forward model gravity (with shared colorbar limits)
    scatter2 = ax2.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=forward_norm,
                           s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5,
                           vmin=vmin_shared, vmax=vmax_shared)
    ax2.set_title(f'Forward Model {"" if not True else f" ({normalization_method})"}')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label(f'Forward Model ({unit_label})')

    # Plot 3: Residuals (observed - predicted) - use symmetric colorbar
    residual_limit = max(abs(np.min(residuals_norm)), abs(np.max(residuals_norm)))
    scatter3 = ax3.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=residuals_norm,
                           s=30, cmap='RdBu_r', alpha=0.8, edgecolors='black', linewidth=0.5,
                           vmin=-residual_limit, vmax=residual_limit)
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

    # Add 1:1 line with shared limits
    lims = [vmin_shared, vmax_shared]
    ax4.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')
    ax4.set_xlim(lims)
    ax4.set_ylim(lims)
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    ax4.set_aspect('equal', adjustable='box')

    # Calculate and display correlation coefficient
    correlation = np.corrcoef(observed_norm, forward_norm)[0, 1]
    ax4.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=12)

    if show:
        plt.show()


def plot_fw_geophysics(fw_values, observed_data: DataFrame, xy_ravel: np.ndarray[tuple[Any, ...], np.dtype], label: str, title: str):
    import matplotlib.pyplot as plt
    scatter = plt.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=fw_values, s=30,
                          cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)

    # Add colorbar for scatter plot
    cbar = plt.colorbar(scatter)
    cbar.set_label(label)

    print(f"Plotting {len(xy_ravel)} actual measurement locations")

    # Always show actual measurement point locations regardless of grid type
    actual_measurement_points = np.column_stack([
            np.array(observed_data.geometry.x.values),
            np.array(observed_data.geometry.y.values)
    ])

    plt.scatter(actual_measurement_points[:, 0], actual_measurement_points[:, 1],
                s=15, c='red', marker='x', alpha=0.5, linewidth=1.5,
                label=f'Actual Measurement Points (n={len(actual_measurement_points)})')

    plt.legend(loc='upper right', framealpha=0.9)
    plt.title(title)
    plt.show()


def plot_gravity_with_uncertainty(gravity_samples: np.ndarray, xy_coords: np.ndarray,
                                  observed_data: np.ndarray = None,
                                  confidence_level: float = 0.95,
                                  title: str = "Gravity Prediction with Uncertainty"):
    """
    Plot spatial distribution of gravity predictions with uncertainty.

    Args:
        gravity_samples: Array of shape (n_samples, n_devices) - MCMC/prior samples
        xy_coords: Array of shape (n_devices, 2 or 3) - spatial coordinates of measurement points
        observed_data: Optional array of shape (n_devices,) - observed gravity values
        confidence_level: Confidence interval level (default 0.95 for 95% CI)
        title: Plot title
    """
    import matplotlib.pyplot as plt

    # Compute statistics
    mean_gravity = np.mean(gravity_samples, axis=0)
    std_gravity = np.std(gravity_samples, axis=0)

    # Compute confidence intervals
    lower_percentile = (1 - confidence_level) / 2 * 100
    upper_percentile = (1 + confidence_level) / 2 * 100
    lower_ci = np.percentile(gravity_samples, lower_percentile, axis=0)
    upper_ci = np.percentile(gravity_samples, upper_percentile, axis=0)

    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(18, 14), layout='constrained')

    # ============ Plot 1: Mean gravity with uncertainty visualization ============
    ax1 = axes[0, 0]

    # Adaptive error bar scaling based on spatial extent
    spatial_extent = max(
        xy_coords[:, 0].max() - xy_coords[:, 0].min(),
        xy_coords[:, 1].max() - xy_coords[:, 1].min()
    )
    mean_std = np.mean(std_gravity)

    # Calculate adaptive scale: make error bars ~2% of spatial extent on average
    target_size = 0.02 * spatial_extent  # Target error bar size
    error_scale = target_size / mean_std if mean_std > 0 else 1.0

    # Cap the error scale to avoid extreme values
    error_scale = np.clip(error_scale, 0.1, 1000)

    # Main scatter plot with mean gravity
    scatter1 = ax1.scatter(xy_coords[:, 0], xy_coords[:, 1],
                           c=mean_gravity, s=100, cmap='viridis_r',
                           edgecolors='black', linewidth=1, zorder=3)

    # Add adaptive error bars
    # Only show if they won't dominate the plot
    if error_scale * mean_std < 0.1 * spatial_extent:
        ax1.errorbar(xy_coords[:, 0], xy_coords[:, 1],
                     xerr=std_gravity * error_scale, yerr=std_gravity * error_scale,
                     fmt='none', ecolor='red', alpha=0.3, capsize=3, zorder=2,
                     label=f'Uncertainty (scaled {error_scale:.1f}×)')
        ax1.legend(fontsize=9, loc='upper right')
    else:
        # If error bars would be too large, show text warning instead
        ax1.text(0.02, 0.98,
                 f'⚠ Large uncertainty\n(mean σ = {mean_std:.1f} μGal)\nSee Plot 2 for details',
                 transform=ax1.transAxes, fontsize=10, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='yellow', alpha=0.7))

    ax1.set_title(f'Mean Gravity ± {confidence_level * 100:.0f}% CI', fontsize=14, fontweight='bold')
    ax1.set_xlabel('X (m)', fontsize=12)
    ax1.set_ylabel('Y (m)', fontsize=12)
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label(r'Mean Gravity ($\mu$Gal)', fontsize=11)
    ax1.grid(True, alpha=0.3)

    # ============ Plot 2: Uncertainty (standard deviation) ============
    ax2 = axes[0, 1]
    scatter2 = ax2.scatter(xy_coords[:, 0], xy_coords[:, 1],
                           c=std_gravity, s=100, cmap='YlOrRd',
                           edgecolors='black', linewidth=1)

    ax2.set_title('Prediction Uncertainty (Std. Dev.)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('X (m)', fontsize=12)
    ax2.set_ylabel('Y (m)', fontsize=12)
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label(r'Standard Deviation ($\mu$Gal)', fontsize=11)
    ax2.grid(True, alpha=0.3)

    # ============ Plot 3: Coefficient of Variation (relative uncertainty) ============
    ax3 = axes[1, 0]
    cv = (std_gravity / np.abs(mean_gravity)) * 100  # in percentage
    scatter3 = ax3.scatter(xy_coords[:, 0], xy_coords[:, 1],
                           c=cv, s=100, cmap='plasma',
                           edgecolors='black', linewidth=1)

    ax3.set_title('Coefficient of Variation (Relative Uncertainty)', fontsize=14, fontweight='bold')
    ax3.set_xlabel('X (m)', fontsize=12)
    ax3.set_ylabel('Y (m)', fontsize=12)
    cbar3 = plt.colorbar(scatter3, ax=ax3)
    cbar3.set_label('CV (%)', fontsize=11)
    ax3.grid(True, alpha=0.3)

    # ============ Plot 4: Observed vs Predicted (if observed data provided) ============
    ax4 = axes[1, 1]

    if observed_data is not None:
        # Scatter plot with error bars
        yerr_lower = np.abs(mean_gravity - lower_ci)
        yerr_upper = np.abs(upper_ci - mean_gravity)

        ax4.errorbar(observed_data, mean_gravity, yerr=[yerr_lower, yerr_upper],
                     fmt='o', alpha=0.6, ecolor='gray', capsize=4, markersize=6,
                     label='Predictions with CI')
        # 1:1 line
        lims = [min(np.min(observed_data), np.min(lower_ci)),
                max(np.max(observed_data), np.max(upper_ci))]
        ax4.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')

        # Calculate metrics
        residuals = observed_data - mean_gravity
        rmse = np.sqrt(np.mean(residuals ** 2))
        correlation = np.corrcoef(observed_data, mean_gravity)[0, 1]

        ax4.set_xlabel('Observed Gravity (μGal)', fontsize=12)
        ax4.set_ylabel('Predicted Gravity (μGal)', fontsize=12)
        ax4.set_title('Observed vs Predicted with Uncertainty', fontsize=14, fontweight='bold')
        ax4.legend(fontsize=10)
        ax4.grid(True, alpha=0.3)

        # Add metrics text box
        textstr = f'R = {correlation:.3f}\nRMSE = {rmse:.2f} μGal'
        ax4.text(0.05, 0.95, textstr, transform=ax4.transAxes, fontsize=11,
                 verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    else:
        # Distribution of predictions at each location (violin-style)
        n_devices = min(10, gravity_samples.shape[1])  # Limit to first 10 devices for clarity
        positions = np.arange(n_devices)

        violin_parts = ax4.violinplot([gravity_samples[:, i] for i in range(n_devices)],
                                      positions=positions, widths=0.7,
                                      showmeans=True, showmedians=True)

        ax4.set_xlabel('Device Index', fontsize=12)
        ax4.set_ylabel('Gravity (μGal)', fontsize=12)
        ax4.set_title(f'Distribution of Predictions (first {n_devices} devices)', fontsize=14, fontweight='bold')
        ax4.set_xticks(positions)
        ax4.grid(True, alpha=0.3, axis='y')

    plt.suptitle(title, fontsize=16, fontweight='bold', y=0.995)
    plt.show()

    # Print summary statistics
    print("\n" + "=" * 60)
    print("UNCERTAINTY SUMMARY STATISTICS")
    print("=" * 60)
    print(f"Mean gravity:        {np.mean(mean_gravity):.2f} ± {np.std(mean_gravity):.2f} μGal")
    print(f"Mean uncertainty:    {np.mean(std_gravity):.2f} μGal")
    print(f"Max uncertainty:     {np.max(std_gravity):.2f} μGal")
    print(f"Min uncertainty:     {np.min(std_gravity):.2f} μGal")
    print(f"Mean CV:             {np.mean(cv):.2f}%")
    if observed_data is not None:
        print(f"RMSE:                {rmse:.2f} μGal")
        print(f"Correlation (R):     {correlation:.3f}")
    print("=" * 60 + "\n")


def plot_gravity_uncertainty_profiles(gravity_samples: np.ndarray, xy_coords: np.ndarray,
                                      observed_data: np.ndarray = None,
                                      n_profiles: int = 4,
                                      confidence_level: float = 0.95):
    """
    Plot line profiles showing gravity with confidence bands.

    Args:
        gravity_samples: Array of shape (n_samples, n_devices)
        xy_coords: Array of shape (n_devices, 2 or 3)
        observed_data: Optional observed gravity values
        n_profiles: Number of profile directions to plot
        confidence_level: Confidence interval level
    """
    import matplotlib.pyplot as plt

    # Compute statistics
    mean_gravity = np.mean(gravity_samples, axis=0)
    lower_percentile = (1 - confidence_level) / 2 * 100
    upper_percentile = (1 + confidence_level) / 2 * 100
    lower_ci = np.percentile(gravity_samples, lower_percentile, axis=0)
    upper_ci = np.percentile(gravity_samples, upper_percentile, axis=0)

    fig, axes = plt.subplots(2, 2, figsize=(16, 10), layout='constrained')
    axes = axes.flatten()

    # Create profiles along different directions
    profile_directions = ['X', 'Y', 'Diagonal1', 'Diagonal2']

    for idx, (ax, direction) in enumerate(zip(axes, profile_directions[:n_profiles])):
        # Sort points based on direction
        if direction == 'X':
            sort_idx = np.argsort(xy_coords[:, 0])
            distance = xy_coords[sort_idx, 0]
            xlabel = 'X coordinate (m)'
        elif direction == 'Y':
            sort_idx = np.argsort(xy_coords[:, 1])
            distance = xy_coords[sort_idx, 1]
            xlabel = 'Y coordinate (m)'
        elif direction == 'Diagonal1':
            # Sort by x + y
            diagonal_coord = xy_coords[:, 0] + xy_coords[:, 1]
            sort_idx = np.argsort(diagonal_coord)
            distance = diagonal_coord[sort_idx]
            xlabel = 'X + Y coordinate (m)'
        else:  # Diagonal2
            # Sort by x - y
            diagonal_coord = xy_coords[:, 0] - xy_coords[:, 1]
            sort_idx = np.argsort(diagonal_coord)
            distance = diagonal_coord[sort_idx]
            xlabel = 'X - Y coordinate (m)'

        # Plot confidence band
        ax.fill_between(distance, lower_ci[sort_idx], upper_ci[sort_idx],
                        alpha=0.3, color='lightblue', label=f'{confidence_level * 100:.0f}% CI')

        # Plot mean
        ax.plot(distance, mean_gravity[sort_idx], 'b-', linewidth=2, label='Mean prediction')

        # Plot individual samples (subset for clarity)
        n_sample_lines = min(20, gravity_samples.shape[0])
        sample_indices = np.linspace(0, gravity_samples.shape[0] - 1, n_sample_lines, dtype=int)
        for sample_idx in sample_indices:
            ax.plot(distance, gravity_samples[sample_idx, sort_idx],
                    'gray', alpha=0.1, linewidth=0.5)

        # Plot observed data if available
        if observed_data is not None:
            ax.scatter(distance, observed_data[sort_idx], c='red', s=30,
                       zorder=5, label='Observed', edgecolors='black', linewidth=0.5)

        ax.set_xlabel(xlabel, fontsize=11)
        ax.set_ylabel('Gravity (μGal)', fontsize=11)
        ax.set_title(f'Profile along {direction}', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9, loc='best')
        ax.grid(True, alpha=0.3)

    plt.suptitle('Gravity Profiles with Uncertainty', fontsize=14, fontweight='bold')
    plt.show()


def plot_gravity_uncertainty_map_interpolated(gravity_samples: np.ndarray, xy_coords: np.ndarray,
                                              observed_data: np.ndarray = None,
                                              grid_resolution: int = 100):
    """
    Create interpolated uncertainty maps for smoother visualization.

    Args:
        gravity_samples: Array of shape (n_samples, n_devices)
        xy_coords: Array of shape (n_devices, 2 or 3)
        observed_data: Optional observed gravity values
        grid_resolution: Resolution of interpolation grid
    """
    import matplotlib.pyplot as plt
    from scipy.interpolate import griddata

    # Compute statistics
    mean_gravity = np.mean(gravity_samples, axis=0)
    std_gravity = np.std(gravity_samples, axis=0)

    # Create regular grid for interpolation
    xi = np.linspace(xy_coords[:, 0].min(), xy_coords[:, 0].max(), grid_resolution)
    yi = np.linspace(xy_coords[:, 1].min(), xy_coords[:, 1].max(), grid_resolution)
    xi_grid, yi_grid = np.meshgrid(xi, yi)

    # Interpolate mean and std
    mean_interp = griddata(xy_coords[:, :2], mean_gravity, (xi_grid, yi_grid), method='cubic')
    std_interp = griddata(xy_coords[:, :2], std_gravity, (xi_grid, yi_grid), method='cubic')

    # Create figure
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), layout='constrained')

    # Plot 1: Mean gravity (interpolated)
    ax1 = axes[0]
    im1 = ax1.contourf(xi_grid, yi_grid, mean_interp, levels=20, cmap='viridis_r')
    ax1.scatter(xy_coords[:, 0], xy_coords[:, 1], c='red', s=20,
                marker='x', linewidth=1.5, label='Measurement locations')
    if observed_data is not None:
        # Add contour lines for observed data
        obs_interp = griddata(xy_coords[:, :2], observed_data, (xi_grid, yi_grid), method='cubic')
        ax1.contour(xi_grid, yi_grid, obs_interp, levels=10, colors='white',
                    linewidths=0.5, alpha=0.5, linestyles='dashed')
    ax1.set_title('Mean Gravity (Interpolated)', fontsize=13, fontweight='bold')
    ax1.set_xlabel('X (m)', fontsize=11)
    ax1.set_ylabel('Y (m)', fontsize=11)
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label(r'Gravity ($\mu$Gal)', fontsize=10)
    ax1.legend(fontsize=9)

    # Plot 2: Uncertainty map (interpolated)
    ax2 = axes[1]
    im2 = ax2.contourf(xi_grid, yi_grid, std_interp, levels=20, cmap='YlOrRd')
    ax2.scatter(xy_coords[:, 0], xy_coords[:, 1], c='blue', s=20,
                marker='x', linewidth=1.5, label='Measurement locations')
    ax2.set_title('Prediction Uncertainty (Std. Dev.)', fontsize=13, fontweight='bold')
    ax2.set_xlabel('X (m)', fontsize=11)
    ax2.set_ylabel('Y (m)', fontsize=11)
    cbar2 = plt.colorbar(im2, ax=ax2)
    cbar2.set_label(r'Std. Dev. ($\mu$Gal)', fontsize=10)
    ax2.legend(fontsize=9)

    plt.show()
