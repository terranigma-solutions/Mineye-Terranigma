from typing import Any

import arviz
import numpy as np
import xarray
from matplotlib import pyplot as plt

import gempy as gp
import gempy_viewer as gpv
from gempy_probability.modules.plot.plot_gempy import plot_gempy
from gempy_probability.modules.plot.plot_posterior import default_blue, default_red
from gempy_probability.modules.fields import fields


def generate_gravity_uncertainty_plots(gravity_samples_norm, observed_gravity_ugal, xy_ravel) -> tuple[str, Any]:
    # Import visualization functions
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_map_interpolated
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_profiles
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_with_uncertainty

    # Convert PyTorch tensor to numpy if needed
    if hasattr(gravity_samples_norm, 'numpy'):
        gravity_samples_norm = gravity_samples_norm.numpy()

    unit_label = 'Aligned Gravity (μGal)'

    print(f"\n{'=' * 60}")
    print(f"EXTRACTED NORMALIZED SAMPLES")
    print(f"{'=' * 60}")
    print(f"Number of samples: {gravity_samples_norm.shape[0]}")
    print(f"Number of devices: {gravity_samples_norm.shape[1]}")
    print(f"Normalization was applied DURING inference (not post-processing)")
    print(f"All samples use consistent normalization parameters from observed data")
    print(f"{'=' * 60}\n")

    # 1. Comprehensive uncertainty visualization (with normalized data)
    plot_gravity_with_uncertainty(
        gravity_samples=gravity_samples_norm,
        xy_coords=xy_ravel,
        observed_data=observed_gravity_ugal,
        confidence_level=0.95,
        title="Gravity Uncertainty Propagation from Dip Uncertainty"
    )

    # 2. Profile plots with confidence bands (with normalized data)
    plot_gravity_uncertainty_profiles(
        gravity_samples=gravity_samples_norm,
        xy_coords=xy_ravel,
        observed_data=(observed_gravity_ugal),
        n_profiles=4,
        confidence_level=0.95
    )

    # 3. Interpolated uncertainty maps (smoother visualization with normalized data)
    plot_gravity_uncertainty_map_interpolated(
        gravity_samples=gravity_samples_norm,
        xy_coords=xy_ravel,
        observed_data=(observed_gravity_ugal),
        grid_resolution=100
    )
    return gravity_samples_norm, unit_label


def gempy_viz(geo_model: gp.data.GeoModel, prior_inference_data: arviz.InferenceData,
              n_samples=20, ve=5):
    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.OCTREE],
        reset=True
    )

    geo_model.geophysics_input = None

    gp.compute_model(gempy_model=geo_model)

    p2d = gpv.plot_2d(
        model=geo_model,
        show_topography=False,
        legend=False,
        show_lith=False,
        show_data=False,
        show=False,
        ve=ve
    )
    plot_gempy(
        geo_model=geo_model,
        n_samples=n_samples,
        samples=(prior_inference_data.prior[r'dips'].values[0, :]),
        update_model_fn=_update_model_for_plotting,
        gempy_plot=p2d
    )

    if hasattr(prior_inference_data, 'posterior') and True:
        plot_gempy(
            geo_model=geo_model,
            n_samples=n_samples,
            samples=(prior_inference_data.posterior[r'dips'].values[0, :]),
            update_model_fn=_update_model_for_plotting,
            gempy_plot=p2d,
            contour_colors=[default_red],
        )

    return p2d


def compute_probability_density_fields(geo_model: gp.data.GeoModel, inference_data: xarray.Dataset,
                                       n_samples=50) -> fields.OnlineProbability:
    from gempy.core.data.grid_modules import RegularGrid
    geo_model.grid.dense_grid = RegularGrid(
        geo_model.grid.extent,
        np.array([50, 50, 50])
    )
    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.DENSE],
        reset=True
    )

    geo_model.geophysics_input = None

    gp.compute_model(gempy_model=geo_model)
    lith = geo_model.solutions.raw_arrays.lith_block

    online_prob = fields.OnlineProbability(
        n_cells=lith.shape[0],
        unique_lithologies=(np.unique(lith))
    )

    for i in np.linspace(0, inference_data.draw.size - 1, n_samples, dtype=int):
        _update_model_for_plotting(
            geo_model=geo_model,
            sample_value=(inference_data[r'dips'].values[0, :][i]),
            sample_idx=i
        )
        gp.compute_model(gempy_model=geo_model)
        online_prob.update(geo_model.solutions.raw_arrays.lith_block)

    return online_prob


def gempy_viz_pro(geo_model: gp.data.GeoModel, prior_inference_data: arviz.InferenceData):
    p2d = gempy_viz(geo_model, prior_inference_data, n_samples=10)

    if len(x_all):
        x_all = np.concatenate(x_all);
        z_all = np.concatenate(z_all)
        _draw_kde(ax, x_all, z_all, gridsize=400, bw=0.03, alpha=0.5, cmap="Blues", zorder=35, lognorm=True)

    p2d.axes[0].set_title("Uncertainty: KDE background + representative realizations")


def plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm):
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))

    # Calculate shared limits once
    vmin_shared = min(np.min(observed_norm), np.min(forward_norm))
    vmax_shared = max(np.max(observed_norm), np.max(forward_norm))
    lims = [vmin_shared, vmax_shared]

    # Sort observed values once
    sorted_indices = np.argsort(observed_norm)
    sorted_observed = observed_norm[sorted_indices]

    # Plot forward models
    for fw in many_forward_norm:
        sorted_fw = fw[sorted_indices]
        ax.scatter(sorted_observed, sorted_fw, alpha=0.7, s=40,
                   edgecolors='black', linewidth=0.5)
        ax.plot(sorted_observed, sorted_fw, alpha=0.3, linewidth=1)

    # Set up plot attributes
    unit_label = r'$\mu$Gal'
    ax.set_xlabel(f'Observed ({unit_label})')
    ax.set_ylabel(f'Forward Model ({unit_label})')
    ax.set_title('Observed vs Forward Model Correlation')

    # Add 1:1 line and set limits
    ax.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')
    ax.set(xlim=lims, ylim=lims, aspect='equal')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Add correlation coefficient
    correlation = np.corrcoef(observed_norm, forward_norm)[0, 1]
    ax.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8),
            fontsize=12)

    fig.show()


def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: float, sample_idx: int):
    # # Modify the surface point
    gp.modify_orientations(
        geo_model=geo_model,
        dip=sample_value,
    )

def plot_probability_heatmap(data, group='prior'):
    # Get the data array from ArviZ InferenceData
    # Standard Shape: (chain, draw, n_points, n_classes)
    # Example: (1, 50, 57, 3)
    probs_da = data[group]['probs_pred']

    # 1. Average over 'chain' and 'draw' dimensions to get mean per point
    # Result Shape: (n_points, n_classes)
    mean_probs = probs_da.mean(dim=["chain", "draw"]).values

    # 2. Transpose for the Heatmap
    # We want Y-axis = Classes, X-axis = Points
    # Result Shape: (n_classes, n_points)
    heatmap_data = mean_probs.T

    plot_heat_map(group, heatmap_data)


def plot_heat_map(group, heatmap_data):
    import seaborn as sns
    # 3. Plot
    plt.figure(figsize=(14, 4))
    sns.heatmap(heatmap_data, cmap="Blues", vmin=0, vmax=1,
                annot=False,  # Set True if you want numbers inside cells
                cbar_kws={'label': 'Probability'})

    plt.title(f"Average Class Probabilities per Point ({group.capitalize()})")
    plt.xlabel("Point Index (0 to 56)")
    plt.ylabel("Rock Unit (Class)")

    # Fix Y-axis labels to be integers (0, 1, 2)
    plt.yticks(ticks=np.arange(heatmap_data.shape[0]) + 0.5,
               labels=np.arange(heatmap_data.shape[0]),
               rotation=0)

    # plt.tight_layout()
    plt.show()
