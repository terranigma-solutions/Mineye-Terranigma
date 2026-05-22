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
    p2d.axes[0].set_title("Uncertainty: representative realizations")


def plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm, unit_label=r'$\mu$Gal'):
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


def plot_probability_heatmap(data, group='prior', slice_idx=None):
    # Get the data array from ArviZ InferenceData
    # Standard Shape: (chain, draw, n_points, n_classes)
    # Example: (1, 50, 57, 3)
    probs_da = data[group]['probs_pred']

    if slice_idx is not None:
        probs_da = probs_da.isel(Joint_Obs_1_dim_0=slice_idx)

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


def generate_paper_quality_figure(geo_model: gp.data.GeoModel, online_prob, topography_path, output_filename="probability_fields_paper.png"):
    import matplotlib.pyplot as plt
    import gempy_viewer as gpv
    from gempy_viewer.modules.plot_2d.visualization_2d import Plot2D
    import numpy as np
    import os

    # Set up publication-style parameters
    plt.rcParams.update({
        'font.size': 11,
        'axes.labelsize': 12,
        'axes.titlesize': 13,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'figure.titlesize': 16,
    })

    fig, axes = plt.subplots(2, 2, figsize=(16, 12), dpi=300, constrained_layout=True)
    ax_baseline = axes[0, 0]
    ax_prob0 = axes[0, 1]
    ax_prob1 = axes[1, 0]
    ax_entropy = axes[1, 1]

    # Helper function to plot onto a specific axis using gempy_viewer
    def plot_to_ax(ax, override_grid=None, show_lith=True, is_entropy=False, title="", cmap='viridis'):
        original_create_figure = Plot2D.create_figure
        try:
            def mock_create_figure(self, *args, **kwargs):
                self.fig = fig
                self.axes = [ax]
                return fig, self.axes
            Plot2D.create_figure = mock_create_figure
            
            kwargs_lith = {'cmap': cmap, 'norm': None} if override_grid is not None else None
            
            gpv.plot_2d(
                geo_model,
                override_regular_grid=override_grid,
                show_data=True,
                show_lith=show_lith,
                show_topography=True,
                ve=5,
                show=False,
                kwargs_lithology=kwargs_lith,
                legend=False
            )
            
            # Post-process axis to make it look clean
            ax.set_title(title, pad=10, fontweight='bold')
            ax.set_xlabel("X coordinate (m)")
            ax.set_ylabel("Z coordinate (m)")
            ax.grid(True, linestyle=":", alpha=0.5)
            
            # Add a colorbar for continuous fields
            if override_grid is not None:
                if len(ax.images) > 0:
                    mappable = ax.images[0]
                    cbar = fig.colorbar(mappable, ax=ax, fraction=0.046, pad=0.04)
                    cbar.set_label("Entropy (Uncertainty)" if is_entropy else "Probability", fontsize=11)
                    if not is_entropy:
                        cbar.mappable.set_clim(0, 1)
        finally:
            Plot2D.create_figure = original_create_figure

    # Panel A: Geological Baseline Section
    plot_to_ax(ax_baseline, title="A) Baseline Geological Cross-Section", show_lith=True)
    
    # Add a custom legend for geologic formations to Panel A to make it publication-ready
    try:
        from matplotlib.patches import Patch
        elements = geo_model.structural_frame.structural_elements
        legend_elements = []
        for element in elements:
            name = element.name
            color = element.color if hasattr(element, 'color') else 'gray'
            legend_elements.append(Patch(facecolor=color, edgecolor='black', label=name, alpha=0.8))
        ax_baseline.legend(handles=legend_elements, loc='lower left', frameon=True, facecolor='white', framealpha=0.9)
    except Exception as e:
        print(f"Could not add custom legend to baseline plot: {e}")

    # Panel B: Probability of Lithology 0 (Plutonites)
    plot_to_ax(ax_prob0, override_grid=online_prob.probability_field[0], title="B) Occurrence Probability: Plutonites", cmap='viridis')

    # Panel C: Probability of Lithology 1 (Sedimentary Host)
    plot_to_ax(ax_prob1, override_grid=online_prob.probability_field[1], title="C) Occurrence Probability: Sedimentary Host", cmap='viridis')

    # Panel D: Information Entropy (Uncertainty)
    plot_to_ax(ax_entropy, override_grid=online_prob.entropy, is_entropy=True, title="D) Information Entropy (Structural Uncertainty)", cmap='magma')

    # Save to file in project root
    project_root = "/home/leguark/PycharmProjects/Mineye-Terranigma"
    save_path = os.path.join(project_root, output_filename)
    fig.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"\n[SUCCESS] Generated paper-quality figure saved to: {save_path}")
    plt.close(fig)


def probability_fields_for(geo_model, inference_data, topography_path):
    online_prob = compute_probability_density_fields(
        geo_model,
        inference_data,
        n_samples=5
    )
    import gempy_viewer as gpv
    import gempy as gp
    import os

    # Set up topography early so that both 2D and 3D plots use it
    simple_geo_model = geo_model
    gp.set_topography_from_file(
        grid=simple_geo_model.grid,
        filepath=topography_path,
        crop_to_extent=[
                simple_geo_model.grid.extent[0],
                simple_geo_model.grid.extent[2],
                simple_geo_model.grid.extent[1],
                simple_geo_model.grid.extent[3]
        ]
    )

    gp.compute_model(geo_model)

    # 1. Generate the paper-quality, multi-panel summary figure
    generate_paper_quality_figure(geo_model, online_prob, topography_path)

    # 2. Original legacy plots (for compatibility/visual feedback)
    if True:
        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.probability_field[0],
            show_data=True,
            ve=5,
            kwargs_lithology={
                    'cmap': 'viridis',
                    'norm': None
            }
        )

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.probability_field[1],
            show_data=True,
            ve=5,
            kwargs_lithology={
                    'cmap': 'viridis',
                    'norm': None
            }
        )

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.entropy,
            show_data=True,
            ve=5,
            kwargs_lithology={
                    'cmap': 'magma',
                    'norm': None
            }
        )

    # * We inject the entropy field into the model
    geo_model.solutions.raw_arrays.scalar_field_matrix[0] = online_prob.entropy
    p3d = gpv.plot_3d(
        model=geo_model,
        active_scalar_field="sf_0",
        show_scalar=True,
        show_lith=False,
        show_topography=True,
        image=True,
        ve=4,
        threshold_kwargs={'value': [0.1, 0.9], 'invert': False},
        kwargs_pyvista_bounds={
                'show_xlabels': False,
                'show_ylabels': False,
                'show_zlabels': False,
        }
    )
    
    # This is a pyvista renderer
    p3d.p
    
    return online_prob
