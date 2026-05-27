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


def generate_gravity_uncertainty_plots(gravity_samples_norm, observed_gravity_ugal, xy_ravel,
                                       unit_label: str = 'Aligned Gravity (μGal)',
                                       response_label: str = 'Gravity',
                                       title_prefix: str = 'Gravity') -> tuple[str, Any]:
    # Import visualization functions
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_map_interpolated
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_profiles
    from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_with_uncertainty

    # Convert PyTorch tensor to numpy if needed
    if hasattr(gravity_samples_norm, 'numpy'):
        gravity_samples_norm = gravity_samples_norm.numpy()

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
        title=f"{title_prefix} Uncertainty Propagation from Dip Uncertainty",
        unit_label=unit_label,
        response_label=response_label,
    )

    # 2. Profile plots with confidence bands (with normalized data)
    plot_gravity_uncertainty_profiles(
        gravity_samples=gravity_samples_norm,
        xy_coords=xy_ravel,
        observed_data=(observed_gravity_ugal),
        n_profiles=4,
        confidence_level=0.95,
        unit_label=unit_label,
        response_label=response_label,
    )

    # 3. Interpolated uncertainty maps (smoother visualization with normalized data)
    plot_gravity_uncertainty_map_interpolated(
        gravity_samples=gravity_samples_norm,
        xy_coords=xy_ravel,
        observed_data=(observed_gravity_ugal),
        grid_resolution=100,
        unit_label=unit_label,
        response_label=response_label,
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
    import os
    from itertools import cycle

    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.patches import Patch

    dense_grid = geo_model.grid.dense_grid
    resolution = np.asarray(dense_grid.resolution, dtype=int)
    extent = np.asarray(dense_grid.extent, dtype=float)
    section_y_idx = int(resolution[1] // 2)
    section_y = dense_grid.y_coord[section_y_idx]
    section_extent = [extent[0], extent[1], extent[4], extent[5]]

    def reshape_field(field: np.ndarray) -> np.ndarray:
        field = np.asarray(field)
        expected_size = int(np.prod(resolution))
        if field.size != expected_size:
            raise ValueError(
                f"Cannot reshape field with {field.size} cells to dense grid resolution "
                f"{tuple(resolution)} ({expected_size} cells)."
            )
        return field.reshape(tuple(resolution), order="C")

    def mid_y_section(field: np.ndarray) -> np.ndarray:
        return reshape_field(field)[:, section_y_idx, :].T

    def topography_profile() -> tuple[np.ndarray, np.ndarray] | None:
        topography = getattr(geo_model.grid, "topography", None)
        if topography is None or topography.values_2d.size == 0:
            return None

        topography_values = np.asarray(topography.values_2d)
        topography_y = topography_values[0, :, 1]
        topography_y_idx = int(np.argmin(np.abs(topography_y - section_y)))
        return topography_values[:, topography_y_idx, 0], topography_values[:, topography_y_idx, 2]

    elements = list(getattr(geo_model.structural_frame, "structural_elements", []))
    default_colors = cycle(plt.get_cmap("tab20").colors)

    def element_color(element):
        color = getattr(element, "color", None)
        return color if color is not None else next(default_colors)

    element_styles = [
            (
                    getattr(element, "name", f"Lithology {idx + 1}"),
                    element_color(element)
            )
            for idx, element in enumerate(elements)
    ]

    topo_profile = topography_profile()

    lithology_ids = np.asarray(online_prob.unique_lithologies, dtype=int)
    lithology_names = []
    lithology_colors = []
    for idx, lithology_id in enumerate(lithology_ids):
        if idx < len(element_styles):
            lithology_name, lithology_color = element_styles[idx]
        else:
            lithology_name, lithology_color = f"Lithology {lithology_id}", next(default_colors)
        lithology_names.append(lithology_name)
        lithology_colors.append(lithology_color)

    baseline = np.round(np.asarray(geo_model.solutions.raw_arrays.lith_block)).astype(int)
    baseline_section = mid_y_section(baseline)
    baseline_index = np.full_like(baseline_section, fill_value=np.nan, dtype=float)
    for idx, lithology_id in enumerate(lithology_ids):
        baseline_index[baseline_section == lithology_id] = idx

    with plt.rc_context({
            "font.size"       : 10,
            "axes.labelsize"  : 10,
            "axes.titlesize"  : 12,
            "xtick.labelsize" : 9,
            "ytick.labelsize" : 9,
            "figure.titlesize": 15,
            "savefig.dpi"     : 300,
    }):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
        axes = axes.ravel()
        fig.suptitle(
            f"Probabilistic geological inversion summary — mid-Y section ({section_y:.0f} m)",
            fontweight="bold"
        )

        lithology_cmap = ListedColormap(lithology_colors)
        lithology_norm = BoundaryNorm(np.arange(len(lithology_ids) + 1) - 0.5, len(lithology_ids))
        axes[0].imshow(
            baseline_index,
            origin="lower",
            extent=section_extent,
            aspect="auto",
            cmap=lithology_cmap,
            norm=lithology_norm,
            interpolation="nearest"
        )
        axes[0].set_title("A) Baseline geological model", loc="left", fontweight="bold")
        legend_handles = [
                Patch(facecolor=color, edgecolor="black", label=name)
                for name, color in zip(lithology_names, lithology_colors)
        ]
        axes[0].legend(
            handles=legend_handles,
            loc="lower left",
            frameon=True,
            facecolor="white",
            framealpha=0.9,
            fontsize=8
        )

        probability_axes = axes[1:3]
        probability_titles = ["B", "C"]
        for plot_idx, ax in enumerate(probability_axes):
            if plot_idx >= online_prob.probability_field.shape[0]:
                ax.axis("off")
                continue
            probability_section = mid_y_section(online_prob.probability_field[plot_idx])
            image = ax.imshow(
                probability_section,
                origin="lower",
                extent=section_extent,
                aspect="auto",
                cmap="viridis",
                vmin=0,
                vmax=1,
                interpolation="bilinear"
            )
            lithology_label = lithology_names[plot_idx] if plot_idx < len(lithology_names) else f"Lithology {plot_idx}"
            ax.contour(
                probability_section,
                levels=[0.5],
                colors="white",
                linewidths=1.0,
                origin="lower",
                extent=section_extent
            )
            ax.set_title(f"{probability_titles[plot_idx]}) Occurrence probability: {lithology_label}", loc="left", fontweight="bold")
            colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.02)
            colorbar.set_label("Probability")

        entropy_section = mid_y_section(online_prob.entropy)
        entropy_image = axes[3].imshow(
            entropy_section,
            origin="lower",
            extent=section_extent,
            aspect="auto",
            cmap="magma",
            vmin=0,
            interpolation="bilinear"
        )
        axes[3].set_title("D) Information entropy (structural uncertainty)", loc="left", fontweight="bold")
        entropy_colorbar = fig.colorbar(entropy_image, ax=axes[3], fraction=0.046, pad=0.02)
        entropy_colorbar.set_label("Entropy (bits)")

        for ax in axes:
            if not ax.axison:
                continue
            if topo_profile is not None:
                ax.plot(
                    topo_profile[0],
                    topo_profile[1],
                    color="black",
                    linewidth=1.0,
                    label="Topography",
                    zorder=5
                )
            ax.set_xlabel("X coordinate (m)")
            ax.set_ylabel("Elevation Z (m)")
            ax.grid(color="white", linestyle=":", linewidth=0.4, alpha=0.45)
            ax.set_xlim(extent[0], extent[1])
            ax.set_ylim(extent[4], extent[5])

        if not os.path.isabs(output_filename):
            output_filename = os.path.join(os.getcwd(), output_filename)
        fig.savefig(output_filename, bbox_inches="tight")
        print(f"\n[SUCCESS] Generated paper-quality figure saved to: {output_filename}")
        plt.close(fig)


def generate_pyvista_paper_quality_figure(geo_model: gp.data.GeoModel, online_prob,
                                          output_filename="probability_fields_paper_pyvista.png"):
    import os
    from itertools import cycle

    from matplotlib.colors import to_hex
    from matplotlib.patches import Patch

    dense_grid = geo_model.grid.dense_grid
    resolution = np.asarray(dense_grid.resolution, dtype=int)
    extent = np.asarray(dense_grid.extent, dtype=float)
    expected_size = int(np.prod(resolution))
    crop_y = extent[2] + 0.58 * (extent[3] - extent[2])

    def normalized_field(field: np.ndarray) -> np.ndarray:
        field = np.asarray(field, dtype=float).reshape(-1)
        if field.size != expected_size:
            raise ValueError(
                f"Cannot render field with {field.size} cells on dense grid "
                f"{tuple(resolution)} ({expected_size} cells)."
            )
        return field

    elements = list(getattr(geo_model.structural_frame, "structural_elements", []))
    default_colors = cycle(plt.get_cmap("tab20").colors)
    lithology_ids = np.asarray(online_prob.unique_lithologies, dtype=int)
    lithology_names = []
    lithology_colors = []
    for idx, lithology_id in enumerate(lithology_ids):
        element = elements[idx] if idx < len(elements) else None
        lithology_names.append(getattr(element, "name", f"Lithology {lithology_id}"))
        color = getattr(element, "color", None) if element is not None else None
        lithology_colors.append(to_hex(color if color is not None else next(default_colors)))

    def positive_field(field: np.ndarray) -> np.ndarray:
        field = normalized_field(field).copy()
        field[field <= 0] = np.nan
        return field

    def configure_camera():
        # Clipped extent: only show up to crop_y in Y direction
        x_mid = float((extent[0] + extent[1]) / 2)
        y_mid = float((extent[2] + crop_y) / 2)
        z_mid = float((extent[4] + extent[5]) / 2)
        
        dx = float(extent[1] - extent[0])
        dz = float(extent[5] - extent[4])
        
        # Use the largest horizontal/vertical span to scale the camera distance
        span = max(dx, dz) * 0.001
        
        camera_position = [
            # 1. Position:
            # Offset X by -0.7 to expose the side (isometric feel)
            # Pull Y back by -1.5 to maintain the primary view toward positive Y
            # Lift Z up by +0.8 to look down from above
            (
                x_mid - 0.7 * span,  
                y_mid - 1.5 * span,  
                z_mid + 0.8 * span   
            ),
            
            # 2. Focal Point: Centered on the model
            (x_mid, y_mid, z_mid),
            
            # 3. View Up: Z-axis is up
            (0, 0, 1),
        ]
        
        return {
            "position": camera_position,
        }
    def apply_scalar_colormap(plotter, cmap, clim):
        if cmap is None:
            return
        try:
            import pyvista as pv

            lookup_table = pv.LookupTable(cmap=cmap, scalar_range=clim)
            lookup_table.nan_opacity = 0.0
            plotter.regular_grid_actor.mapper.lookup_table = lookup_table
        except Exception as exc:
            print(f"[WARNING] Could not apply PyVista colormap '{cmap}': {exc}")

    def render_panel(field, scalar_name, clim, opacity=0.95, mask_zero=False, cmap=None, show_lith=False):
        scalar_bar_args = {
                "title": scalar_name.replace("_", " ").title(),
                "title_font_size": 18,
                "label_font_size": 14,
                "font_family": "arial",
                "vertical": False,
                "height": 0.055,
                "width": 0.52,
                "position_x": 0.24,
                "position_y": 0.055,
        }
        camera_args = configure_camera()

        if show_lith:
            p3d = gpv.plot_3d(
                model=geo_model,
                show_scalar=False,
                show_lith=True,
                show_topography=True,
                image=True,
                show=False,
                ve=4,
                kwargs_plot_structured_grid={
                        "opacity": opacity,
                },
                kwargs_pyvista_bounds={
                        'show_xlabels': False,
                        'show_ylabels': False,
                        'show_zlabels': False,
                },
                kwargs_pyvista_camera= camera_args,
            )
        else:
            scalar_field = positive_field(field) if mask_zero else normalized_field(field)
            original_scalar_field = np.asarray(geo_model.solutions.raw_arrays.scalar_field_matrix[0]).copy()
            geo_model.solutions.raw_arrays.scalar_field_matrix[0] = scalar_field
            from gempy_viewer.modules.plot_3d import drawer_structured_grid_3d

            original_plot_structured_grid = drawer_structured_grid_3d.plot_structured_grid

            def plot_structured_grid_with_cmap(*args, **kwargs):
                if cmap is not None:
                    kwargs["cmap"] = cmap
                return original_plot_structured_grid(*args, **kwargs)

            drawer_structured_grid_3d.plot_structured_grid = plot_structured_grid_with_cmap
            try:
                p3d = gpv.plot_3d(
                    model=geo_model,
                    active_scalar_field="sf_0",
                    show_scalar=True,
                    show_lith=False,
                    show_topography=True,
                    image=True,
                    show=False,
                    ve=4,
                    threshold_kwargs={'value': [0.05, 10], 'invert': False},
                    kwargs_plot_topography={
                            "opacity": 0.,
                    },
                    kwargs_plot_structured_grid={
                            "clim": clim,
                            "opacity": opacity,
                            "nan_opacity": 0.0,
                            "scalar_bar_args": scalar_bar_args,
                    },
                    kwargs_pyvista_bounds={
                            'show_xlabels': False,
                            'show_ylabels': False,
                            'show_zlabels': False,
                    },
                    kwargs_pyvista_camera= camera_args,
                )
            finally:
                drawer_structured_grid_3d.plot_structured_grid = original_plot_structured_grid
                geo_model.solutions.raw_arrays.scalar_field_matrix[0] = original_scalar_field

        plotter = p3d.p
        apply_scalar_colormap(p3d, cmap, clim)

        # Clip the structured grid vertically parallel to Y (across the C-direction)
        if p3d.regular_grid_mesh is not None:
            clip_bounds = [float(extent[0]), float(extent[1]),
                           float(extent[2]), float(crop_y),
                           float(extent[4]), float(extent[5])]
            try:
                clipped_mesh = p3d.regular_grid_mesh.clip_box(bounds=clip_bounds, invert=False)
                if show_lith:
                    new_actor = plotter.add_mesh(
                        clipped_mesh, show_scalar_bar=False,
                        interpolate_before_map=True, opacity=opacity
                    )
                else:
                    new_actor = plotter.add_mesh(
                        clipped_mesh, cmap=cmap, show_scalar_bar=True,
                        interpolate_before_map=True, opacity=opacity,
                        clim=clim, nan_opacity=0.0,
                        scalar_bar_args=scalar_bar_args,
                    )
                plotter.remove_actor(p3d.regular_grid_actor)
                p3d.regular_grid_actor = new_actor
                p3d.regular_grid_mesh = clipped_mesh
            except Exception as clip_exc:
                print(f"[WARNING] Could not clip grid: {clip_exc}")

        image = plotter.screenshot(return_img=True, transparent_background=False)
        plotter.close()
        return image

    second_probability = online_prob.probability_field[1] if online_prob.probability_field.shape[0] > 1 else online_prob.probability_field[0]
    entropy_max = float(np.nanmax(online_prob.entropy)) if np.isfinite(online_prob.entropy).any() else 1.0
    panel_specs = [
            (
                    "A) Baseline geological model",
                    None,
                    "lithology",
                    (-0.5, len(lithology_ids) - 0.5),
                    0.72,
                    False,
                    None,
                    True,
            ),
            (
                    f"B) Occurrence probability: {lithology_names[0] if lithology_names else 'Lithology 0'}",
                    online_prob.probability_field[0],
                    "probability",
                    (0, 1),
                    0.86,
                    True,
                    "viridis",
                    False,
            ),
            (
                    f"C) Occurrence probability: {lithology_names[1] if len(lithology_names) > 1 else 'Lithology 1'}",
                    second_probability,
                    "probability",
                    (0, 1),
                    0.86,
                    True,
                    "viridis",
                    False,
            ),
            (
                    "D) Information entropy (structural uncertainty)",
                    online_prob.entropy,
                    "entropy",
                    (0, entropy_max),
                    0.92,
                    True,
                    "magma",
                    False,
            ),
    ]

    try:
        rendered_panels = [render_panel(*spec[1:]) for spec in panel_specs]
    except Exception as exc:
        print(f"[WARNING] Could not generate PyVista paper figure: {exc}")
        return None

    with plt.rc_context({
            "font.size"       : 10,
            "axes.titlesize"  : 12,
            "figure.titlesize": 15,
            "savefig.dpi"     : 300,
    }):
        fig, axes = plt.subplots(2, 2, figsize=(14, 10), constrained_layout=True)
        axes = axes.ravel()
        fig.suptitle("Probabilistic geological inversion summary — cropped 3D PyVista view", fontweight="bold")
        for ax, (title, *_), image in zip(axes, panel_specs, rendered_panels):
            ax.imshow(image)
            ax.set_title(title, loc="left", fontweight="bold")
            ax.axis("off")

        if not os.path.isabs(output_filename):
            output_filename = os.path.join(os.getcwd(), output_filename)
        fig.savefig(output_filename, bbox_inches="tight")
        print(f"\n[SUCCESS] Generated PyVista paper-quality figure saved to: {output_filename}")
        plt.close(fig)
        return output_filename


def probability_fields_for(geo_model, inference_data, topography_path):
    import gempy_viewer as gpv
    online_prob = setup_probability_fields(geo_model, inference_data, topography_path)

    # 1. Generate the paper-quality, multi-panel summary figure
    generate_paper_quality_figure(geo_model, online_prob, topography_path)
    generate_pyvista_paper_quality_figure(geo_model, online_prob)

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
        threshold_kwargs={'value': [0.2, 0.9], 'invert': False},
        kwargs_pyvista_bounds={
                'show_xlabels': False,
                'show_ylabels': False,
                'show_zlabels': False,
        }
    )

    return online_prob


def setup_probability_fields(geo_model, inference_data, topography_path):
    online_prob = compute_probability_density_fields(
        geo_model,
        inference_data,
        n_samples=5
    )
    import gempy as gp

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
    return online_prob
