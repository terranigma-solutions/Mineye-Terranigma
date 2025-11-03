from functools import partial
from typing import Any

import arviz
import numpy as np
from pyro.distributions import Distribution

import gempy as gp
from gempy_engine.core.backend_tensor import BackendTensor
from gempy_engine.core.data.interpolation_input import InterpolationInput
from gempy_probability.modules.plot.plot_gempy import plot_gempy
import gempy_viewer as gpv


def plot(gravity_samples_norm, observed_gravity_ugal, xy_ravel) -> tuple[str, Any]:
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


def baseline(geo_model) -> Any:
    print("\n" + "=" * 60)
    print("COMPUTING BASELINE FORWARD MODEL")
    print("=" * 60)
    print("Computing gravity with mean/initial prior parameters...")

    baseline_fw_gravity = geo_model.solutions.gravity

    # Convert to numpy for normalization parameter computation
    if hasattr(baseline_fw_gravity, 'numpy'):
        baseline_fw_gravity_np = baseline_fw_gravity.numpy()
    else:
        baseline_fw_gravity_np = np.array(baseline_fw_gravity)

    print(f"Baseline forward model statistics:")
    print(f"  Mean: {np.mean(baseline_fw_gravity_np):.2f} μGal")
    print(f"  Std:  {np.std(baseline_fw_gravity_np):.2f} μGal")
    print(f"  Min:  {np.min(baseline_fw_gravity_np):.2f} μGal")
    print(f"  Max:  {np.max(baseline_fw_gravity_np):.2f} μGal")
    print("=" * 60 + "\n")
    return -baseline_fw_gravity_np


def normalize(baseline_fw_gravity_np, observed_gravity):
    # Convert observed gravity from mGal to μGal for comparison
    from mineye.GeoModel.geophysics import compute_alignment_params
    norm_params = compute_alignment_params(
        observed=observed_gravity,
        baseline_forward=baseline_fw_gravity_np
    )

    return norm_params


def setup_geomodel(gravity_data, simple_geo_model: gp.data.GeoModel):
    geo_model: gp.data.GeoModel = simple_geo_model
    BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

    xy_ravel = np.column_stack([
            np.array(gravity_data.geometry.x.values),
            np.array(gravity_data.geometry.y.values),
            np.full(len(gravity_data), 0)  # Set Z to surface elevation
    ])

    _gravity_precomputations(density_plutonites=2.9, density_sedimentary_host=2.3, xy_ravel=xy_ravel, simple_geo_model=geo_model)

    geo_model.interpolation_options.mesh_extraction = False

    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.CENTERED],
        reset=True
    )
    gp.compute_model(geo_model)
    return geo_model, xy_ravel


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


def gempy_viz(geo_model: gp.data.GeoModel, prior_inference_data: arviz.InferenceData):
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
        show=False
    )

    plot_gempy(
        geo_model=geo_model,
        n_samples=20,
        samples=(prior_inference_data.prior[r'dips'].values[0, :]),
        update_model_fn=_update_model_for_plotting,
        gempy_plot=p2d
    )


def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: float, sample_idx: int):
    # # Modify the surface point
    gp.modify_orientations(
        geo_model=geo_model,
        dip=sample_value,
    )


def create_orientation_modifier(key: str):
    """Factory function that creates orientation modifier functions for different keys."""
    return partial(modify_orientations, key=key)


def modify_orientations(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
        key: str
) -> InterpolationInput:
    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame
    from gempy.modules.data_manipulation.manipulate_points import (
        compute_adp_from_gradients,
        convert_orientation_to_pole_vector
    )

    interp_input: InterpolationInput = interpolation_input_from_structural_frame(geo_model)
    samples_value = samples[key]

    og_gradients = interp_input.orientations.dip_gradients
    azimuth, dip, polarity = compute_adp_from_gradients(
        G_x=og_gradients[:, 0],
        G_y=og_gradients[:, 1],
        G_z=og_gradients[:, 2],
    )

    gradients = convert_orientation_to_pole_vector(
        azimuth=azimuth,
        dip=samples_value,
        polarity=polarity
    )

    interp_input.orientations.dip_gradients = gradients
    return interp_input

