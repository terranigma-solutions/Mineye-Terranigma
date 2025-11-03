from typing import Any

import numpy as np

import gempy as gp
from gempy_engine.core.backend_tensor import BackendTensor


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


def setup_geomodel(self, gravity_data, simple_geo_model: gp.data.GeoModel):
    geo_model: gp.data.GeoModel = simple_geo_model
    BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

    xy_ravel = np.column_stack([
            np.array(gravity_data.geometry.x.values),
            np.array(gravity_data.geometry.y.values),
            np.full(len(gravity_data), 0)  # Set Z to surface elevation
    ])

    self.gravity_precomputations(density_plutonites=2.9, density_sedimentary_host=2.3, xy_ravel=xy_ravel, simple_geo_model=geo_model)

    geo_model.interpolation_options.mesh_extraction = False

    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.CENTERED],
        reset=True
    )
    gp.compute_model(geo_model)
    return geo_model, xy_ravel


def gravity_precomputations(density_plutonites: float, density_sedimentary_host: float,
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
