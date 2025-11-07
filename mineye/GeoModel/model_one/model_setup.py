import os
from typing import Any

import numpy as np

import gempy as gp
from gempy_engine.core.backend_tensor import BackendTensor

import geopandas as gpd


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


def setup_geomodel(gravity_data, simple_geo_model: gp.data.GeoModel):
    geo_model: gp.data.GeoModel = simple_geo_model
    BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

    xy_ravel = np.column_stack([
            np.array(gravity_data.geometry.x.values),
            np.array(gravity_data.geometry.y.values),
            np.full(len(gravity_data), 0)  # Set Z to surface elevation
    ])

    _gravity_precomputations(density_plutonites=2.9, density_sedimentary_host=2.3, xy_ravel=xy_ravel, simple_geo_model=geo_model)
    import torch
    geo_model.geophysics_input.tz = torch.tensor(geo_model.geophysics_input.tz)
    geo_model.interpolation_options.mesh_extraction = False

    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.CENTERED],
        reset=True
    )
    gp.compute_model(geo_model)
    return geo_model, xy_ravel


def read_gravity(geophysical_dir):
    gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
    observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

    # Take a spatially distributed subset of measurements
    n_points = 20  # Adjust this number to control how many points you want
    xy_coords = gravity_data.geometry.apply(lambda p: (p.x, p.y)).to_list()
    xy_array = np.array(xy_coords)

    # Use K-means clustering to get well-distributed points
    from sklearn.cluster import KMeans
    kmeans = KMeans(n_clusters=n_points, random_state=42)
    kmeans.fit(xy_array)

    # Find the closest points to cluster centers
    from scipy.spatial.distance import cdist
    centers = kmeans.cluster_centers_
    distances = cdist(centers, xy_array)
    indices = [np.argmin(dist) for dist in distances]

    # Filter the gravity data
    observed_gravity = observed_gravity[indices]
    gravity_data = gravity_data.iloc[indices]

    observed_gravity_ugal = observed_gravity * 1000
    return gravity_data, observed_gravity_ugal


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
        densities=np.array([density_plutonites, density_sedimentary_host])  # kg/m³ for different formations,
    )
