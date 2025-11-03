import os
import time

import gempy as gp
import gempy_viewer as gpv

import numpy as np
import geopandas as gpd

from mineye.GeoModel.plotting.probabilistic_analysis import _plot_comparison, _plot_fw_gravity


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
