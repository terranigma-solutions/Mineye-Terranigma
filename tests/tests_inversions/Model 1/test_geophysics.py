import os
import time

from matplotlib import pyplot as plt

import gempy as gp
import gempy_viewer as gpv

import numpy as np
import geopandas as gpd

from mineye.GeoModel.plotting.probabilistic_analysis import plot_comparison, plot_fw_geophysics


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
    
    gpv.plot_2d(simple_geo_model)
    gpv.plot_3d(simple_geo_model, ve=5, image=True)


    if PLOT:=True:
        # For actual measurement locations, show as scatter plot with color-coded gravity values
        plot_fw_geophysics(
            fw_values=grav,
            observed_data=gravity_data,
            xy_ravel=xy_ravel,
            label=r'Forward Model Gravity ($\mu$gal)',
            title='Forward Model Gravity Results'
        )
        plot_comparison(observed_gravity, grav, xy_ravel)
        
        
def test_simple_model_magnetics(simple_geo_model, geophysical_dir):
    filtered_magnetic = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson'))
    # Choose 20 random values
    sampled_magnetic = filtered_magnetic.sample(n=20)
    
    
    observed_magnetics = filtered_magnetic['MAG']

    # Create two subplots side by side
    fig, (ax1) = plt.subplots(1, 1, figsize=(24, 10))

    # Plot the filtered magnetic data
    filtered_magnetic.plot(
        column='MAG',
        ax=ax1,
        cmap='RdYlBu_r',
        legend=True,
        legend_kwds={'label': 'Magnetic Intensity (nT)', 'orientation': 'vertical'}
    )

    ax1.set_title('Filtered Magnetic Data - Tharsis AOI 1')
    ax1.set_xlabel('Easting (m)')
    ax1.set_ylabel('Northing (m)')
    ax1.set_aspect('equal')

    print(f"Magnetic data loaded successfully!")
    print(f"Data shape: {filtered_magnetic.shape}")
    print(f"Columns: {filtered_magnetic.columns.tolist()}")
    print(f"Magnetic range: {filtered_magnetic['MAG'].min():.2f} to {filtered_magnetic['MAG'].max():.2f} nT")
    plt.show()

    xy_ravel = np.column_stack([
            np.array(sampled_magnetic.geometry.x.values),
            np.array(sampled_magnetic.geometry.y.values),
            np.full(len(sampled_magnetic), 0)  # Set Z to surface elevation
    ])


    if PLOT:=True:
        # For actual measurement locations, show as scatter plot with color-coded gravity values
        fw_magnetics = np.zeros(len(xy_ravel))
        # plot_fw_geophysics(fw_magnetics, sampled_magnetic, xy_ravel)
        # plot_comparison(sampled_magnetic['MAG'].values, fw_magnetics, xy_ravel)


        plot_fw_geophysics(
            fw_values=fw_magnetics,
            observed_data=sampled_magnetic,
            xy_ravel=xy_ravel,
            label=r'Forward Model Magnetics (nT)',
            title='Forward Model Magnetics Results'
        )
        plot_comparison(
            observed=sampled_magnetic['MAG'].values,
            fw_values=fw_magnetics,
            xy_ravel=xy_ravel
        )


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
