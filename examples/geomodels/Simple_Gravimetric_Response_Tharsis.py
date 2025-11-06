import sys
import os
import numpy as np
import gempy as gp
import geopandas as gpd
from mineye.config import paths
from mineye.config.example_parameters import TharsisGravityConfig
from mineye.GeoModel import helper_plotter
from examples.geomodels.Simple_Model_Tharsis import create_simple_model

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

def load_gravity_data():
    """Load gravity data from GeoJSON file."""
    gravity_data = gpd.read_file(paths.get_gravity_data_path())
    observed_gravity = gravity_data[TharsisGravityConfig.GRAVITY_FIELD_NAME].values  # in mGal
    return gravity_data, observed_gravity

def create_measurement_grid(geo_model, gravity_data):
    """Create measurement grid for gravity calculations."""
    extent = geo_model.grid.regular_grid.extent
    min_x, max_x, min_y, max_y, min_z, max_z = extent

    if TharsisGravityConfig.USE_ACTUAL_MEASUREMENT_LOCATIONS:
        print("Using actual gravity measurement locations...")
        xy_ravel = np.column_stack([
            np.array(gravity_data.geometry.x.values),
            np.array(gravity_data.geometry.y.values),
            np.full(len(gravity_data), max_z)
        ])
        print(f"Using {len(xy_ravel)} actual measurement points")
    else:
        print("Using regular grid for measurements...")
        grav_res = TharsisGravityConfig.GRAVITY_RESOLUTION
        X = np.linspace(min_x, max_x, grav_res)
        Y = np.linspace(min_y, max_y, grav_res)
        Z = max_z
        xyz = np.meshgrid(X, Y, Z)
        xy_ravel = np.vstack(list(map(np.ravel, xyz))).T
        print(f"Created regular grid: {grav_res}x{grav_res} = {len(xy_ravel)} points")

    return xy_ravel


def compute_forward_gravity_model(geo_model, xy_ravel):
    """Compute forward gravity model."""
    print("Computing forward gravity model...")

    # Set centered grid
    print("Setting up centered grid...")
    gp.set_centered_grid(
        grid=geo_model.grid,
        centers=xy_ravel,
        resolution=TharsisGravityConfig.CENTERED_GRID_RESOLUTION,
        radius=TharsisGravityConfig.CENTERED_GRID_RADIUS
    )

    # Calculate gravity gradient
    print("Calculating gravity gradient...")
    gravity_gradient = gp.calculate_gravity_gradient(geo_model.grid.centered_grid)

    # Configure geophysics input
    print("Configuring geophysics input...")
    geo_model.geophysics_input = gp.data.GeophysicsInput(
        tz=gravity_gradient,
        densities=np.array([
            TharsisGravityConfig.DENSITY_PLUTONITES,
            TharsisGravityConfig.DENSITY_SEDIMENTARY_HOST
        ])
    )

    # Compute forward model
    print("Computing forward model...")
    geo_model.interpolation_options.mesh_extraction = False
    sol = gp.compute_model(
        gempy_model=geo_model,
        engine_config=gp.data.GemPyEngineConfig(
            backend=gp.data.AvailableBackends.numpy,
            dtype='float32'
        )
    )

    print("Forward modeling complete!")
    return sol.gravity


def normalize_data(observed_ugal, forward_model, method):
    """Apply normalization to gravity data."""
    print(f"Applying {method} normalization...")

    if method == 'zscore':
        obs_mean, obs_std = np.mean(observed_ugal), np.std(observed_ugal)
        fwd_mean, fwd_std = np.mean(forward_model), np.std(forward_model)
        observed_norm = (observed_ugal - obs_mean) / obs_std
        forward_norm = (forward_model - fwd_mean) / fwd_std
        unit_label = 'Z-score'

    elif method == 'minmax':
        obs_min, obs_max = np.min(observed_ugal), np.max(observed_ugal)
        fwd_min, fwd_max = np.min(forward_model), np.max(forward_model)
        observed_norm = (observed_ugal - obs_min) / (obs_max - obs_min)
        forward_norm = (forward_model - fwd_min) / (fwd_max - fwd_min)
        unit_label = 'Normalized [0-1]'

    elif method == 'mean_center':
        observed_norm = observed_ugal - np.mean(observed_ugal)
        forward_norm = forward_model - np.mean(forward_model)
        unit_label = 'Mean-centered (μGal)'

    elif method == 'relative':
        obs_range = np.max(observed_ugal) - np.min(observed_ugal)
        fwd_range = np.max(forward_model) - np.min(forward_model)
        observed_norm = observed_ugal / obs_range
        forward_norm = forward_model / fwd_range
        unit_label = 'Relative to range'

    else:
        raise ValueError(f"Unknown normalization method: {method}")

    print(f"  Observed stats (normalized): mean={np.mean(observed_norm):.3f}, std={np.std(observed_norm):.3f}")
    print(f"  Forward stats (normalized):  mean={np.mean(forward_norm):.3f}, std={np.std(forward_norm):.3f}")

    return observed_norm, forward_norm, unit_label

def test_gravity_forward_model():
    """Test gravity forward model using Tharsis data."""
    # Create geological model
    geo_model = create_simple_model()

    # Load gravity data
    gravity_data, observed_gravity = load_gravity_data()

    # Create measurement grid
    xy_ravel = create_measurement_grid(geo_model, gravity_data)

    # Compute forward model
    grav = compute_forward_gravity_model(geo_model, xy_ravel)

    # Plot forward model
    if TharsisGravityConfig.SHOW_FORWARD_MODEL:
        helper_plotter.plot_forward_gravity_model(
            xy_ravel,
            grav,
            gravity_data,
            gravity_resolution=TharsisGravityConfig.GRAVITY_RESOLUTION,
            use_actual_locations=TharsisGravityConfig.USE_ACTUAL_MEASUREMENT_LOCATIONS
        )

    # Comparison plots
    if TharsisGravityConfig.SHOW_COMPARISON_PLOTS and TharsisGravityConfig.USE_ACTUAL_MEASUREMENT_LOCATIONS:
        print("\n=== Observed vs Predicted Comparison ===")

        # Convert units: observed is in mGal, predicted in μGal
        observed_ugal = observed_gravity * 1000
        forward_model = grav.copy()

        # Apply normalization if enabled
        if TharsisGravityConfig.NORMALIZE_DATA:
            observed_norm, forward_norm, unit_label = normalize_data(
                observed_ugal, forward_model, TharsisGravityConfig.NORMALIZATION_METHOD
            )
        else:
            observed_norm = observed_ugal
            forward_norm = forward_model
            unit_label = 'μGal'

        residuals_norm = observed_norm - forward_norm

        # Plot comparison
        helper_plotter.plot_gravity_comparison(
            xy_ravel,
            observed_norm,
            forward_norm,
            residuals_norm,
            unit_label,
            normalize_data=TharsisGravityConfig.NORMALIZE_DATA,
            normalization_method=TharsisGravityConfig.NORMALIZATION_METHOD
        )

