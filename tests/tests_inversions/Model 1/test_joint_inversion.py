import os
import torch
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import pyro
from pyro import distributions as dist
import pytest

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig

from mineye.GeoModel.model_one.model_setup import setup_geomodel, read_gravity, baseline
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.joint_probabilistic_model import joint_set_priors, generate_joint_likelihood

from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality
from mineye.GeoModel.model_one.visualization import (generate_gravity_uncertainty_plots, 
                                                     gempy_viz, 
                                                     plot_many_observed_vs_forward,
                                                     plot_joint_inversion_results,
                                                     compute_probability_density_fields)
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue

# noinspection PyUnusedImports
from tests.tests_inversions.conftest import simple_geo_model, topography_dir, base_dir, geophysical_dir

class TestJointInversion:
    def test_joint_inversion(self, simple_geo_model, geophysical_dir, base_dir, n_samples=50,
                             arviz_data_filename="arviz_data_joint_Feb05_2026.nc"):
        """Test joint inversion of gravity and EnMap data."""
        
        # 1. Setup Gravity Data
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        
        # 2. Setup EnMap Data
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')
        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")
        
        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        labels_enmap[labels_enmap == 2] = 1  # Normalize labels
        
        # 3. Setup GeoModel with both grids (GemPy supports one active grid, so we might need to be careful)
        # Usually, for joint inversion, we want both observations at their respective locations.
        # In GemPy, we can use a CustomGrid for EnMap and the standard Grid for Gravity if it's on a grid,
        # but here gravity is at specific station locations (which are also treated as a custom grid in GemPy usually).
        
        # In test_gravity_inversion.py: geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        # setup_geomodel sets the custom grid to the gravity station locations.
        
        # If we have two sets of custom locations, we should combine them into one custom grid.
        
        # Combine EnMap points and Gravity points?
        # Actually, EnMap likelihood needs the scalar field at xyz_enmap.
        # Gravity likelihood needs the gravity at gravity_data locations.
        
        # Let's combine all custom points:
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        
        combined_custom_points = np.vstack([gravity_xyz, xyz_enmap])
        
        # Now set the custom grid
        gp.set_custom_grid(simple_geo_model.grid, combined_custom_points)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=True
        )
        
        # We need to tell the likelihood functions where their data is in the custom grid.
        # Gravity is first len(gravity_xyz) points.
        # EnMap is the rest.
        
        # 4. Normalization for Gravity
        norm_params = normalize(
            baseline_fw_gravity_np=(baseline(simple_geo_model)),
            observed_gravity=observed_gravity_ugal,
            method="align_to_reference",
            extrapolation_buffer=0.3
        )
        
        # 5. Define Priors
        model_priors = {
            'dips': dist.Normal(
                loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
                scale=torch.tensor(10, dtype=torch.float64),
                validate_args=True
            ),
            'density': dist.Normal(
                loc=(torch.tensor([2.9 - 2.67, 2.3 - 2.67])),
                scale=torch.tensor(0.15),
            ).to_event(1)
        }
        
        # 6. Create Probabilistic Model
        likelihood_fn = generate_joint_likelihood(norm_params, n_gravity_points=len(gravity_xyz))
        
        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=joint_set_priors,
            likelihood_fn=likelihood_fn,
            obs_name="Joint_Obs"
        )
        
        # 7. Prepare observed data
        gravity_obs_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
        enmap_obs_tensor = torch.tensor(labels_enmap, dtype=torch.float64)
        
        # Joint observations as a list (supported by our modified make_gempy_pyro_model)
        joint_obs = [gravity_obs_tensor, enmap_obs_tensor]
        
        # 8. Run Prior Predictive
        print("Running joint prior predictive...")
        prior_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=simple_geo_model,
            y_obs_list=joint_obs,
            n_samples=20,
            plot_trace=True
        )
        
        # 9. Run MCMC Inference
        print("Running joint NUTS inference...")
        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=simple_geo_model,
            y_obs_list=joint_obs,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                num_samples=100,
                warmup_steps=100,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )
        
        data.extend(prior_data)
        save_path = os.path.join(os.path.dirname(__file__), arviz_data_filename)
        data.to_netcdf(save_path)
        print(f"Saved joint inversion results to {save_path}")

        # 10. Visualization
        plot_joint_inversion_results(
            data=data,
            observed_gravity=observed_gravity_ugal,
            xy_gravity=gravity_xyz[:, :2],  # xy_ravel expects (n, 2)
            observed_enmap=labels_enmap,
            geo_model=simple_geo_model,
            n_gravity_points=len(gravity_xyz)
        )
        plt.show()

    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        # Run comprehensive diagnostics with plots
        check_mcmc_quality(data, verbose=True, plot=True)

    def test_run_predictive_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        az.plot_trace(data.prior)
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        
        # We need the gravity locations for plotting
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        xy_ravel = gravity_xyz[:, :2]

        # Prepare data
        observed_norm = observed_gravity_ugal
        response = r'$\mu_{gravity}$'
        forward_norm = data.prior[response].mean(axis=1)
        
        # Handle different structures of prior/posterior_predictive if needed
        if hasattr(data, 'prior') and response in data.prior:
            many_forward_norm = data.prior[response].values[0, -20:]
            plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)

    def test_run_kde_sections(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        # Reconstruct the geomodel setup as in test_joint_inversion
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        xyz_enmap = np.load(xyz_path)
        combined_custom_points = np.vstack([gravity_xyz, xyz_enmap])
        
        gp.set_custom_grid(simple_geo_model.grid, combined_custom_points)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=True
        )

        gempy_viz(simple_geo_model, data, n_samples=100, ve=3)

    def test_run_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))

        az.plot_posterior(data, var_names=["dips", "density"])
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        xy_ravel = gravity_xyz[:, :2]

        plt.rcParams['figure.dpi'] = 72  # Lower DPI for faster rendering

        axes = az.plot_density(
            data=[data, data.prior],
            var_names=["dips", "density"],
            filter_vars="like",
            hdi_prob=0.9999,
            shade=.2,
            data_labels=["Posterior", "Prior"],
            colors=[default_red, default_blue],
        )
        plt.show()
        
        response = r'$\mu_{gravity}$'
        if hasattr(data, 'prior') and response in data.prior:
            plot_geophysics_comparison(forward_norm=data.prior[response].mean(axis=1), 
                                       normalization_method='align_to_reference', 
                                       observed_ugal=observed_gravity_ugal, 
                                       xy_ravel=xy_ravel)

        if hasattr(data.posterior_predictive, response):
            plot_geophysics_comparison(forward_norm=data.posterior_predictive[response].mean(axis=1), 
                                       normalization_method='align_to_reference', 
                                       observed_ugal=observed_gravity_ugal, 
                                       xy_ravel=xy_ravel)

        # * 9) Analysis inference
        if hasattr(data, 'prior') and response in data.prior:
            generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[response].values[0, :],
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        # * 9) Analysis Gempy Model
        gempy_viz(simple_geo_model, data)

    def test_probability_plots(self, simple_geo_model, geophysical_dir, base_dir, topography_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))

        # Reconstruct the geomodel setup as in test_joint_inversion
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        xyz_enmap = np.load(xyz_path)
        combined_custom_points = np.vstack([gravity_xyz, xyz_enmap])
        
        gp.set_custom_grid(simple_geo_model.grid, combined_custom_points)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=True
        )

        self._probability_fields_for(simple_geo_model, data.prior, topography_dir)
        self._probability_fields_for(simple_geo_model, data.posterior, topography_dir)

    @staticmethod
    def _probability_fields_for(geo_model, inference_data, topography_dir):
        online_prob = compute_probability_density_fields(
            geo_model,
            inference_data,
            n_samples=20
        )
        import gempy_viewer as gpv

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.probability_field[0],
            show_data=True,
            ve=5,
            kwargs_lithology={'cmap': 'viridis', 'norm': None}
        )

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.probability_field[1],
            show_data=True,
            ve=5,
            kwargs_lithology={'cmap': 'viridis', 'norm': None}
        )

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.entropy,
            show_data=True,
            ve=5,
            kwargs_lithology={'cmap': 'magma', 'norm': None}
        )

        topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')
        gp.set_topography_from_file(
            grid=geo_model.grid,
            filepath=topography_path,
            crop_to_extent=[
                geo_model.grid.extent[0],
                geo_model.grid.extent[2],
                geo_model.grid.extent[1],
                geo_model.grid.extent[3]
            ]
        )

        gp.compute_model(geo_model)

        # * We inject the entropy field into the model
        geo_model.solutions.raw_arrays.scalar_field_matrix[0] = online_prob.entropy
        gpv.plot_3d(
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
