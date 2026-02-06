import os

import arviz as az
import matplotlib.pyplot as plt
import numpy as np
import pytest
import torch
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality, check_likelihood_balance
from mineye.GeoModel.model_one.joint_probabilistic_model import joint_set_priors
from mineye.GeoModel.model_one.model_setup import setup_geomodel, read_gravity, baseline
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_hierarchical_per_station, enmap_likelihood_fn
from mineye.GeoModel.model_one.visualization import (generate_gravity_uncertainty_plots,
                                                     gempy_viz,
                                                     plot_many_observed_vs_forward,
                                                     probability_fields_for, plot_probability_heatmap, plot_heat_map)
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison
# noinspection PyUnusedImports
from tests.tests_inversions.conftest import simple_geo_model, topography_dir, base_dir, geophysical_dir


class TestJointInversion:
    
    def test_likelihood_balance(self, simple_geo_model, geophysical_dir, base_dir, n_samples=50):
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        geo_model.interpolation_options.sigmoid_slope = 100

        # 2. Setup EnMap Data
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')
        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")

        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        labels_enmap[labels_enmap == 2] = 1  # Normalize labels

        simple_geo_model.interpolation_options.mesh_extraction = False
        simple_geo_model.interpolation_options.evaluation_options.number_octree_levels = 1
        gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=False
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
                'dips'   : dist.Normal(
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
        # likelihood_fn = generate_joint_likelihood(norm_params)
        gravity_dist = generate_multigravity_likelihood_hierarchical_per_station(norm_params)
        enmap_dist = enmap_likelihood_fn

        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=joint_set_priors,
            likelihood_fn=[gravity_dist, enmap_dist],
            obs_name="Joint_Obs"
        )
        # 7. Prepare observed data
        gravity_obs_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
        enmap_obs_tensor = torch.tensor(labels_enmap, dtype=torch.float64)
        check_likelihood_balance(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list =  [gravity_obs_tensor, enmap_obs_tensor]
        )


    
    def test_joint_inversion(self, simple_geo_model, geophysical_dir, base_dir,
                             arviz_data_filename="arviz_data_joint_Feb05_2026_v2.nc"):
        """Test joint inversion of gravity and EnMap data."""

        # 1. Setup Gravity Data
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        geo_model.interpolation_options.sigmoid_slope = 100

        # 2. Setup EnMap Data
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')
        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")

        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        labels_enmap[labels_enmap == 2] = 1  # Normalize labels

        simple_geo_model.interpolation_options.mesh_extraction = False
        simple_geo_model.interpolation_options.evaluation_options.number_octree_levels = 1
        gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=False
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
                'dips'   : dist.Normal(
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
        # likelihood_fn = generate_joint_likelihood(norm_params)
        gravity_dist = generate_multigravity_likelihood_hierarchical_per_station(norm_params)
        enmap_dist = enmap_likelihood_fn

        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=joint_set_priors,
            likelihood_fn=[gravity_dist, enmap_dist],
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
            n_samples=200,
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
                num_samples=1000,
                warmup_steps=1000,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        data.extend(prior_data)
        save_path = os.path.join(os.path.dirname(__file__), arviz_data_filename)
        data.to_netcdf(save_path)
        print(f"Saved joint inversion results to {save_path}")


    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        # Run comprehensive diagnostics with plots
        check_mcmc_quality(data, verbose=True, plot=True)

    def test_run_predictive_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        az.plot_trace(data.prior)
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)

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
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

        gempy_viz(geo_model, data, n_samples=20, ve=3)

    def test_run_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))

        az.plot_posterior(data, var_names=["dips", "density"])
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

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

        if hasattr(data, 'posterior') and response in data.prior:
            gravity_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.posterior_predictive[response].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        # * 9) Analysis Gempy Model
        gempy_viz(simple_geo_model, data)

    def test_probability_plots(self, simple_geo_model, geophysical_dir, base_dir, topography_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))


        topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')

        probability_fields_for(simple_geo_model, data.prior, topography_path)
        probability_fields_for(simple_geo_model, data.posterior, topography_path)


    def test_check_trace(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_Feb05_2026.nc"))
        az.plot_trace(data.posterior)
        plt.show()

        plot_probability_heatmap(data, 'prior')
        plot_probability_heatmap(data, 'posterior_predictive')

        plot_heat_map(
            group='posterior_predictive',
            heatmap_data=(data.posterior_predictive['probs_pred'].isel(chain=0, draw=-1).values.T)
        )

        plot_heat_map(
            group='posterior_predictive',
            heatmap_data=(data.posterior_predictive['probs_pred'].isel(chain=0, draw=-10).values.T)
        )

        plot_heat_map(
            group='posterior_predictive',
            heatmap_data=(data.posterior_predictive['probs_pred'].isel(chain=0, draw=-20).values.T)
        )
