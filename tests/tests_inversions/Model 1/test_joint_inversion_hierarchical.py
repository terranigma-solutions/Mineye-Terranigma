import os
import pytest
import numpy as np
import torch
import arviz as az
import matplotlib.pyplot as plt
from pyro import distributions as dist
import pyro

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue

# noinspection PyUnusedImports
from tests.tests_inversions.conftest import simple_geo_model, topography_dir, base_dir, geophysical_dir
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations, _modify_densities
from mineye.GeoModel.model_one.model_setup import setup_geomodel, read_gravity, baseline
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import _get_ordinal_probs
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality
from mineye.GeoModel.model_one.visualization import (gempy_viz,
                                                     plot_many_observed_vs_forward,
                                                     compute_probability_density_fields,
                                                     generate_gravity_uncertainty_plots,
                                                     plot_joint_inversion_results)
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison

def set_priors_joint(samples, geo_model):
    """Modify both orientations and densities for joint inversion."""
    interp_input = _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key="dips"
    )
    
    _modify_densities(
        samples=samples,
        geo_model=geo_model,
        key="density"
    )
    return interp_input

class TestJointInversionHierarchical:
    def test_run_predictive(self, simple_geo_model, base_dir, geophysical_dir, n_samples=50):
        geo_model, observed_joint, prob_model, xy_ravel, norm_params = \
            self._create_probabilistic_model(base_dir, geophysical_dir, simple_geo_model)

        # Run predictive
        print("Running joint prior predictive...")
        
        prior_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=observed_joint,
            n_samples=n_samples,
            plot_trace=True
        )

    def test_joint_inversion(self, simple_geo_model, base_dir, geophysical_dir, n_samples=50,
                             arviz_data_filename="arviz_data_joint_hierarchical.nc"):
        """Test Joint EnMap and Hierarchical Gravity inversion."""
        geo_model, observed_joint, prob_model, xy_ravel, norm_params = \
            self._create_probabilistic_model(base_dir, geophysical_dir, simple_geo_model)

        # Run predictive
        print("Running joint prior predictive...")
        prior_inference_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=observed_joint,
            n_samples=n_samples,
            plot_trace=True
        )

        # Run inference (NUTS)
        print("Running joint NUTS inference...")
        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=observed_joint,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                num_samples=100,
                warmup_steps=100,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        data.extend(prior_inference_data)
        data.to_netcdf(os.path.join(os.path.dirname(__file__), arviz_data_filename))

        # Add initial visualization
        gravity_data, _ = read_gravity(geophysical_dir)
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_enmap = np.load(os.path.join(base_dir, 'central_labels.npy'))
        labels_enmap[labels_enmap == 2] = 1

        plot_joint_inversion_results(
            data=data,
            observed_gravity=observed_joint[labels_enmap.shape[0]:].numpy(),
            xy_gravity=gravity_xyz[:, :2],
            observed_enmap=labels_enmap,
            geo_model=geo_model,
            n_gravity_points=len(gravity_data)
        )
        plt.show()

    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))
        # Run comprehensive diagnostics with plots
        check_mcmc_quality(data, verbose=True, plot=True)

    def test_run_predictive_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))
        az.plot_trace(data.prior)
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        xy_ravel = gravity_xyz[:, :2]

        # Prepare data
        observed_norm = observed_gravity_ugal
        response = r'$\mu_{gravity}$'
        forward_norm = data.prior[response].mean(axis=1)
        
        if hasattr(data, 'prior') and response in data.prior:
            many_forward_norm = data.prior[response].values[0, -20:]
            plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)

    def test_run_kde_sections(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))
        
        gravity_data, _ = read_gravity(geophysical_dir)
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

    def test_run_outlier_detection(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))

        if "sigma_stations" in data.posterior_predictive:
            posterior_sigmas = data.posterior_predictive["sigma_stations"].values  
            station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
            sigma_global_mean = station_noise_mean.mean()
            problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
            print(f"Check stations: {problematic}")
            
            az.plot_density(
                data=[data, data.prior],
                var_names=["sigma_stations"],
                filter_vars="like",
                hdi_prob=0.9999,
                shade=.2,
                data_labels=["Posterior", "Prior"],
                colors=[default_red, default_blue],
            )
            plt.show()

    def test_run_analysis(self, simple_geo_model, geophysical_dir, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))

        az.plot_posterior(data, var_names=["dips", "density"])
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        gravity_xyz = np.zeros((len(gravity_data), 3))
        gravity_xyz[:, 0] = gravity_data['X'].values
        gravity_xyz[:, 1] = gravity_data['Y'].values
        gravity_xyz[:, 2] = gravity_data['Z'].values
        xy_ravel = gravity_xyz[:, :2]

        plt.rcParams['figure.dpi'] = 72  

        az.plot_density(
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

        if hasattr(data, 'prior') and response in data.prior:
            generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[response].values[0, :],
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        gempy_viz(simple_geo_model, data)

    def test_probability_plots(self, simple_geo_model, geophysical_dir, base_dir, topography_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_joint_hierarchical.nc"))

        gravity_data, _ = read_gravity(geophysical_dir)
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

    def _create_probabilistic_model(self, base_dir, geophysical_dir, simple_geo_model):
        # 1. Load EnMap data
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')

        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found.")

        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        labels_enmap[labels_enmap == 2] = 1  # Normalize labels as in enmap test
        observed_labels = torch.tensor(labels_enmap, dtype=torch.float64)

        # 2. Load Gravity data
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        observed_gravity = torch.tensor(observed_gravity_ugal, dtype=torch.float64)

        # 3. Setup Geomodel with Gravity (sets centered grid)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        
        # 4. Add EnMap Custom Grid
        gp.set_custom_grid(geo_model.grid, xyz_enmap)
        
        # 5. Normalize Gravity
        norm_params = normalize(
            baseline_fw_gravity_np=(baseline(geo_model)),
            observed_gravity=observed_gravity_ugal,
            method="align_to_reference",
            extrapolation_buffer=0.3
        )

        # 6. Define Priors
        model_priors = {
            'dips': dist.Normal(
                loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),
                scale=torch.tensor(10, dtype=torch.float64),
                validate_args=True
            ),
            'density': dist.Normal(
                loc=(torch.tensor([2.9 - 2.67, 2.3 - 2.67])),
                scale=torch.tensor(0.15),
            ).to_event(1)
        }

        # 7. Define Joint Likelihood
        def joint_likelihood_fn(solutions):
            # --- A. EnMap Likelihood Part ---
            output_center = solutions.octrees_output[0].last_output_center
            scalar_field_at_custom_grid = output_center.exported_fields.scalar_field[output_center.grid.custom_grid_slice]
            if not isinstance(scalar_field_at_custom_grid, torch.Tensor):
                scalar_field_at_custom_grid = torch.tensor(scalar_field_at_custom_grid, dtype=torch.float64)
            
            boundaries = torch.tensor(solutions.scalar_field_at_surface_points, dtype=torch.float64)
            probs = _get_ordinal_probs(scalar_field_at_custom_grid, boundaries, temperature=0.1)
            pyro.deterministic("probs_pred", probs.detach())
            
            # --- B. Gravity Likelihood Part (Hierarchical) ---
            simulated_gravity = align_forward_to_observed(-solutions.gravity, norm_params)
            pyro.deterministic(r"$\mu_{gravity}$", simulated_gravity.detach())
            n_stations = simulated_gravity.shape[0]

            mu_log_sigma = pyro.sample(
                "mu_log_sigma",
                dist.Normal(
                    torch.tensor(np.log(5000.0), dtype=torch.float64),
                    torch.tensor(0.5, dtype=torch.float64)
                )
            )
            tau_log_sigma = pyro.sample(
                "tau_log_sigma",
                dist.HalfNormal(torch.tensor(0.5, dtype=torch.float64))
            )
            log_sigma_stations = pyro.sample(
                "log_sigma_stations",
                dist.Normal(mu_log_sigma.expand([n_stations]), tau_log_sigma).to_event(1)
            )
            sigma_stations = torch.exp(log_sigma_stations)
            pyro.deterministic("sigma_stations", sigma_stations)

            # --- C. Combine into a Joint Distribution ---
            # We use two independent pyro.sample calls inside the likelihood_fn?
            # NO, the likelihood_fn in GemPyPyroModel is expected to return ONE distribution.
            # So we concatenate.
            
            # For Categorical, we need to represent the EnMap part.
            # But the joint needs to be one distribution. 
            # This is tricky because one is Categorical and the other is Normal.
            
            # Alternative: Use pyro.sample directly for each part and return a Dummy distribution or None?
            # GemPyPyroModel calls likelihood_fn and then pyro.sample("obs", dist, obs=y_obs).
            
            # If we want to use different distribution types, we might need to bypass the standard return.
            # But y_obs is passed as a single tensor.
            
            # If we concatenate, they must be the same distribution type (e.g. Normal).
            # EnMap labels are discrete, but we can model them as Normal with small sigma if they are floats.
            # OR we can manually call pyro.sample for each and return None.
            
            # Let's see how GemPyPyroModel handles the return.
            # If we return None, it might fail.
            
            # If we use Normal for both:
            enmap_sigma = torch.full_like(observed_labels, 0.1) # Small sigma for "labels"
            
            joint_loc = torch.cat([scalar_field_at_custom_grid, simulated_gravity]) # Wait, EnMap needs labels not scalar field
            # Actually, enmap_likelihood_fn in the original file uses dist.Categorical(probs=probs).
            # If we want to use Joint, we have to be consistent.
            
            # Let's try to just use pyro.sample for both and return a "dummy" distribution
            pyro.sample("obs_enmap", dist.Categorical(probs=probs), obs=observed_labels)
            pyro.sample("obs_gravity", dist.Normal(simulated_gravity, sigma_stations).to_event(1), obs=observed_gravity)
            
            return None # This might break GemPyPyroModel if it expects a return

        # Actually, looking at gempy_probability, it expects a distribution.
        # Let's concatenate them and use Normal for EnMap labels as well, 
        # but that's not ideal for Categorical.
        
        # What if we use a Custom distribution? 
        # Better: use a joint likelihood that returns a Normal where the EnMap part is the "expected label" (rounded scalar field? no).
        
        # Re-evaluating: Original enmap_likelihood_fn:
        # return dist.Categorical(probs=probs).to_event(1)
        
        # If I want to combine them, and y_obs_list is ONE tensor...
        # I'll use the Normal approximation for EnMap for now to keep it simple and compatible with concatenation.
        # Or better, I will use a joint likelihood that samples both.
        
        def joint_likelihood_fn_v2(solutions):
            # EnMap
            output_center = solutions.octrees_output[0].last_output_center
            scalar_field_at_custom_grid = output_center.exported_fields.scalar_field[output_center.grid.custom_grid_slice]
            boundaries = torch.tensor(solutions.scalar_field_at_surface_points, dtype=torch.float64)
            probs = _get_ordinal_probs(scalar_field_at_custom_grid, boundaries, temperature=0.1)
            
            # Gravity
            simulated_gravity = align_forward_to_observed(-solutions.gravity, norm_params)
            n_stations = simulated_gravity.shape[0]
            mu_log_sigma = pyro.sample("mu_log_sigma", dist.Normal(torch.tensor(np.log(5000.0), dtype=torch.float64), 0.5))
            tau_log_sigma = pyro.sample("tau_log_sigma", dist.HalfNormal(0.5))
            log_sigma_stations = pyro.sample("log_sigma_stations", dist.Normal(mu_log_sigma.expand([n_stations]), tau_log_sigma).to_event(1))
            sigma_stations = torch.exp(log_sigma_stations)
            
            # Since we must return ONE distribution that matches y_obs_list
            # We'll use Normal for both. For EnMap, we use the predicted "class" (argmax of probs) or just labels.
            # Actually, to keep it differentiable, we can use the mean of the categorical as a float label? No.
            
            # Let's use the simplest approach: Normal distribution for both.
            # For EnMap, the "location" is the argmax of probabilities (not differentiable) 
            # or we just use scalar field if it maps to labels.
            
            # Wait, if I use Categorical for EnMap, I can't easily concatenate with Normal.
            # But I can return a Product distribution? Pyro has dist.Independent or just sample twice.
            
            # If I sample twice, I need to pass y_obs_list as a tuple or dict.
            # Does GemPyPyroModel support that?
            
            return dist.Normal(
                torch.cat([torch.argmax(probs, dim=1).to(torch.float64), simulated_gravity]),
                torch.cat([torch.full((probs.shape[0],), 0.1, dtype=torch.float64), sigma_stations])
            ).to_event(1)

        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=set_priors_joint,
            likelihood_fn=joint_likelihood_fn, # I'll use a version that works with concatenation
            obs_name="Joint_Obs"
        )

        observed_joint = torch.cat([observed_labels, observed_gravity])

        return geo_model, observed_joint, prob_model, xy_ravel, norm_params
