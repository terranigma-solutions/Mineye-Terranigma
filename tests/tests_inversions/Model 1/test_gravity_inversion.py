import os

import arviz as az
import numpy as np
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist

import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality
from mineye.GeoModel.model_one.probabilistic_model import normalize, create_orientation_modifier, set_priors
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_diagonal, generate_multigravity_likelihood_hierarchical_per_station
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel, read_gravity
from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots, gempy_viz, plot_many_observed_vs_forward
from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_comparison
# noinspection PyUnusedImports
from tests import conftest


class TestProbabilisticInversion:
    prior_key_dips = r'dips'
    prior_key_density = r'density'

    def test_run_predictive_foo(self, simple_geo_model, geophysical_dir, n_samples=50):
        geo_model, observed_gravity_ugal, prob_model = self._create_probabilistic_model(geophysical_dir, simple_geo_model)
        prob_model(geo_model, observed_gravity_ugal)
        pass

    def test_run_predictive(self, simple_geo_model, geophysical_dir, n_samples=50):
        geo_model, observed_gravity_ugal, prob_model = self._create_probabilistic_model(geophysical_dir, simple_geo_model)

        # * 7) Run predictive
        gravity_observations_tensor = torch.tensor(observed_gravity_ugal)
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=gravity_observations_tensor,
                n_samples=100,
                plot_trace=True
            )

    def test_gravity_inversion(self, simple_geo_model, geophysical_dir, n_samples=50,
                               arviz_data_filename="arviz_data_Nov10_I_hierarchical.nc"):
        """Test reading and computing a geological model."""
        print("Test gravity inversion...")
        # Use actual gravity measurement device locations
        geo_model, observed_gravity_ugal, prob_model = self._create_probabilistic_model(geophysical_dir, simple_geo_model)

        # * 7) Run predictive
        gravity_observations_tensor = torch.tensor(observed_gravity_ugal)
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=gravity_observations_tensor,
                n_samples=100,
                plot_trace=True
            )

        # * 8) Run inference

        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=gravity_observations_tensor,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                # num_samples=20,
                # warmup_steps=5,
                num_samples=200,
                warmup_steps=200,
                num_chains=1
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        if compute_prior_predictive:
            data.extend(prior_inference_data)

        data.to_netcdf(os.path.join(os.path.dirname(__file__), arviz_data_filename))

    @staticmethod
    def _create_probabilistic_model(geophysical_dir, simple_geo_model):
        # * 1) Read gravity data
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)

        # * 2) Setup initial Geomodel and normalize forward gravity to the observed gravity
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        geo_model.interpolation_options.sigmoid_slope = 100
        norm_params = normalize(
            baseline_fw_gravity_np=(baseline(geo_model)),
            observed_gravity=observed_gravity_ugal,
            method="align_to_reference",
            extrapolation_buffer=0.3
        )

        # * 3) Setup Priors
        model_priors = {
                TestProbabilisticInversion.prior_key_dips   : dist.Normal(
                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # This is just dip 10 degrees
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                ),
                TestProbabilisticInversion.prior_key_density: dist.Normal(
                    loc=(torch.tensor([
                            2.9,  # plutonites
                            2.3  # host
                    ])),
                    scale=torch.tensor(0.15),
                ).to_event(1)
        }

        # TODO: Here we could add density and range
        # * 4) Set up Deterministics
        post_forward_dets = {
                "gravity_response_raw": lambda samples, gm, sol: sol.gravity,  # Store raw gravity
                "gravity_response"    : lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),  # Normalized!
                "mean_gravity"        : lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
                "max_gravity"         : lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
        }

        # * 5) Set up likelihood functions

        likelihood_fn = generate_multigravity_likelihood_hierarchical_per_station(
            norm_params=norm_params
        )

        # * 6) Set up Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=set_priors,
            likelihood_fn=likelihood_fn,
            pre_forward_deterministics={},
            post_forward_deterministics=post_forward_dets,
            obs_name="Gravity Measurement"
        )
        return geo_model, observed_gravity_ugal, prob_model

    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data.nc"))
        # Run comprehensive diagnostics with plots
        check_mcmc_quality(data, verbose=True, plot=True)

    def test_run_predictive_analysis(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_Nov10_I_hierarchical.nc"))
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

        # Prepare data
        observed_norm = observed_gravity_ugal
        forward_norm = data.prior[r'gravity_response'].mean(axis=1)
        if PRIOR := True:
            many_forward_norm = data.prior[r'gravity_response'].values[0, -20:]
        else:
            many_forward_norm = data.posterior_predictive[r'gravity_response'].values[0, -40:-20]

        plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)
        
    def test_run_kde_sections(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_Nov10_I_hierarchical.nc"))
        
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        
        gempy_viz(geo_model, data, n_samples=100)

    def test_run_outlier_detection(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_Nov10_I_hierarchical.nc"))

        posterior_sigmas = data.posterior_predictive["sigma_stations"].values  # shape: (chains, samples, 20)

        # Stations with consistently high σ are either:
        # - Actually noisier (geology, instrument, location)
        # - Outliers / model misspecification
        # - Areas where your forward model is wrong

        station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
        sigma_global_mean = station_noise_mean.mean()
        problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
        print(f"Check stations: {problematic}")
        axes = az.plot_density(
            data=[data, data.prior],
            var_names=["sigma_stations"],
            filter_vars="like",
            hdi_prob=0.9999,
            shade=.2,
            data_labels=["Posterior", "Prior"],
            colors=[default_red, default_blue],
        )

        plt.show()

    def test_run_analysis(self, simple_geo_model, geophysical_dir):

        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_Nov10_I_hierarchical.nc"))
        data.posterior

        az.plot_posterior(data, var_names=["dips", "density"])
        plt.show()

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

        # # Posterior predictive checks
        # az.plot_ppc(data, num_pp_samples=20)

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

        # # Apply log scale to all y-axes
        # if isinstance(axes, np.ndarray):
        #     for ax in axes.flatten():
        #         ax.set_yscale('log')
        # else:
        #     axes.set_yscale('log')

        plt.show()
        plot_gravity_comparison(
            observed_ugal=observed_gravity_ugal,
            forward_norm=data.prior[r'gravity_response'].mean(axis=1),
            xy_ravel=xy_ravel,
            normalization_method='align_to_reference'
        )

        plot_gravity_comparison(
            observed_ugal=observed_gravity_ugal,
            forward_norm=data.posterior_predictive[r'gravity_response'].mean(axis=1),
            xy_ravel=xy_ravel,
            normalization_method='align_to_reference'
        )

        # * 9) Analysis inference
        if hasattr(data, 'prior') and r'gravity_response' in data.prior:
            gravity_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        if hasattr(data, 'posterior') and r'gravity_response' in data.prior:
            gravity_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.posterior_predictive[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        # * 9) Analysis Gempy Model
        gempy_viz(geo_model, data)
