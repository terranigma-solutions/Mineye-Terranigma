import os

import arviz as az
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist, optim
from pyro.distributions import Distribution
from pyro.infer import SVI, Trace_ELBO
from pyro.infer.autoguide import AutoIAFNormal
from pyro.optim import Adam

import gempy_probability as gpp
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize, create_orientation_modifier, _modify_orientations, _modify_densities, set_priors
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_diagonal, generate_multigravity_likelihood_hierarchical_per_station, generate_multigravity_likelihood_per_station_stable
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel, read_gravity
from mineye.GeoModel.model_one.visualization import plot, gempy_viz
from mineye.GeoModel.plotting.probabilistic_analysis import _plot_comparison, plot_gravity_comparison
# noinspection PyUnusedImports
from tests import conftest
import numpy as np

import gempy as gp


class TestProbabilisticInversionVI:
    prior_key_dips = r'dips'
    prior_key_density = r'density'

    def test_gravity_inversion_vi(self, simple_geo_model, geophysical_dir, n_samples=50):
        """Test reading and computing a geological model with Variational Inference using Normalizing Flows."""
        print("Test gravity inversion with VI...")
        # Use actual gravity measurement device locations
        # * 1) Read gravity data
        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)

        # * 2) Setup initial Geomodel and normalize forward gravity to the observed gravity
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        geo_model.interpolation_options.sigmoid_slope = 100
        baseline_fw_gravity_np = baseline(geo_model)
        norm_params = normalize(
            baseline_fw_gravity_np=(baseline(geo_model)),
            observed_gravity=observed_gravity_ugal,
            method="align_to_reference",
            extrapolation_buffer=0.3
        )

        # * 3) Setup Priors
        model_priors = {
                TestProbabilisticInversionVI.prior_key_dips   : dist.Normal(

                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # This is just dip 10 degrees
                    scale=torch.tensor(20, dtype=torch.float64),
                    validate_args=True
                ).to_event(1),
                TestProbabilisticInversionVI.prior_key_density: dist.Normal(
                    loc=(torch.tensor([
                            2.9,  # plutonites
                            2.3  # host
                    ])),
                    scale=torch.tensor(0.3),
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
        # likelihood_fn = generate_multigravity_likelihood_hierarchical_per_station(
        #     norm_params=norm_params
        # )

        likelihood_fn = generate_multigravity_likelihood_per_station_stable(
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

        # Exploring model dependencies
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

            if False:
                self._plot_prior_predictive(geo_model, observed_gravity_ugal, prior_inference_data, xy_ravel)
                return 

        # * 8) Run Variational Inference with Normalizing Flows

        guide, losses = self._run_vi(geo_model, gravity_observations_tensor, prob_model)

        self._plot_loss_curve(losses)
        # * 9) Sample from the variational posterior
        data = self._get_posterior_predictive(
            geo_model,
            gravity_observations_tensor,
            guide,
            model_priors,
            observed_gravity_ugal,
            post_forward_dets,
            pre_forward_dets={},
            prob_model=prob_model
        )

        if compute_prior_predictive:
            data.extend(prior_inference_data)

        data.to_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_vi_Nov10_I_hierarchical.nc"))

        print("VI inference completed and saved!")

        return data

    def _run_vi(self, geo_model, gravity_observations_tensor, prob_model):
        # Set up the guide (variational distribution) using normalizing flows
        guide = AutoIAFNormal(
            prob_model,
            [64],  # Hidden layer dimensions for the flow
            num_transforms=3
        )

        # Set up optimizer
        # optimizer = optim.ClippedAdam({"lr": 1e-4, "clip_norm": 10.0})
        optimizer = optim.ClippedAdam({
                "lr"       : 5e-3,  # Start higher, then anneal
                "clip_norm": 1.0,  # Tighter clipping
                "betas"    : (0.9, 0.999),
                "lrd"      : 0.9996  # Learning rate decay
        })

        # Set up SVI
        svi = SVI(
            model=prob_model,
            guide=guide,
            optim=optimizer,
            loss=Trace_ELBO(num_particles=5)
        )

        # Run SVI optimization
        num_iterations = 5000
        losses = []
        print("Starting VI optimization...")
        for i in range(num_iterations):
            loss = svi.step(geo_model, gravity_observations_tensor)
            if torch.isnan(torch.tensor(loss)) or torch.isinf(torch.tensor(loss)):
                print(f"NaN/Inf at iteration {i}")
                break
            losses.append(loss)
            if i % 500 == 0:
                print(f"Iteration {i}/{num_iterations}, Loss: {loss:.4f}")

        print(f"Final loss: {losses[-1]:.4f}")
        return guide, losses

    @staticmethod
    def _get_posterior_predictive(geo_model, gravity_observations_tensor, guide, model_priors,
                                  observed_gravity_ugal, post_forward_dets,
                                  pre_forward_dets, prob_model):
        from pyro.infer import Predictive

        num_samples = 1000
        predictive = Predictive(
            model=prob_model,
            guide=guide,
            num_samples=num_samples
        )

        posterior_samples = predictive(geo_model, gravity_observations_tensor)

        # Convert to ArviZ InferenceData
        posterior_dict = {}
        for k, v in posterior_samples.items():
            if k in model_priors or k in pre_forward_dets or k in post_forward_dets:
                if v.dim() > 1:
                    posterior_dict[k] = v.detach().cpu().numpy()[None, :, :]
                else:
                    posterior_dict[k] = v.detach().cpu().numpy()[None, :]

        data = az.from_dict(
            posterior=posterior_dict,
            observed_data={"Gravity Measurement": observed_gravity_ugal}
        )
        return data

    @staticmethod
    def _plot_prior_predictive(geo_model, observed_gravity_ugal, data, xy_ravel):
        gravity_samples_norm = data.prior[r'gravity_response'].values[0, :]  # (n_samples, n_devices)

        plot_gravity_comparison(
            observed_ugal=observed_gravity_ugal,
            forward_norm=data.prior[r'gravity_response'].mean(axis=1),
            xy_ravel=xy_ravel,
            normalization_method='align_to_reference'
        )

        gravity_samples_norm, unit_label = plot(
            gravity_samples_norm=gravity_samples_norm,
            observed_gravity_ugal=observed_gravity_ugal,
            xy_ravel=xy_ravel
        )
        # gempy_viz(geo_model, prior_inference_data)

        baseline_fw_gravity_np = baseline(geo_model)

    @staticmethod
    def _plot_loss_curve(losses):
        plt.figure(figsize=(10, 6))

        # Filter out extreme values for better visualization
        loss_array = np.array(losses)
        upper_bound = np.percentile(loss_array, 99)
        filtered_losses = np.where(loss_array < upper_bound, loss_array, upper_bound)

        plt.plot(filtered_losses)
        plt.xlabel('Iteration')
        plt.ylabel('ELBO Loss (filtered)')
        plt.title('Variational Inference Convergence')
        plt.grid(True)
        plt.show()

    def test_run_diagnostics_vi(self):
        """Run diagnostics on VI results."""
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_vi.nc"))

        # Note: Traditional MCMC diagnostics (like R-hat, ESS) are not applicable to VI
        # Instead, we can look at posterior distributions and posterior predictive checks
        print("VI posterior summary:")
        print(az.summary(data))

    def test_run_analysis_vi(self, simple_geo_model, geophysical_dir):
        """Analyze VI results."""
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_vi_Nov10_I_hierarchical.nc"))

        gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

        # Posterior predictive checks
        if hasattr(data, 'posterior_predictive'):
            az.plot_ppc(data, num_pp_samples=20)
            plt.show()

        plt.rcParams['figure.dpi'] = 72  # Lower DPI for faster rendering

        # Compare prior and posterior
        if hasattr(data, 'prior'):
            axes = az.plot_density(
                data=[data, data.prior],
                var_names=["dips"],
                filter_vars="like",
                hdi_prob=0.999,
                shade=.2,
                data_labels=["Posterior (VI)", "Prior"],
                colors=[default_red, default_blue],
            )

            # Apply log scale to all y-axes
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
            forward_norm=data.posterior[r'gravity_response'].mean(axis=1),
            xy_ravel=xy_ravel,
            normalization_method='align_to_reference'
        )
        
        # Analysis inference
        if hasattr(data, 'prior') and r'gravity_response' in data.prior:
            gravity_samples_norm, unit_label = plot(
                gravity_samples_norm=data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        if hasattr(data, 'posterior') and r'gravity_response' in data.prior:
            gravity_samples_norm, unit_label = plot(
                gravity_samples_norm=data.posterior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        if False:
            # Analysis Gempy Model
            gempy_viz(geo_model, data)


