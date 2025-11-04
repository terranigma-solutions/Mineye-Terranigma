import io
import os

import arviz as az
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist, optim
from pyro.infer import SVI, Trace_ELBO
from pyro.infer.autoguide import AutoNormalizingFlow, AutoIAFNormal
from pyro.optim import Adam

import gempy_probability as gpp
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize, create_orientation_modifier
from mineye.GeoModel.model_one.probabilistic_model_diagnostics import trace_pyro_model
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_diagonal
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel, read_gravity
from mineye.GeoModel.model_one.visualization import plot, gempy_viz
# noinspection PyUnusedImports
from tests import conftest
import numpy as np


class TestProbabilisticInversionVI:
    prior_key = r'dips'

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
        norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)

        # * 3) Setup Priors
        model_priors = {
                r'dips': dist.Normal(
                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # This is just dip 10 degrees
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                ).to_event(1)
        }

        # TODO: Here we could add density and range
        # * 4) Set up Deterministics
        pre_forward_dets = {
                "dips_degrees": lambda samples, gm: samples["dips"],  # Just pass through
        }

        post_forward_dets = {
                "gravity_response_raw": lambda samples, gm, sol: sol.gravity,  # Store raw gravity
                "gravity_response"    : lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),  # Normalized!
                "mean_gravity"        : lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
                "max_gravity"         : lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
        }

        # * 5) Set up likelihood functions
        likelihood_fn = generate_multigravity_likelihood_diagonal(
            norm_params=norm_params
        )

        # * 6) Set up Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=create_orientation_modifier(key=TestProbabilisticInversionVI.prior_key),
            likelihood_fn=likelihood_fn,
            pre_forward_deterministics=pre_forward_dets,
            post_forward_deterministics=post_forward_dets,
            obs_name="Gravity Measurement"
        )

        # Exploring model dependencies
        
        if False:
            self._graph(prob_model, geo_model, torch.tensor(observed_gravity_ugal))

        # * 7) Run predictive
        gravity_observations_tensor = torch.tensor(observed_gravity_ugal)
        trace = trace_pyro_model(prob_model, geo_model, torch.tensor(observed_gravity_ugal, dtype=torch.float64))
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=gravity_observations_tensor,
                n_samples=10,
                plot_trace=True
            )

        # * 8) Run Variational Inference with Normalizing Flows

        # Set up the guide (variational distribution) using normalizing flows
        guide = AutoIAFNormal(
            prob_model,
            [10],  # Hidden layer dimensions for the flow
            num_transforms=3
        )

        # Set up optimizer
        optimizer = optim.ClippedAdam({"lr": 1e-3, "clip_norm": 10.0})

        # Set up SVI
        svi = SVI(
            model=prob_model,
            guide=guide,
            optim=optimizer,
            loss=Trace_ELBO()
        )

        # Run SVI optimization
        num_iterations = 1000
        losses = []
        print("Starting VI optimization...")
        for i in range(num_iterations):
            loss = svi.step(geo_model, gravity_observations_tensor)
            losses.append(loss)
            if i % 500 == 0:
                print(f"Iteration {i}/{num_iterations}, Loss: {loss:.4f}")

        print(f"Final loss: {losses[-1]:.4f}")

        # Plot loss curve
        plt.figure(figsize=(10, 6))
        plt.plot(losses)
        plt.xlabel('Iteration')
        plt.ylabel('ELBO Loss')
        plt.title('Variational Inference Convergence')
        plt.grid(True)
        plt.show()

        # * 9) Sample from the variational posterior
        from pyro.infer import Predictive

        num_samples = 1000
        predictive = Predictive(
            model=prob_model,
            guide=guide,
            num_samples=num_samples
        )

        posterior_samples = predictive(geo_model, gravity_observations_tensor)

        # Convert to ArviZ InferenceData
        data = az.from_dict(
            posterior={k: v.detach().cpu().numpy()[None, :, :] if v.dim() > 1
            else v.detach().cpu().numpy()[None, :]
                       for k, v in posterior_samples.items()
                       if k in model_priors or k in pre_forward_dets or k in post_forward_dets},
            observed_data={"Gravity Measurement": observed_gravity_ugal}
        )

        if compute_prior_predictive:
            data.extend(prior_inference_data)

        data.to_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_vi.nc"))

        print("VI inference completed and saved!")

        return data

    def _graph(self, model, geo_model, y_obs_list=None):
        # ! This is not working well, the geo_model dependency is not properly picked up
        
        from pyro.infer.inspect import get_dependencies
        import pyro
        dependencies = get_dependencies(model, model_args=(geo_model, y_obs_list[:1]))
        dependencies

        # %%
        graph = pyro.render_model(
            model=model,
            model_args=(geo_model, y_obs_list,),
            render_params=True,
            render_distributions=True,
            render_deterministic=True
        )

        graph.attr(dpi='300')
        # Convert the graph to a PNG image format
        s = graph.pipe(format='png')

        # Open the image with PIL
        from PIL import Image
        image = Image.open(io.BytesIO(s))

        # Plot the image with matplotlib
        plt.figure(figsize=(10, 4))
        plt.imshow(image)
        plt.axis('off')  # Turn off axis
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
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_vi.nc"))

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
                hdi_prob=0.95,
                shade=.2,
                data_labels=["Posterior (VI)", "Prior"],
                colors=[default_red, default_blue],
            )

            # Apply log scale to all y-axes
            if isinstance(axes, np.ndarray):
                for ax in axes.flatten():
                    ax.set_yscale('log')
            else:
                axes.set_yscale('log')

            plt.show()

        # Analysis inference
        if hasattr(data, 'prior') and r'gravity_response' in data.prior:
            gravity_samples_norm, unit_label = plot(
                gravity_samples_norm=data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
                observed_gravity_ugal=observed_gravity_ugal,
                xy_ravel=xy_ravel
            )

        # Analysis Gempy Model
        gempy_viz(geo_model, data)
