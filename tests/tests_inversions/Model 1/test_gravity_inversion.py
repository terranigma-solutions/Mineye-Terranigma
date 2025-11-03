import os
from functools import partial

import arviz as az
import geopandas as gpd
import pyro
import torch
from pyro import distributions as dist

import gempy_probability as gpp
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_1_aux import setup_geomodel, baseline, normalize, plot, gempy_viz, create_orientation_modifier
# noinspection PyUnusedImports
from tests import conftest


class TestProbabilisticInversion:
    prior_key = r'dips'

    def test_gravity_inversion(self, simple_geo_model, geophysical_dir, n_samples=50):
        """Test reading and computing a geological model."""

        # Use actual gravity measurement device locations
        # * 1) Read gravity data
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal
        observed_gravity_ugal = observed_gravity * 1000

        # * 2) Setup initial Geomodel and normalize forward gravity to the observed gravity
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        baseline_fw_gravity_np = baseline(geo_model)
        norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)

        # * 3) Setup Priors
        model_priors = {
                r'dips': dist.Normal(
                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # This is just dip 10 degrees
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                )
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
        length_scale_prior = torch.tensor(1_000.0)
        variance_prior = torch.tensor(25.0 ** 2)
        covariance_matrix = gaussian_kernel(xy_ravel[:,:2], length_scale_prior, variance_prior)
        likelihood_fn = generate_multigravity_likelihood(covariance_matrix)
        
        # * 6) Set up Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=create_orientation_modifier(key=TestProbabilisticInversion.prior_key),
            likelihood_fn=likelihood_fn,
            pre_forward_deterministics=pre_forward_dets,
            post_forward_deterministics=post_forward_dets,
            obs_name="Gravity Measurement"
        )

        # * 7) Run predictive
        prior_inference_data: az.InferenceData = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=observed_gravity_ugal,
            n_samples=n_samples,
            plot_trace=True
        )

        # * 8) Run inference

        # * 9) Analysis inference
        gravity_samples_norm, unit_label = plot(
            gravity_samples_norm=prior_inference_data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
            observed_gravity_ugal=observed_gravity_ugal,
            xy_ravel=xy_ravel
        )

        # * 9) Analysis Gempy Model

        gempy_viz(geo_model, prior_inference_data)


import gempy as gp
import pyro.distributions as dist


def generate_multigravity_likelihood(covariance_matrix):
    return partial(multigravity_likelihood, covariance_matrix=covariance_matrix)


def multigravity_likelihood(solutions: gp.data.Solutions, covariance_matrix) -> dist:
    simulated_geophysics = solutions.gravity
    pyro.deterministic(r'$\mu_{gravity}$', simulated_geophysics)
    normal = dist.MultivariateNormal(simulated_geophysics, covariance_matrix)
    return normal


def gaussian_kernel(locations, length_scale, variance):
    import torch
    # Compute the squared Euclidean distance between each pair of points
    locations = torch.tensor(locations)
    distance_squared = torch.cdist(locations, locations, p=2).pow(2)
    # Compute the covariance matrix using the Gaussian kernel
    covariance_matrix = variance * torch.exp(-0.5 * distance_squared / length_scale ** 2)
    return covariance_matrix
