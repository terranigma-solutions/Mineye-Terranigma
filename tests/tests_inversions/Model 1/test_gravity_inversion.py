import os
from functools import partial

import arviz as az
import geopandas as gpd
import numpy as np
import pyro
import torch

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
        # length_scale_prior = torch.tensor(1_000.0)
        # variance_prior = torch.tensor(25.0 ** 2)

        length_scale = pyro.sample(
            "length_scale",
            dist.LogNormal(
                loc=torch.log(torch.tensor(2000.0)),  # Median = 2 km
                scale=1.0  # 68% interval: [~700m, ~5.5km], 95%: [~250m, ~16km]
            )
        )
        # Option A: Inverse-Gamma (conjugate, traditional)
        variance = pyro.sample(
            "variance",
            dist.InverseGamma(
                concentration=3.0,  # Shape parameter
                rate=75000.0        # Scale parameter
            )
            # This gives: Mean ≈ 37,500, Mode ≈ 25,000
            # Roughly covers your 15k-40k range
        )
        covariance_matrix = gaussian_kernel(xy_ravel[:,:2], length_scale, variance)
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
        # After MCMC
        # print(f"Divergences: {data.sample_stats.diverging.sum().item()}")
        # print(f"Max tree depth: {(data.sample_stats.tree_depth == 10).sum().item()}")
        # print(f"ESS: {az.ess(data)}")
        # print(f"R-hat: {az.rhat(data)}")  # Should be < 1.01
        # 
        # # Posterior predictive checks
        # az.plot_ppc(data, num_pp_samples=100)
        # * 9) Analysis inference
        gravity_samples_norm, unit_label = plot(
            gravity_samples_norm=prior_inference_data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
            observed_gravity_ugal=observed_gravity_ugal,
            xy_ravel=xy_ravel
        )

        # * 9) Analysis Gempy Model

        gempy_viz(geo_model, prior_inference_data)


    def test_gravity_duplicates(self, geophysical_dir):
        """Test reading and computing a geological model."""

        # Use actual gravity measurement device locations
        # * 1) Read gravity data
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

        # * 1b) CHECK FOR DUPLICATES
        xy_coords = gravity_data[['geometry']].apply(lambda row: (row.geometry.x, row.geometry.y), axis=1)
        xy_array = np.array(xy_coords.tolist())

        # Find unique locations (within 1m tolerance)
        from scipy.spatial.distance import pdist, squareform
        distances = squareform(pdist(xy_array))
        np.fill_diagonal(distances, np.inf)  # Ignore self-distances

        duplicates = np.any(distances < 1.0, axis=1)  # Locations within 1m
        if np.any(duplicates):
            print(f"WARNING: Found {duplicates.sum()} stations with duplicates/near-duplicates")
            print("Removing duplicates...")

            # Keep first occurrence of each duplicate set
            keep_mask = np.ones(len(xy_array), dtype=bool)
            for i in np.where(duplicates)[0]:
                if keep_mask[i]:
                    # Find all duplicates of this point
                    dups = distances[i] < 1.0
                    # Keep first, remove rest
                    dup_indices = np.where(dups)[0]
                    if len(dup_indices) > 0:
                        keep_mask[dup_indices[1:]] = False

            # Filter data
            gravity_data = gravity_data[keep_mask]
            observed_gravity = gravity_data['VALU_BOU267'].values
            observed_gravity_ugal = observed_gravity * 1000
            print(f"Kept {keep_mask.sum()} unique stations out of {len(keep_mask)}")



import gempy as gp
import pyro.distributions as dist


def generate_multigravity_likelihood(covariance_matrix):
    return partial(multigravity_likelihood, covariance_matrix=covariance_matrix)


def multigravity_likelihood(solutions: gp.data.Solutions, covariance_matrix) -> dist:
    simulated_geophysics = solutions.gravity
    pyro.deterministic(r'$\mu_{gravity}$', simulated_geophysics)

    nu = pyro.sample(
        name="nu",
        fn=(dist.Exponential(0.2))  # Mean = 5, favors moderate heavy tails
    ) + 2.0  # Ensures nu > 2 (so variance exists)
    # Student-t for heavy tails
    likelihood = dist.MultivariateStudentT(
        df=nu,  # degrees of freedom (sample as hyperparameter)
        loc=simulated_geophysics,
        scale_tril=torch.linalg.cholesky(covariance_matrix)
    )
    return likelihood


def gaussian_kernel_(locations, length_scale, variance):
    import torch
    # Compute the squared Euclidean distance between each pair of points
    locations = torch.tensor(locations)
    distance_squared = torch.cdist(locations, locations, p=2).pow(2)
    # Compute the covariance matrix using the Gaussian kernel
    covariance_matrix = variance * torch.exp(-0.5 * distance_squared / length_scale ** 2)
    return covariance_matrix


def gaussian_kernel(locations, length_scale, variance, nugget=None):
    """
    Numerically stable Gaussian kernel with automatic jitter.

    Parameters
    ----------
    locations : torch.Tensor or np.ndarray
        Shape (n_stations, 2) - station coordinates
    length_scale : float or torch.Tensor
        Correlation length scale (same units as locations)
    variance : float or torch.Tensor
        Signal variance (µGal²)
    nugget : float or torch.Tensor, optional
        Nugget effect (independent measurement noise).
        If None, uses 0.1% of variance as minimum

    Returns
    -------
    covariance_matrix : torch.Tensor
        Positive-definite (n_stations, n_stations) covariance
    """
    import torch

    # Type safety
    locations = torch.tensor(locations, dtype=torch.float64)
    length_scale = torch.as_tensor(length_scale, dtype=torch.float64)
    variance = torch.as_tensor(variance, dtype=torch.float64)

    # Default nugget: 0.1% of signal variance (common practice)
    if nugget is None:
        nugget = 0.001 * variance
    else:
        nugget = torch.as_tensor(nugget, dtype=torch.float64)

    n_stations = locations.shape[0]

    # Compute distances
    distance_squared = torch.cdist(locations, locations, p=2).pow(2)

    # Stabilized exponential: avoid underflow for large distances
    # exp(-x) ≈ 0 for x > 30, so clip argument
    exponent = -0.5 * distance_squared / (length_scale ** 2 + 1e-10)  # Avoid division by zero
    exponent = torch.clamp(exponent, min=-30.0)  # Prevent underflow

    # Kernel
    K = variance * torch.exp(exponent)

    # Add nugget (diagonal noise)
    K = K + torch.eye(n_stations, dtype=torch.float64) * nugget

    # Final safety check - if still not PD, add more jitter
    try:
        torch.linalg.cholesky(K)
    except RuntimeError:
        # Gradually increase jitter until PD
        extra_jitter = 1e-6 * variance
        max_attempts = 10
        for attempt in range(max_attempts):
            K = K + torch.eye(n_stations, dtype=torch.float64) * extra_jitter
            try:
                torch.linalg.cholesky(K)
                print(f"Warning: Added extra jitter {extra_jitter:.2e} to ensure positive-definiteness")
                break
            except RuntimeError:
                extra_jitter *= 10
        else:
            raise RuntimeError("Could not make covariance matrix positive-definite even with jitter")

    return K

def gaussian_kernel__(locations, length_scale, variance):
    import torch

    # Ensure consistent dtype
    locations = torch.tensor(locations, dtype=torch.float64)
    length_scale = torch.as_tensor(length_scale, dtype=torch.float64)
    variance = torch.as_tensor(variance, dtype=torch.float64)

    # Compute squared Euclidean distances
    distance_squared = torch.cdist(locations, locations, p=2).pow(2)

    # Gaussian kernel
    covariance_matrix = variance * torch.exp(-0.5 * distance_squared / length_scale ** 2)

    # DIAGNOSTIC: Print matrix properties
    print(f"Covariance matrix shape: {covariance_matrix.shape}")
    print(f"Diagonal min/max: {covariance_matrix.diag().min():.2e} / {covariance_matrix.diag().max():.2e}")
    print(f"Off-diagonal min/max: {covariance_matrix[~torch.eye(covariance_matrix.shape[0], dtype=bool)].min():.2e} / {covariance_matrix[~torch.eye(covariance_matrix.shape[0], dtype=bool)].max():.2e}")
    print(f"Condition number: {torch.linalg.cond(covariance_matrix):.2e}")
    print(f"Min eigenvalue: {torch.linalg.eigvalsh(covariance_matrix)[0]:.2e}")

    # Check for duplicates
    unique_locs = torch.unique(locations, dim=0)
    if unique_locs.shape[0] < locations.shape[0]:
        print(f"WARNING: Found {locations.shape[0] - unique_locs.shape[0]} duplicate locations!")

    return covariance_matrix