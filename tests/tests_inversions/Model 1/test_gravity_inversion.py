import os
from functools import partial

import arviz as az
import geopandas as gpd
import numpy as np
import pyro
import torch

import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
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
        likelihood_fn = generate_multigravity_likelihood(
            covariance_matrix=(gaussian_kernel(xy_ravel[:, :2], length_scale, variance)),
            norm_params=norm_params
        )
        
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
        gravity_observations_tensor = torch.tensor(observed_gravity_ugal)
        trace = trace_pyro_model(prob_model, geo_model, torch.tensor(observed_gravity_ugal, dtype=torch.float64))
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=gravity_observations_tensor,
                n_samples=n_samples,
                plot_trace=True
            )

        # * 8) Run inference

        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=gravity_observations_tensor,
            config=NUTSConfig(
                step_size=0.1,
                adapt_step_size=True,
                target_accept_prob=0.9,
                max_tree_depth=10,
                init_strategy='auto',
                num_samples=200,
                warmup_steps=50,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )
        
        if compute_prior_predictive:
            data.extend(prior_inference_data)

        # After MCMC
        print(f"Divergences: {data.sample_stats.diverging.sum().item()}")
        print(f"Max tree depth: {(data.sample_stats.tree_depth == 10).sum().item()}")
        print(f"ESS: {az.ess(data)}")
        print(f"R-hat: {az.rhat(data)}")  # Should be < 1.01
        # 
        # # Posterior predictive checks
        az.plot_ppc(data, num_pp_samples=100)
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


def generate_multigravity_likelihood(covariance_matrix, norm_params):
    return partial(multigravity_likelihood, covariance_matrix=covariance_matrix, norm_params=norm_params)


def multigravity_likelihood(solutions: gp.data.Solutions, covariance_matrix, norm_params) -> dist:
    simulated_geophysics = align_forward_to_observed(-solutions.gravity, norm_params)
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

def debug_gradient_flow(likelihood_fn, solutions, obs_data):
    """
    Systematically check where gradient flow breaks.
    
    Call this BEFORE running NUTS to diagnose issues.
    """
    import torch

    print("\n" + "="*80)
    print("GRADIENT FLOW DIAGNOSTIC")
    print("="*80)

    # Get the distribution
    lik_dist = likelihood_fn(solutions)

    # Compute log probability
    try:
        log_prob = lik_dist.log_prob(obs_data)
        print(f"✓ Log probability computed: {log_prob.item():.4f}")
    except Exception as e:
        print(f"✗ Failed to compute log_prob: {e}")
        return

    # Check if we can compute gradients
    print(f"\nLog prob requires_grad: {log_prob.requires_grad}")
    print(f"Log prob grad_fn: {log_prob.grad_fn}")

    if not log_prob.requires_grad:
        print("\n⚠️  WARNING: log_prob does not require gradients!")
        print("This means NUTS cannot compute gradients for sampling.")

        # Trace back through the computation graph
        print("\nTracing back through computation:")
        print(f"  Distribution type: {type(lik_dist).__name__}")

        # Check distribution parameters
        if hasattr(lik_dist, 'loc'):
            print(f"  loc requires_grad: {lik_dist.loc.requires_grad}")
            print(f"  loc grad_fn: {lik_dist.loc.grad_fn}")

        if hasattr(lik_dist, 'scale_tril'):
            print(f"  scale_tril requires_grad: {lik_dist.scale_tril.requires_grad}")
            print(f"  scale_tril grad_fn: {lik_dist.scale_tril.grad_fn}")

        if hasattr(lik_dist, 'covariance_matrix'):
            print(f"  covariance_matrix requires_grad: {lik_dist.covariance_matrix.requires_grad}")
            print(f"  covariance_matrix grad_fn: {lik_dist.covariance_matrix.grad_fn}")

        if hasattr(lik_dist, 'df'):
            print(f"  df requires_grad: {lik_dist.df.requires_grad}")
            print(f"  df grad_fn: {lik_dist.df.grad_fn}")
    else:
        print("\n✓ Gradient flow is intact!")

        # Try computing gradients
        try:
            # Create dummy parameters to test backward pass
            test_grad = torch.autograd.grad(
                outputs=log_prob,
                inputs=[p for p in lik_dist.parameters() if p.requires_grad],
                allow_unused=True
            )
            print(f"✓ Backward pass successful, got {len(test_grad)} gradients")
        except Exception as e:
            print(f"✗ Backward pass failed: {e}")

    print("="*80 + "\n")


# Usage in your test:
def test_gravity_inversion(self, simple_geo_model, geophysical_dir, n_samples=50):
    # ... setup code ...

    # Before running NUTS, test gradient flow:
    from gempy_probability.modules.forwards import run_gempy_forward
    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame

    # Create a test forward pass
    test_interp_input = interpolation_input_from_structural_frame(geo_model)
    test_solutions = run_gempy_forward(test_interp_input, geo_model)

    # Debug the likelihood
    debug_gradient_flow(likelihood_fn, test_solutions, torch.tensor(observed_gravity_ugal, dtype=torch.float64))

    # ... continue with NUTS ...


def trace_pyro_model(prob_model, geo_model, obs_data, print_full=True):
    """
    Use Pyro's trace functionality to see all sample sites and deterministics.

    This shows you EXACTLY what Pyro sees.
    """
    from pyro import poutine
    import torch

    print("\n" + "=" * 80)
    print("PYRO MODEL TRACE")
    print("=" * 80 + "\n")

    # Trace the model execution
    trace = poutine.trace(prob_model).get_trace(geo_model, obs_data)

    print(f"{'Site Name':<30} | {'Type':<15} | {'Shape':<15} | {'Grad?':<10} | {'Grad Fn'}")
    print("-" * 100)

    for name, node in trace.nodes.items():
        if node['type'] == 'sample':
            value = node['value']

            if isinstance(value, torch.Tensor):
                shape_str = str(tuple(value.shape))
                grad_str = "✓" if value.requires_grad else "✗"
                grad_fn_str = str(value.grad_fn)[:30] if value.grad_fn else "None"

                print(f"{name:<30} | {node['type']:<15} | {shape_str:<15} | {grad_str:<10} | {grad_fn_str}")

                # Flag problematic ones
                if not value.requires_grad and name not in ['obs', 'Gravity Measurement']:
                    print(f"  ⚠️  WARNING: Sample site '{name}' has no gradient flow!")
            else:
                print(f"{name:<30} | {node['type']:<15} | {str(type(value)):<15} | {'N/A':<10} | N/A")

    if print_full:
        print("\n" + "=" * 80)
        print("FULL TRACE DETAILS")
        print("=" * 80 + "\n")
        print(trace.format_shapes())

    print("\n" + "=" * 80 + "\n")

    return trace