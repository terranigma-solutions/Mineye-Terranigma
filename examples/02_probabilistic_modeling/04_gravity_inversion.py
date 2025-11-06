"""
Gravity Inversion with NUTS Sampling
=====================================

This example demonstrates Bayesian inversion using gravity data.
We use NUTS (No-U-Turn Sampler) to infer geological structure parameters from observed gravity.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import gempy as gp
import gempy_probability as gpp

import torch
import pyro
import pyro.distributions as dist

import arviz as az
from gempy_probability.core.samplers_data import NUTSConfig

# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(1234)

# %%
# Import helper functions
from mineye.config import paths
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel, read_gravity
from mineye.GeoModel.model_one.probabilistic_model import normalize, create_orientation_modifier
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_diagonal
from mineye.GeoModel.model_one.probabilistic_model_diagnostics import trace_pyro_model

# %%
# Define Model Extent and Resolution
# -----------------------------------

min_x = -709521
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
extent = [min_x, max_x, min_y, max_y, model_depth, max_z]

# Model resolution: use octree with refinement level 5
resolution = None
refinement = 5

# %%
# Get Data Paths
# --------------

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
geophysical_dir = paths.get_geophysical_dir()

print(f"Orientations: {mod_or_path}")
print(f"Points: {mod_pts_path}")
print(f"Geophysical data: {geophysical_dir}")

# %%
# Create Initial GemPy Model
# ---------------------------

simple_geo_model = gp.create_geomodel(
    project_name='gravity_inversion',
    extent=extent,
    refinement=refinement,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }
)

# %%
# Load Gravity Observations
# --------------------------

gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
print(f"\nGravity observations:")
print(f"  Number of measurements: {len(observed_gravity_ugal)}")
print(f"  Range: {observed_gravity_ugal.min():.1f} to {observed_gravity_ugal.max():.1f} µGal")
print(f"  Mean: {observed_gravity_ugal.mean():.1f} µGal")

# %%
# Setup Geomodel with Gravity Configuration
# ------------------------------------------

geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100

print(f"\n✓ Geomodel configured with {len(xy_ravel)} gravity measurement locations")

# %%
# Compute Baseline Forward Model
# -------------------------------
# Calculate forward gravity with initial parameters

baseline_fw_gravity_np = baseline(geo_model)
print(f"\nBaseline forward gravity:")
print(f"  Range: {baseline_fw_gravity_np.min():.1f} to {baseline_fw_gravity_np.max():.1f} µGal")
print(f"  Mean: {baseline_fw_gravity_np.mean():.1f} µGal")

# %%
# Normalize Forward Model to Observations
# ----------------------------------------

norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)
print(f"\nNormalization parameters:")
print(f"  Method: {norm_params['method']}")
print(f"  Parameters: {norm_params}")

# %%
# Define Prior Distribution
# --------------------------

n_orientations = geo_model.orientations_copy.xyz.shape[0]
model_priors = {
    r'dips': dist.Normal(
        loc=(torch.ones(n_orientations) * 10.0),
        scale=torch.tensor(10.0, dtype=torch.float64),
        validate_args=True
    )
}

print(f"\nPrior on dips:")
print(f"  Number of orientations: {n_orientations}")
print(f"  Mean: 10°")
print(f"  Std: 10°")

# %%
# Set Up Deterministics
# ----------------------
# Track intermediate values during inference

pre_forward_dets = {
    "dips_degrees": lambda samples, gm: samples["dips"],
}

post_forward_dets = {
    "gravity_response_raw": lambda samples, gm, sol: sol.gravity,
    "gravity_response": lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),
    "mean_gravity": lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
    "max_gravity": lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
}

# %%
# Set Up Likelihood Function
# ---------------------------

likelihood_fn = generate_multigravity_likelihood_diagonal(norm_params=norm_params)
print("✓ Likelihood function created")

# %%
# Create Probabilistic Model
# ---------------------------

prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
    priors=model_priors,
    set_interp_input_fn=create_orientation_modifier(key=r'dips'),
    likelihood_fn=likelihood_fn,
    pre_forward_deterministics=pre_forward_dets,
    post_forward_deterministics=post_forward_dets,
    obs_name="Gravity Measurement"
)

print("✓ Probabilistic model created")

# %%
# Trace Pyro Model
# ----------------
# Verify model structure

gravity_observations_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
trace = trace_pyro_model(prob_model, geo_model, gravity_observations_tensor)
print("✓ Model trace verified")

# %%
# Run Prior Predictive Check
# ---------------------------

print("\nRunning prior predictive sampling (10 samples)...")
prior_inference_data: az.InferenceData = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=gravity_observations_tensor,
    n_samples=10,
    plot_trace=True
)

print("✓ Prior predictive sampling complete")

# %%
# Run NUTS Inference
# ------------------
# Perform Bayesian inversion using NUTS sampler

print("\nRunning NUTS inference...")
print("  Warmup: 200 steps")
print("  Sampling: 200 samples")
print("  Chains: 2")

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
        num_samples=200,
        warmup_steps=200,
        num_chains=2
    ),
    plot_trace=True,
    run_posterior_predictive=True
)

print("✓ NUTS inference complete")

# %%
# Combine Prior and Posterior
# ----------------------------

data.extend(prior_inference_data)
print("✓ Prior and posterior combined")

# %%
# Summary Statistics
# ------------------

posterior_dips = data.posterior['dips'].values
print(f"\nPosterior dip statistics:")
print(f"  Shape: {posterior_dips.shape}")
print(f"  Mean: {posterior_dips.mean():.2f}°")
print(f"  Std: {posterior_dips.std():.2f}°")
print(f"  Median: {np.median(posterior_dips):.2f}°")

prior_dips = data.prior['dips'].values
print(f"\nPrior dip statistics:")
print(f"  Mean: {prior_dips.mean():.2f}°")
print(f"  Std: {prior_dips.std():.2f}°")

# %%
# Gravity Fit Statistics
# -----------------------

posterior_gravity = data.posterior_predictive['Gravity Measurement'].values
print(f"\nPosterior predictive gravity:")
print(f"  Shape: {posterior_gravity.shape}")
print(f"  Mean: {posterior_gravity.mean():.1f} µGal")
print(f"  Std: {posterior_gravity.std():.1f} µGal")

print(f"\nObserved gravity:")
print(f"  Mean: {observed_gravity_ugal.mean():.1f} µGal")
print(f"  Std: {observed_gravity_ugal.std():.1f} µGal")

residuals = posterior_gravity.mean(axis=(0, 1)) - observed_gravity_ugal
print(f"\nResiduals (posterior mean - observed):")
print(f"  Mean: {residuals.mean():.2f} µGal")
print(f"  RMS: {np.sqrt((residuals**2).mean()):.2f} µGal")
