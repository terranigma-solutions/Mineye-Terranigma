"""
Bayesian Joint Inversion: Hierarchical Gravity and EnMap
========================================================

This tutorial extends the joint inversion workflow by introducing a **hierarchical 
likelihood** for the gravity data.

**Why Hierarchical Joint Inversion?**

In real-world scenarios, different measurement stations might have varying levels of 
noise due to local conditions, instrument calibration, or terrain correction errors. 
A hierarchical model allows the inference process to:

1. **Estimate per-station noise**: Instead of assuming a single noise level for all 
   gravity stations, we infer the noise at each location.
2. **Robustness to outliers**: Stations with high noise are automatically downweighted,
   preventing them from biasing the geological model.
3. **Information sharing**: Through hyper-priors, stations "borrow strength" from 
   each other to better estimate their respective noise levels.

**Combining with EnMap**

While the gravity noise is modeled hierarchically, the EnMap data continues to provide 
categorical surface constraints, ensuring the resulting model is consistent with 
mapped lithology.
"""

import os
import torch
import numpy as np
import arviz as az
import matplotlib.pyplot as plt
import pyro
from pyro import distributions as dist
import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue

# Import helper functions from Mineye
from mineye.config import paths
from mineye.GeoModel.model_one.model_setup import setup_geomodel, read_gravity, baseline
from mineye.GeoModel.model_one.probabilistic_model import normalize, _modify_orientations, _modify_densities
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import _get_ordinal_probs
from mineye.GeoModel.model_one.visualization import (
    gempy_viz, 
    plot_joint_inversion_results
)

# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(1234)

# %%
# Define Parameter Mapping Function
# ---------------------------------
# We define how sampled parameters update the GemPy model.

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

# %%
# Load and Prepare Data
# ---------------------

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
geophysical_dir = paths.get_geophysical_dir()
base_dir = os.path.dirname(mod_or_path)

# 1. Load EnMap data
xyz_enmap = np.load(os.path.join(base_dir, 'central_xyz.npy'))
labels_enmap = np.load(os.path.join(base_dir, 'central_labels.npy'))
labels_enmap[labels_enmap == 2] = 1  # Normalize labels
observed_labels = torch.tensor(labels_enmap, dtype=torch.float64)

# 2. Load Gravity data
gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
observed_gravity = torch.tensor(observed_gravity_ugal, dtype=torch.float64)

# %%
# Setup Geomodel
# --------------

extent = [-707521, -675558, 4526832, 4551949, -500, 505]
simple_geo_model = gp.create_geomodel(
    project_name='joint_hierarchical',
    extent=extent,
    refinement=5,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)
gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={"Tournaisian_Plutonites": ["Tournaisian Plutonites"]}
)

# Configure grids
geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
gp.set_custom_grid(geo_model.grid, xyz_enmap)

# %%
# Normalization and Priors
# ------------------------

norm_params = normalize(
    baseline_fw_gravity_np=(baseline(geo_model)),
    observed_gravity=observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)

model_priors = {
    'dips': dist.Normal(
        loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),
        scale=torch.tensor(10, dtype=torch.float64),
        validate_args=True
    ).to_event(1),
    'density': dist.Normal(
        loc=(torch.tensor([2.9 - 2.67, 2.3 - 2.67])),
        scale=torch.tensor(0.15),
    ).to_event(1)
}

# %%
# Define Hierarchical Joint Likelihood
# ------------------------------------
#
# This is the core of the hierarchical approach. We sample `mu_log_sigma` and 
# `tau_log_sigma` as global hyperparameters, which then constrain the 
# `log_sigma_stations` for each individual measurement.

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
    # We invert the gravity sign if needed to match observations
    simulated_gravity = align_forward_to_observed(-solutions.gravity, norm_params)
    pyro.deterministic(r"$\mu_{gravity}$", simulated_gravity.detach())
    n_stations = simulated_gravity.shape[0]

    # Hyper-priors for log-noise
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
    
    # Per-station noise
    log_sigma_stations = pyro.sample(
        "log_sigma_stations",
        dist.Normal(mu_log_sigma.expand([n_stations]), tau_log_sigma).to_event(1)
    )
    sigma_stations = torch.exp(log_sigma_stations)
    pyro.deterministic("sigma_stations", sigma_stations)

    # --- C. Manual Sampling for Joint Data ---
    # Since we have mixed distribution types (Categorical and Normal), we call
    # pyro.sample directly for each component.
    pyro.sample("obs_enmap", dist.Categorical(probs=probs), obs=observed_labels)
    pyro.sample("obs_gravity", dist.Normal(simulated_gravity, sigma_stations).to_event(1), obs=observed_gravity)
    
    return [] # Return an empty list to satisfy make_gempy_pyro_model's loop without sampling more

# %%
# Create Probabilistic Model
# --------------------------

prob_model = gpp.make_gempy_pyro_model(
    priors=model_priors,
    set_interp_input_fn=set_priors_joint,
    likelihood_fn=joint_likelihood_fn,
    obs_name="Joint_Obs"
)

# For hierarchical models where we sample multiple observations in the likelihood,
# we can pass a combined list of tensors to match the empty list if needed,
# but since we sample manually, y_obs_list is just to satisfy the API.
observed_joint = [observed_labels, observed_gravity]

# %%
# Run Inference
# -------------

print("Running joint hierarchical NUTS inference...")
data = gpp.run_nuts_inference(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=observed_joint,
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

# %%
# Results and Outlier Detection
# -----------------------------

# Visualize joint results
gravity_xyz = np.zeros((len(gravity_data), 3))
gravity_xyz[:, 0] = gravity_data['X'].values
gravity_xyz[:, 1] = gravity_data['Y'].values
gravity_xyz[:, 2] = gravity_data['Z'].values

plot_joint_inversion_results(
    data=data,
    observed_gravity=observed_gravity_ugal,
    xy_gravity=gravity_xyz[:, :2],
    observed_enmap=labels_enmap,
    geo_model=geo_model,
    n_gravity_points=len(gravity_data)
)

# Outlier Detection
if "sigma_stations" in data.posterior_predictive:
    posterior_sigmas = data.posterior_predictive["sigma_stations"].values  
    station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
    sigma_global_mean = station_noise_mean.mean()
    problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
    print(f"Potential outlier stations (noise > 2x mean): {problematic}")

plt.show()
