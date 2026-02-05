"""
Bayesian Joint Inversion: Gravity and EnMap Data
================================================

This tutorial demonstrates how to perform a Bayesian joint inversion using GemPy-Probability.
We combine two different types of geophysical and surface data:

1. **Gravity Data**: Subsurface density constraints from gravitational measurements.
2. **EnMap Data**: Surface lithological constraints (hyperspectral-derived).

By integrating these datasets, we can achieve more robust geological models that are
consistent with both surface observations and subsurface physical properties.

**What are we doing?**

We're performing Bayesian inference to estimate:
- **Dips**: Orientation angles of geological layers.
- **Densities**: Rock densities for different units.

The joint approach leverages the strengths of each data type: EnMap provides high-resolution
surface mapping, while gravity data constrains the 3D geometry of bodies at depth.

**The Joint Bayesian Framework**

In joint inversion, we seek the posterior distribution of parameters θ given multiple
datasets y₁ (gravity) and y₂ (EnMap):

.. math::

    p(θ|y_1, y_2) \propto p(y_1|θ) p(y_2|θ) p(θ)

Assuming independence between the observation errors of different datasets, the joint
likelihood is the product of individual likelihoods.
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
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.joint_probabilistic_model import joint_set_priors, generate_joint_likelihood
from mineye.GeoModel.model_one.visualization import (
    gempy_viz, 
    plot_many_observed_vs_forward,
    plot_joint_inversion_results
)
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison

# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(1234)

# %%
# Load Data and Model Setup
# -------------------------
#
# First, we need to load the geological model and the observation data.
# We use a base GemPy model and then configure it for the inversion.

# Get paths from configuration
mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
geophysical_dir = paths.get_geophysical_dir()
base_dir = os.path.dirname(mod_or_path) # Assuming data is nearby

# Define extent (matching model 1)
extent = [-707521, -675558, 4526832, 4551949, -500, 505]
refinement = 5

# Create GemPy model
simple_geo_model = gp.create_geomodel(
    project_name='joint_inversion',
    extent=extent,
    refinement=refinement,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={"Tournaisian_Plutonites": ["Tournaisian Plutonites"]}
)

# %%
# Setup Gravity and EnMap Data
# ----------------------------

# 1. Load Gravity Data
gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)

# 2. Load EnMap Data (previously preprocessed)
xyz_enmap = np.load(os.path.join(base_dir, 'central_xyz.npy'))
labels_enmap = np.load(os.path.join(base_dir, 'central_labels.npy'))
labels_enmap[labels_enmap == 2] = 1  # Normalize labels

print(f"Gravity stations: {len(observed_gravity_ugal)}")
print(f"EnMap observations: {len(labels_enmap)}")

# %%
# Configure Model Grid for Joint Observations
# -------------------------------------------
#
# In GemPy, forward models are computed at specific grid locations. For joint inversion,
# we must ensure that the model computes gravity at station locations and lithology 
# at EnMap pixel locations. We combine these into a Custom Grid.

gravity_xyz = np.zeros((len(gravity_data), 3))
gravity_xyz[:, 0] = gravity_data['X'].values
gravity_xyz[:, 1] = gravity_data['Y'].values
gravity_xyz[:, 2] = gravity_data['Z'].values

combined_custom_points = np.vstack([gravity_xyz, xyz_enmap])

# Set the custom grid
gp.set_custom_grid(simple_geo_model.grid, combined_custom_points)
gp.set_active_grid(
    grid=simple_geo_model.grid,
    grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
    reset=True
)

# %%
# Normalization and Priors
# ------------------------

# Normalize gravity baseline
norm_params = normalize(
    baseline_fw_gravity_np=(baseline(simple_geo_model)),
    observed_gravity=observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)

# Define Priors
model_priors = {
    'dips': dist.Normal(
        loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
        scale=torch.tensor(10, dtype=torch.float64),
        validate_args=True
    ).to_event(1),
    'density': dist.Normal(
        loc=(torch.tensor([2.9 - 2.67, 2.3 - 2.67])),
        scale=torch.tensor(0.15),
    ).to_event(1)
}

# %%
# Define Joint Likelihood
# -----------------------
#
# The joint likelihood function handles both data types. It extracts the 
# relevant parts of the forward model solution and compares them to the 
# observations.

likelihood_fn = generate_joint_likelihood(norm_params, n_gravity_points=len(gravity_xyz))

prob_model = gpp.make_gempy_pyro_model(
    priors=model_priors,
    set_interp_input_fn=joint_set_priors,
    likelihood_fn=likelihood_fn,
    obs_name="Joint_Obs"
)

# Prepare observed data tensors
# Note: joint_obs is passed as a list of tensors
gravity_obs_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
enmap_obs_tensor = torch.tensor(labels_enmap, dtype=torch.float64)
joint_obs = [gravity_obs_tensor, enmap_obs_tensor]

# %%
# Prior Predictive Checks
# -----------------------
# Verify that the model and priors can reasonably represent the data.

print("Running joint prior predictive...")
prior_data = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=simple_geo_model,
    y_obs_list=joint_obs,
    n_samples=20,
    plot_trace=True
)

# %%
# MCMC Inference
# --------------
# We use NUTS to sample from the posterior distribution.

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
        num_samples=100,
        warmup_steps=100,
    ),
    plot_trace=True,
    run_posterior_predictive=True
)

data.extend(prior_data)

# %%
# Visualization of Results
# ------------------------
# Finally, we visualize the results of the joint inversion, comparing 
# observed vs. modeled values for both datasets.

plot_joint_inversion_results(
    data=data,
    observed_gravity=observed_gravity_ugal,
    xy_gravity=gravity_xyz[:, :2],
    observed_enmap=labels_enmap,
    geo_model=simple_geo_model,
    n_gravity_points=len(gravity_xyz)
)

# Show cross-sections with uncertainty
gempy_viz(simple_geo_model, data, n_samples=20)

plt.show()
