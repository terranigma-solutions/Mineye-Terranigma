"""
Joint Bayesian Inversion: Gravity and EnMap
===========================================

This tutorial demonstrates a joint Bayesian inversion combining two disparate data types:
1. **Gravity Data**: Volumetric measurements sensitive to subsurface density contrasts.
2. **EnMap Satellite Data**: Surface lithology classifications sensitive to the intersection 
   of geological units with topography.

Joint inversion is a powerful technique to reduce non-uniqueness in geophysical models. By 
combining datasets that have different sensitivities, we can better constrain the 
subsurface architecture.

**Key Concepts in this Tutorial:**

1. **Multi-Grid Configuration**: How to handle different observation grids (Centered grid for 
   gravity, Custom grid for EnMap) within a single GemPy model.
2. **Joint Likelihood Formulation**: Defining a Pyro model with multiple likelihood functions.
3. **Likelihood Balance Check**: A critical diagnostic step to ensure one dataset doesn't 
   statistically "drown out" the other due to differences in data volume or noise assumptions.
4. **Joint Parameter Inference**: Simultaneously estimating layer dips and rock densities.

**The Likelihood Balance Problem**

In Bayesian joint inversion, we seek the posterior:

.. math::

    p(θ|y_{grav}, y_{enmap}) ∝ p(y_{grav}|θ) p(y_{enmap}|θ) p(θ)

In log-space:

.. math::

    \log p(θ|y_{total}) = \log p(y_{grav}|θ) + \log p(y_{enmap}|θ) + \log p(θ)

If the magnitude of :math:`\log p(y_{grav}|θ)` is significantly larger than 
:math:`\log p(y_{enmap}|θ)` (e.g., -10,000 vs -100), the optimizer/sampler will primarily 
focus on satisfying the gravity data, effectively ignoring the EnMap data. This "imbalance" 
can happen because:
- One dataset has many more observation points.
- The noise standard deviation (:math:`σ`) for one dataset is assumed to be too small.
- The measurement units or normalization scales are inconsistent.

We'll demonstrate how to check and address this balance before running full inference.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import pyro
import pyro.distributions as dist
import arviz as az
from pathlib import Path
import inspect

import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv
from gempy_probability.core.samplers_data import NUTSConfig

# Import Mineye-specific helpers
from mineye.config import paths
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel, read_gravity
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.joint_probabilistic_model import joint_set_priors
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import (
    generate_multigravity_likelihood_hierarchical_per_station,
    enmap_likelihood_fn
)
from mineye.GeoModel.model_one.inference_diagnostics import check_likelihood_balance
from mineye.GeoModel.model_one.visualization import (
    generate_gravity_uncertainty_plots,
    gempy_viz,
    plot_probability_heatmap
)

# Set random seeds for reproducibility
seed = 42
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# %%
# 1. Setup Model Extent and Load Data
# ----------------------------------

extent = [-707521, -675558, 4526832, 4551949, -500, 505]
refinement = 5

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
geophysical_dir = paths.get_geophysical_dir()

# Create Initial GemPy Model
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
# 2. Setup Gravity Observations
# ----------------------------
gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100

# %%
# 3. Setup EnMap Observations (Custom Grid)
# ----------------------------------------
# EnMap data is provided as XYZ points and categorical rock-type labels.
# We add these to the GemPy model as a 'Custom Grid'.

# For this example, we'll assume the files exist in the same folder as this script
current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
xyz_path = os.path.join(current_dir, 'central_xyz.npy')
labels_path = os.path.join(current_dir, 'central_labels.npy')

if not os.path.exists(xyz_path):
    # Fallback to dummy data if real data isn't found for the sphinx build
    print("Warning: Real EnMap data not found. Using dummy data for demonstration.")
    xyz_enmap = np.array([[-690000, 4540000, 500]])
    labels_enmap = np.array([1])
else:
    xyz_enmap = np.load(xyz_path)
    labels_enmap = np.load(labels_path)
    labels_enmap[labels_enmap == 2] = 1  # Normalize labels to match model units

# Configure Custom Grid for EnMap
gp.set_custom_grid(geo_model.grid, xyz_enmap)

# Activate both Centered (Gravity) and Custom (EnMap) grids
gp.set_active_grid(
    grid=geo_model.grid,
    grid_type=[geo_model.grid.GridTypes.CENTERED, geo_model.grid.GridTypes.CUSTOM],
    reset=True
)

# %%
# 4. Normalization and Likelihood Setup
# ------------------------------------

# Normalize gravity based on a baseline model run
norm_params = normalize(
    baseline_fw_gravity_np=baseline(geo_model),
    observed_gravity=observed_gravity_ugal,
    method="align_to_reference"
)

# Define Priors
model_priors = {
        'dips'   : dist.Normal(
            loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),
            scale=torch.tensor(10, dtype=torch.float64)
        ),
        'density': dist.Normal(
            loc=(torch.tensor([2.9 - 2.67, 2.3 - 2.67])),
            scale=torch.tensor(0.15),
        ).to_event(1)
}

# Define Joint Likelihood
gravity_likelihood = generate_multigravity_likelihood_hierarchical_per_station(norm_params)
enmap_likelihood = enmap_likelihood_fn

# Create the Pyro model
prob_model = gpp.make_gempy_pyro_model(
    priors=model_priors,
    set_interp_input_fn=joint_set_priors,
    likelihood_fn=[gravity_likelihood, enmap_likelihood],
    obs_name="Joint_Obs"
)

# Prepare observed data tensors
gravity_obs_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
enmap_obs_tensor = torch.tensor(labels_enmap, dtype=torch.float64)
joint_obs = [gravity_obs_tensor, enmap_obs_tensor]

# %%
# 5. Check Likelihood Balance
# --------------------------
# This is a crucial step! We calculate the log-probabilities of each dataset 
# at the starting parameters to see if one dominates.

print("\n--- Performing Likelihood Balance Check ---")
check_likelihood_balance(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=joint_obs
)

# %%
# 6. Run Joint Inference (NUTS)
# ----------------------------


RUN_SIMULATION = False

if RUN_SIMULATION:

    print("Running joint prior predictive...")
    prior_data = gpp.run_predictive(
        prob_model=prob_model,
        geo_model=simple_geo_model,
        y_obs_list=joint_obs,
        n_samples=200,
        plot_trace=True
    )
    
    print("\nRunning Joint NUTS Inference...")
    # Reduced samples for tutorial speed; increase for production
    data = gpp.run_nuts_inference(
        prob_model=prob_model,
        geo_model=geo_model,
        y_obs_list=joint_obs,
        config=NUTSConfig(
            num_samples=200,
            warmup_steps=200,
            target_accept_prob=0.8
        ),
        plot_trace=True,
        run_posterior_predictive=True
    )

else:

    # Get the directory of the current file
    current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
    data_path = current_dir / "arviz_data_joint_feb2026.nc"

    if not data_path.exists():
        raise FileNotFoundError(
            f"Data file not found at {data_path}. "
            f"Please run the simulation first with RUN_SIMULATION=True"
        )

    # Read the data file
    data = az.from_netcdf(str(data_path))
    print(f" Loaded inference results from {data_path}")

# %%
# 7. Analyze Results
# -----------------

# Trace Plot
az.plot_trace(data, var_names=["dips", "density"])
plt.show()

# Posterior Distributions
az.plot_posterior(data, var_names=["density"])
plt.title("Posterior Densities (Joint Inversion)")
plt.show()

# Gravity Uncertainty Plot
# This shows how well the posterior models the observed gravity data
# response = r'$\mu_{gravity}$'
# generate_gravity_uncertainty_plots(
#     gravity_samples_norm=data.posterior_predictive[response].values[0, :],
#     observed_gravity_ugal=observed_gravity_ugal,
#     xy_ravel=xy_ravel
# )

# Probability Heatmap (EnMap predictive check)
plot_probability_heatmap(data, 'posterior_predictive')
plt.show()

# 3D Visualization of the inferred model
# gempy_viz(geo_model, data)

# %%
# Summary of the Joint Likelihood Function
# ----------------------------------------
# Internally, when we pass a list of likelihood functions to `make_gempy_pyro_model`,
# it creates multiple observation sites. 
#
# If `obs_name` is "Joint_Obs", the sites will be named:
# 1. `Joint_Obs_0`: Receives the output of the first likelihood (Gravity).
# 2. `Joint_Obs_1`: Receives the output of the second likelihood (EnMap).
#
# The Likelihood Balance Check looks for these specific nodes in the Pyro trace
# and compares their total `log_prob`. If the ratio is outside [0.01, 100], 
# it's a sign that you might need to:
# - Increase the noise (:math:`σ`) of the dominant dataset.
# - Decrease the noise of the ignored dataset.
# - Re-examine your normalization strategy.
