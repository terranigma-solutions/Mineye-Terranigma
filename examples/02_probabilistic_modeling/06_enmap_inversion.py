"""
Bayesian EnMap Inversion: Categorical Likelihood and Ordinal Probabilities
==========================================================================

This tutorial demonstrates how to perform Bayesian inversion using surface lithological 
classifications (e.g., from EnMap satellite data) as observations. Unlike gravity 
inversion which uses continuous values, EnMap inversion deals with discrete categorical labels.

**What are we doing?**

We're inferring geological parameters (like layer orientations) by comparing the 
predicted lithology at the surface with observed lithological labels derived from 
remote sensing data.

**The Likelihood: Categorical with Ordinal Probabilities**

Since the observed data consists of discrete rock type labels (e.g., 0 for Unit A, 1 for Unit B), 
 we use a **Categorical likelihood**. However, GemPy's forward model is differentiable 
and produces a continuous scalar field. To bridge this gap, we use an **Ordinal Probabilities** 
approach.

1.  **Scalar Field to Probabilities**: We take the continuous scalar field value at the surface.
2.  **Sigmoid Boundaries**: We define boundaries (thresholds) between units. We use 
    sigmoid functions centered at these boundaries to compute the probability of a 
    point belonging to each unit.
    
    .. math::
        P(\text{below boundary}) = \sigma\left(\frac{\text{boundary} - \text{scalar}}{\text{temperature}}\right)

3.  **Softmax-like behavior**: The difference between cumulative probabilities at 
    successive boundaries gives the probability mass for each discrete class.
4.  **Temperature**: A 'temperature' parameter controls the sharpness of these 
    probabilistic transitions. Lower temperature means sharper, more deterministic boundaries.

**Why this approach?**

- **Differentiability**: The sigmoid-based mapping is fully differentiable, which 
  is crucial for Gradient-based MCMC methods like **NUTS** (No-U-Turn Sampler).
- **Uncertainty**: It naturally accounts for classification uncertainty near 
  geological interfaces.

"""

import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import pyro
import pyro.distributions as dist
import arviz as az
import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue

# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# %%
# Import helper functions and data paths
# --------------------------------------

from mineye.config import paths
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations
from mineye.GeoModel.model_one.visualization import (
    gempy_viz,
    plot_probability_heatmap,
    compute_probability_density_fields
)
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import enmap_likelihood_fn

# %%
# Step 1: Load Data and Setup Model
# ---------------------------------
# We load the orientations and points for our base geological model, 
# and the EnMap extracted surface points and labels.

# Define paths to find the central_xyz files in the same folder as this script
base_dir = os.path.dirname(os.path.abspath(__file__))
xyz_path = os.path.join(base_dir, 'central_xyz.npy')
labels_path = os.path.join(base_dir, 'central_labels.npy')

# Load EnMap data
xyz_enmap = np.load(xyz_path)
labels_enmap = np.load(labels_path)
labels_enmap[labels_enmap == 2] = 1  # Normalize labels to binary for this example

print(f"Loaded {len(xyz_enmap)} EnMap observation points.")

# Define Model Extent (matching the gravity example area)
extent = [-707521, -675558, 4526832, 4551949, -500, 505]

# Create Initial GemPy Model
simple_geo_model = gp.create_geomodel(
    project_name='enmap_inversion',
    extent=extent,
    refinement=5,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=paths.get_orientations_path(),
        path_to_surface_points=paths.get_points_path(),
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={"Tournaisian_Plutonites": ["Tournaisian Plutonites"]}
)

# %%
# Step 2: Configure Custom Grid for Observations
# ----------------------------------------------
# We set a custom grid at the exact XYZ locations where we have EnMap labels.

gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
gp.set_active_grid(
    grid=simple_geo_model.grid,
    grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
    reset=True
)

# %%
# Step 3: Define Priors and Probabilistic Model
# ---------------------------------------------
# We set a prior distribution on the dip of the geological layers.

model_priors = {
        'dips': dist.Normal(
            loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
            scale=torch.tensor(10, dtype=torch.float64)
        )
}


def set_priors_enmap(samples, geo_model):
    return _modify_orientations(samples=samples, geo_model=geo_model, key="dips")


# Create the Pyro model
prob_model = gpp.make_gempy_pyro_model(
    priors=model_priors,
    set_interp_input_fn=set_priors_enmap,
    likelihood_fn=enmap_likelihood_fn,
    obs_name="EnMap Labels"
)

# Prepare observed data as a torch tensor
labels_enmap_tensor = torch.tensor(labels_enmap, dtype=torch.float64)

# %%
# Step 4: Prior Predictive Check
# ------------------------------
# Before running inference, we check if our priors allow for the observed rock types.

print("Running prior predictive...")
prior_data = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=simple_geo_model,
    y_obs_list=[labels_enmap_tensor],
    n_samples=100
)

# Plotting the prior distribution of labels vs observed
pred_labels = prior_data.prior['EnMap Labels'][0]
plt.figure(figsize=(8, 4))
plt.hist([pred_labels.values.flatten(), labels_enmap],
         label=['Prior Predicted', 'Observed Ground Truth'],
         bins=np.arange(3) - 0.5, density=True, rwidth=0.8)
plt.xticks([0, 1])
plt.xlabel("Rock Type Label")
plt.ylabel("Density")
plt.title("Prior Predictive Check")
plt.legend()
plt.show()

# %%
# Step 5: Bayesian Inference using NUTS
# -------------------------------------
# We use the No-U-Turn Sampler (NUTS) to sample from the posterior distribution.

RUN_SIMULATION = False
if RUN_SIMULATION:
    print("Running NUTS inference...")
    inference_data = gpp.run_nuts_inference(
        prob_model=prob_model,
        geo_model=simple_geo_model,
        y_obs_list=labels_enmap_tensor,
        config=NUTSConfig(
            num_samples=200,
            warmup_steps=200,
            target_accept_prob=0.65,
            max_tree_depth=5
        ),
        plot_trace=True,
        run_posterior_predictive=True
    )
else:
    from pathlib import Path
    import inspect

    # Get the directory of the current file
    current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
    data_path = current_dir / "arviz_data_enmap_feb2026.nc"

    if not data_path.exists():
        raise FileNotFoundError(
            f"Data file not found at {data_path}. "
            f"Please run the simulation first with RUN_SIMULATION=True"
        )

    # Read the data file
    data = az.from_netcdf(str(data_path))
    print(f" Loaded inference results from {data_path}")

# %%
# Step 6: Posterior Analysis
# --------------------------
# We visualize how the data has informed our knowledge of the layer dips.

az.plot_density(
    data=[inference_data, prior_data],
    var_names=["dips"],
    data_labels=["Posterior", "Prior"],
    colors=[default_red, default_blue],
    shade=0.2
)
plt.title("Update of Dip Distributions")
plt.show()

# Visualize the resulting geological model uncertainty
gempy_viz(simple_geo_model, inference_data)

# %%
# Explanation of the EnMap Likelihood Function
# -------------------------------------------
#
# The `enmap_likelihood_fn` used above performs the following:
#
# 1.  **Extraction**: It gets the scalar field values from GemPy at the custom grid locations.
# 2.  **Boundaries**: It identifies the scalar values corresponding to the interfaces 
#     (isovalues).
# 3.  **Ordinal Proportions**: It calls `_get_ordinal_probs` which:
#     - Uses a sigmoid function to create 'soft' transitions between units.
#     - Computes $P(Unit_i) = P(\text{below boundary}_{i+1}) - P(\text{below boundary}_i)$.
#     - Ensures probabilities sum to 1.
# 4.  **Categorical Distribution**: Returns a `pyro.distributions.Categorical` 
#     distribution parameterized by these probabilities.
#
# This allows us to calculate the likelihood of the discrete EnMap labels given 
# the continuous geological model.
