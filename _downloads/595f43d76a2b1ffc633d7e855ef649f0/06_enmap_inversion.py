r"""
Bayesian EnMap Inversion: Categorical Likelihood and Ordinal Probabilities
==========================================================================

This tutorial demonstrates how to perform Bayesian inversion using surface lithological 
classifications (e.g., from EnMap satellite data) as observations. Unlike gravity 
inversion which uses continuous values, EnMap inversion deals with discrete categorical labels.

**What are we doing?**

We're inferring geological parameters (like layer orientations) by comparing the 
predicted lithology at the surface with observed lithological labels derived from 
remote sensing data.

**Why Bayesian Inference?**

Bayesian inference allows us to:
1.  **Quantify Uncertainty**: Instead of a single "best-fit" model, we get a posterior 
    distribution representing all models consistent with both our prior knowledge and 
    the EnMap observations.
2.  **Incorporate Prior Knowledge**: We can formally include geological information 
    (e.g., typical layer dips) through prior distributions.
3.  **Handle Classification Errors**: The probabilistic approach naturally accounts 
    for potential misclassifications in the remote sensing data.

**The Bayesian Framework**

We seek the posterior distribution of geological parameters :math:`\theta` given 
the observed EnMap labels :math:`y`:

.. math::

    p(\theta|y) = \frac{p(y|\theta) p(\theta)}{p(y)}

Where:
- :math:`p(\theta)` is the **prior**, encoding our initial knowledge.
- :math:`p(y|\theta)` is the **likelihood**, describing the probability of the 
  labels given a specific geological configuration.
- :math:`p(y)` is the **evidence** (marginal likelihood).

**The Forward Model**

The relationship between parameters and observations is given by:

.. math::

    y = f(\theta) + \epsilon

Where :math:`f(\theta)` is the GemPy geological model that predicts rock types 
at the surface, and :math:`\epsilon` represents the observation noise or 
classification uncertainty.

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

# %%
# Import Libraries
# ----------------
#
# We load the required scientific packages, GemPy ecosystem, and Pyro for 
# probabilistic programming.

import os
import numpy as np
import matplotlib.pyplot as plt
import torch
import pyro
import pyro.distributions as dist
import arviz as az
import inspect
from pathlib import Path

# GemPy ecosystem
import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv

# Configuration and helpers
from mineye.config import paths
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations
from mineye.GeoModel.model_one.visualization import (
    gempy_viz,
    plot_probability_heatmap,
    compute_probability_density_fields,
    plot_many_observed_vs_forward
)
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import enmap_likelihood_fn

# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(seed)

# %%
# Step 1: Model Configuration
# ---------------------------
#
# We define the spatial extent of our study area and the resolution of the 
# geological model.
#
# **Model Extent**:
# - X: -707521 to -675558
# - Y: 4526832 to 4551949
# - Z: -500 to 505 (meters)

# Define Model Extent (matching the gravity example area)
extent = [-707521, -675558, 4526832, 4551949, -500, 505]

# %%
# Step 2: Load EnMap Observations
# -------------------------------
#
# **What is EnMap data?**
#
# EnMap (Environmental Mapping and Analysis Program) provides hyperspectral 
# satellite imagery. In this example, we use a lithological classification 
# derived from such data.
#
# **Units**: Discrete categorical labels (0, 1, 2...).

# Define paths to find the central_xyz files in the same folder as this script
current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
xyz_path = os.path.join(current_dir, 'central_xyz.npy')
labels_path = os.path.join(current_dir, 'central_labels.npy')

# Load EnMap data
xyz_enmap = np.load(xyz_path)
labels_enmap = np.load(labels_path)
labels_enmap[labels_enmap == 2] = 1  # Normalize labels to binary for this example

print(f"\nEnMap observations:")
print(f"  Number of measurements: {len(xyz_enmap)}")
print(f"  Classes present: {np.unique(labels_enmap)}")
print(f"  Class 0 count: {np.sum(labels_enmap == 0)}")
print(f"  Class 1 count: {np.sum(labels_enmap == 1)}")

# %%
# Step 3: Initial Geological Model Setup
# --------------------------------------
#
# We create the base GemPy model using geological orientations and points.

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

# Compute initial model
gp.compute_model(simple_geo_model)

# Visualize initial model
gpv.plot_2d(simple_geo_model, direction='z', show_data=True)
plt.title("Initial Geological Model (Surface)")
plt.show()

# %%
# Step 4: Configure Custom Grid for Observations
# ----------------------------------------------
#
# We set a custom grid at the exact XYZ locations where we have EnMap labels. 
# This ensures GemPy predicts rock types exactly where we have observations.

gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
gp.set_active_grid(
    grid=simple_geo_model.grid,
    grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
    reset=True
)

# %%
# Step 5: Baseline Forward Model
# ------------------------------
#
# We compute the model at the custom grid locations to establish a baseline 
# comparison with the observations.

# Compute model on custom grid
sol = gp.compute_model(simple_geo_model)

# The results at custom grid are stored in solutions
output_center = sol.octrees_output[0].last_output_center
scalar_field_at_custom_grid = output_center.exported_fields.scalar_field[output_center.grid.custom_grid_slice]

print(f"\nBaseline Forward Model:")
print(f"  Scalar field range: {scalar_field_at_custom_grid.min():.2f} to {scalar_field_at_custom_grid.max():.2f}")

# %%
# Step 6: Define Prior Distributions
# -----------------------------------
#
# **Prior Knowledge in Bayesian Inference**
#
# Priors encode our geological knowledge *before* seeing the data. Here we 
# specify our uncertainty about the dip of the geological layers.
#
# **Prior on Dips**
#
# We use a Normal distribution:
#
# .. math::
#
#     \text{dips} \sim \mathcal{N}(10^\circ, 10^\circ)
#
# **Interpretation**:
# - **Mean (μ)**: 10 degrees. Our best guess for the regional dip.
# - **Std (σ)**: 10 degrees. Represents significant uncertainty in our knowledge.

model_priors = {
        'dips': dist.Normal(
            loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
            scale=torch.tensor(10, dtype=torch.float64)
        ).to_event(1)
}

def set_priors_enmap(samples, geo_model):
    return _modify_orientations(samples=samples, geo_model=geo_model, key="dips")

# %%
# Step 7: Inspecting the Likelihood Function
# ------------------------------------------
#
# **Understanding the Likelihood Function**
#
# The likelihood function is crucial for connecting model predictions to 
# observations. For EnMap data, we use an ordinal classification likelihood.
# Let's inspect what it does internally:

# Display the function signature and docstring
help(enmap_likelihood_fn)

# %%
# **Full Source Code**
#
# For complete transparency, here's the full implementation of the 
# EnMap likelihood function:

print("Likelihood function source code:")
print("=" * 60)
print(inspect.getsource(enmap_likelihood_fn))

# %%
# Step 8: Probabilistic Model Creation
# ------------------------------------
#
# We combine the priors, the geological model update function, and the 
# likelihood function into a single Pyro model.
#
# **Internal Workflow**:
# 1. Sample :math:`\theta` (dips) from priors.
# 2. Update GemPy model with sampled :math:`\theta`.
# 3. Compute forward model (predict scalar fields at surface).
# 4. Apply `enmap_likelihood_fn` to compute :math:`p(data|\theta)`.

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
# Step 9: Prior Predictive Check
# ------------------------------
#
# Before running expensive inference, we check if our priors allow for the 
# observed rock types. If the prior predictions don't even cover the 
# observations, the model might be poorly specified.

print("Running prior predictive...")
prior_data = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=simple_geo_model,
    y_obs_list=None,
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
# Step 10: MCMC Inference using NUTS
# ----------------------------------
#
# We use the No-U-Turn Sampler (NUTS) to sample from the posterior distribution. 
# NUTS is an efficient Hamiltonian Monte Carlo algorithm that avoids the need 
# for manual tuning of step sizes.

RUN_SIMULATION = False
if RUN_SIMULATION:
    print("Running NUTS inference...")
    data = gpp.run_nuts_inference(
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

    data.extend(prior_data)
    data.to_netcdf("arviz_data_enmap.nc")
else:
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
    print(f"Loaded inference results from {data_path}")

# %%
# Step 11: Posterior Analysis
# ---------------------------
#
# We visualize how the data has informed our knowledge of the layer dips and 
# assess the reduction in uncertainty.

# Posterior Density Plot
az.plot_density(
    data=[data, data.prior],
    var_names=["dips"],
    data_labels=["Posterior", "Prior"],
    colors=[default_red, default_blue],
    shade=0.2
)
plt.title("Update of Dip Distributions")
plt.show()

# Uncertainty Reduction Calculation
prior_std = data.prior["dips"].values.std()
post_std = data.posterior["dips"].values.std()
reduction = (1 - post_std / prior_std) * 100
print(f"\nUncertainty Reduction in Dips: {reduction:.1f}%")

# Visualize the resulting geological model uncertainty
print("\nVisualizing geological uncertainty...")
gempy_viz(simple_geo_model, data)

# %%
# Step 12: Summary and Conclusions
# --------------------------------
#
# **Workflow Recap**
#
# 1. **Prior Definition**: Normal distribution on layer dips.
# 2. **Likelihood Construction**: Ordinal classification using sigmoids.
# 3. **Model Integration**: Pyro/GemPy bridge.
# 4. **Inference**: NUTS sampling.
# 5. **Validation**: Posterior density and geological visualization.
#
# **Key Insights**:
# - Categorical data can be used for inversion by mapping continuous fields to 
#   probabilities via sigmoid functions.
# - Differentiability is preserved, allowing for efficient gradient-based sampling.
# - Remote sensing observations can significantly reduce uncertainty in subsurface 
#   geometrical parameters.
# Final workflow visualization (ASCII art)
# 
# ┌─────────────────────────────────────────────────────────────────┐
# │ 1. PRIOR DEFINITION                                             │
# │    dips ~ Normal(10, 10)                                        │
# └─────────────────────────────────────────────────────────────────┘
#                                ↓
# ┌─────────────────────────────────────────────────────────────────┐
# │ 2. PROBABILISTIC MODEL CONSTRUCTION                             │
# │    Bridge GemPy and Pyro with enmap_likelihood_fn               │
# └─────────────────────────────────────────────────────────────────┘
#                                ↓
# ┌─────────────────────────────────────────────────────────────────┐
# │ 3. MCMC INFERENCE (NUTS)                                        │
# │    Sample posterior p(dips | EnMap Labels)                      │
# └─────────────────────────────────────────────────────────────────┘
#                                ↓
# ┌─────────────────────────────────────────────────────────────────┐
# │ 4. POSTERIOR ANALYSIS                                           │
# │    Compare prior vs posterior distributions                     │
# └─────────────────────────────────────────────────────────────────┘
# sphinx_gallery_thumbnail_number = 2