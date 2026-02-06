"""
Bayesian Joint Inversion: Gravity and EnMap
===========================================

This tutorial demonstrates a joint Bayesian inversion combining two disparate data types:
1. **Gravity Data**: Volumetric measurements sensitive to subsurface density contrasts.
2. **EnMap Satellite Data**: Surface lithology classifications sensitive to the intersection 
   of geological units with topography.

Joint inversion is a powerful technique to reduce non-uniqueness in geophysical models. By 
combining datasets that have different sensitivities, we can better constrain the 
subsurface architecture.

**What are we doing?**

We're performing Bayesian inference on a geological model using both gravity measurements 
and EnMap surface classifications. The goal is to estimate geological parameters such as:

- **Dips**: Orientation angles of geological layers (degrees).
- **Densities**: Rock densities for different lithological units (g/cm³).

**Why Bayesian Inference?**

Bayesian inference provides several advantages over traditional deterministic approaches:

1. **Posterior probability distributions**: Not just point estimates, but full probability 
   distributions that quantify uncertainty in our parameter estimates.
2. **Uncertainty quantification**: Explicit characterization of what we know and don't know 
   about the subsurface.
3. **Prior knowledge incorporation**: Ability to incorporate geological expertise and 
   physical constraints into the inversion.
4. **Data Integration**: A natural framework for combining multiple datasets with different 
   noise characteristics and measurement scales.

**The Bayesian Framework**

In Bayesian joint inversion, we seek the posterior distribution of parameters θ given 
multiple datasets $y_{grav}$ and $y_{enmap}$ using Bayes' theorem:

# .. math::
#
#     p(\theta | y_{grav}, y_{enmap}) \propto p(y_{grav} | \theta) p(y_{enmap} | \theta) p(\theta)
#
# In log-space, this becomes an additive process:
#
# .. math::
#
#     \log p(\theta | y_{total}) = \log p(y_{grav} | \theta) + \log p(y_{enmap} | \theta) + \log p(\theta)

**The Forward Model**

We have two distinct forward models mapping the same parameters θ to different observation spaces:

1. **Gravity Forward Model**: $y_{grav} = f_{grav}(\theta) + \epsilon_{grav}$
   Maps density distributions to gravitational acceleration.
2. **EnMap Forward Model**: $y_{enmap} = f_{enmap}(\theta) + \epsilon_{enmap}$
   Maps unit boundaries to surface lithology classifications.

**Key Concepts in this Tutorial:**

1. **Multi-Grid Configuration**: How to handle different observation grids (Centered grid for 
   gravity, Custom grid for EnMap) within a single GemPy model.
2. **Joint Likelihood Formulation**: Defining a Pyro model with multiple likelihood functions.
3. **Likelihood Balance Check**: A critical diagnostic step to ensure one dataset doesn't 
   statistically "drown out" the other.
4. **Joint Parameter Inference**: Simultaneously estimating layer dips and rock densities.
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
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality, check_likelihood_balance
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
#
# **Model Configuration**
#
# We define the spatial extent of our geological model and the resolution (refinement).
# Higher refinement leads to more accurate results but increases computation time.

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
#
# **What is gravity data?**
#
# Gravity data measures small variations in the Earth's gravitational field caused by 
# density differences in the subsurface.
#
# **Units**: μGal (microgal)

gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
print(f"\nGravity observations:")
print(f"  Number of measurements: {len(observed_gravity_ugal)}")
print(f"  Range: {observed_gravity_ugal.min():.1f} to {observed_gravity_ugal.max():.1f} μGal")
print(f"  Mean: {observed_gravity_ugal.mean():.1f} μGal")

geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100

# %%
# 3. Setup EnMap Observations (Custom Grid)
# ----------------------------------------
#
# **What is EnMap data?**
#
# EnMap (Environmental Mapping and Analysis Program) is a hyperspectral satellite mission. 
# In this context, we use it as a categorical surface classification of rock types.
#
# **Data Representation**: Categorical labels mapping to geological units.

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

print(f"\nEnMap observations:")
print(f"  Number of pixels: {len(labels_enmap)}")
print(f"  Categories: {np.unique(labels_enmap)}")

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
#
# **Gravity Normalization**
#
# Since gravity data often has a regional component, we align our forward model to 
# the observations using a baseline model run.

norm_params = normalize(
    baseline_fw_gravity_np=baseline(geo_model),
    observed_gravity=observed_gravity_ugal,
    method="align_to_reference"
)

# %%
# **Step 5: Define Prior Distributions**
# -----------------------------------
#
# **Prior Knowledge in Bayesian Inference**
#
# Priors encode our geological knowledge *before* seeing the data. We define priors 
# for layer dips and rock densities.
#
# **Prior on Dips**
#
# We use a Normal distribution for orientations:
#
# .. math::
#
#     \text{dips} \sim \mathcal{N}(\mu=10^\circ, \sigma=10^\circ)
#
# **Prior on Density**
#
# We define density contrasts relative to a background density (2.67 g/cm³):
#
# .. math::
#
#     \text{density} \sim \mathcal{N}(\mu, \sigma=0.15)

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

# %%
# **Step 6: Define Joint Likelihood**
# ----------------------------------
#
# The joint likelihood is the product of individual likelihoods:
#
# .. math::
#
#     p(y_{grav}, y_{enmap}|\theta) = p(y_{grav}|\theta) \cdot p(y_{enmap}|\theta)
#
# **Gravity Likelihood**: Hierarchical per-station noise model.
# **EnMap Likelihood**: Categorical distribution for surface units.

gravity_likelihood = generate_multigravity_likelihood_hierarchical_per_station(norm_params)
enmap_likelihood = enmap_likelihood_fn

# %%
# **Inspecting the Likelihood Functions**
#
# To understand exactly what happens inside, we can inspect the source code:

print("Gravity Likelihood Function:")
print("=" * 30)
print(inspect.getsource(generate_multigravity_likelihood_hierarchical_per_station))

print("\nEnMap Likelihood Function:")
print("=" * 30)
print(inspect.getsource(enmap_likelihood_fn))

# %%
# **Probabilistic Model Creation**
# -------------------------------
#
# We combine everything using `make_gempy_pyro_model`. 
#
# **Internal Workflow Diagram**:
#
# .. code-block:: text
#
#     Sample θ (Dips, Densities)
#            ↓
#     joint_set_priors(θ) → Update GeoModel
#            ↓
#     Forward Model Computation (Gravity + Custom Grid)
#            ↓
#     [Likelihood 1 (Gravity), Likelihood 2 (EnMap)]
#            ↓
#     p(data|\theta) → MCMC Proposes Next \theta

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
#
# We use the No-U-Turn Sampler (NUTS), a state-of-the-art MCMC algorithm.
#
# **Configuration Parameters**:
# - `num_samples`: Number of posterior samples to collect.
# - `warmup_steps`: Iterations for the sampler to adapt its parameters.
# - `target_accept_prob`: Target acceptance rate (0.8 is usually a good balance).
# - `step_size`: Initial integration step size.

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
    print(f" Loaded inference results from {data_path}")

# %%
# 7. Analyze Results
# -----------------
#
# **Posterior Analysis**
#
# We evaluate the quality of our inference by looking at the trace plots and 
# posterior distributions.
#
# **Trace Plot**: Should show "fuzzy caterpillar" chains indicating good mixing.
# **Posterior Distributions**: Should be narrower than the priors, indicating 
# that the data has informed our knowledge.

# Trace Plot
az.plot_trace(data, var_names=["dips", "density"])
plt.show()

# Posterior Distributions
az.plot_posterior(data, var_names=["density"])
plt.title("Posterior Densities (Joint Inversion)")
plt.show()

# %%
# **Predictive Checks**
#
# We check how well our posterior models can replicate the observed data.

# Probability Heatmap (EnMap predictive check)
# This shows the probability of each unit at the surface.
plot_probability_heatmap(data, 'posterior_predictive')
plt.show()

# %%
# **Summary and Conclusions**
# --------------------------
#
# In this tutorial, we've successfully:
# 1. Configured a GemPy model with multiple grids (Centered for Gravity, Custom for EnMap).
# 2. Defined a joint Bayesian model with multiple likelihood functions.
# 3. Performed a Likelihood Balance Check to ensure fair data integration.
# 4. Inferred subsurface dips and densities from joint observations.

# **Joint Inversion Workflow Recap:**
#
#
#     1. Data Preparation (Gravity + EnMap)
#     2. Grid Setup (Multi-grid)
#     3. Prior Definition (Geological constraints)
#     4. Likelihood Formulation (Noise models)
#     5. Balance Diagnostic (Crucial for joint!)
#     6. MCMC Inference (NUTS)
#     7. Posterior Analysis & Visualization
# sphinx_gallery_thumbnail_number = 3
