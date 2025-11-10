"""
Bayesian Gravity Inversion: Complete Workflow
==============================================

This tutorial provides a comprehensive walkthrough of Bayesian gravity inversion using
the GemPy-Probability framework. We'll explore both the theoretical foundations and
practical implementation of inferring geological parameters from gravity observations.

**What are we doing?**

We're performing Bayesian inference on a geological model using gravity measurements.
The goal is to estimate geological parameters such as:

- **Dips**: Orientation angles of geological layers (degrees)
- **Densities**: Rock densities for different lithological units (g/cm³)

**Why Bayesian Inference?**

Bayesian inference provides several advantages over traditional deterministic approaches:

1. **Posterior probability distributions**: Not just point estimates, but full probability
   distributions that quantify uncertainty in our parameter estimates.

2. **Uncertainty quantification**: Explicit characterization of what we know and don't know
   about the subsurface.

3. **Prior knowledge incorporation**: Ability to incorporate geological expertise, physical
   constraints, and previous studies into the inversion.

4. **Outlier detection**: Identification of problematic measurements or model misspecification
   through hierarchical modeling.

**The Bayesian Framework**

In Bayesian inference, we update our beliefs about parameters θ (e.g., dips, densities)
given observed data y (gravity measurements) using Bayes' theorem:

.. math::

    p(θ|y) = \\frac{p(y|θ) p(θ)}{p(y)}

Where:

- **p(θ|y)**: Posterior distribution - our updated beliefs after seeing the data
- **p(y|θ)**: Likelihood - probability of observing the data given parameters
- **p(θ)**: Prior distribution - our initial beliefs before seeing data
- **p(y)**: Evidence - normalizing constant (marginal likelihood)

**The Forward Model**

Unlike linear regression where the relationship between parameters and observations is
simple (e.g., y = θ₁x + θ₂), here we have a complex forward model:

.. math::

    y = f(θ) + ε

Where:

- **y**: Observed gravity at measurement stations
- **f(θ)**: Complex forward model combining geological interpolation and gravity computation
- **θ**: Parameters (dips, densities, orientations)
- **ε**: Measurement noise (we'll model this hierarchically)

The forward model f involves:

1. Geological interpolation (potential field method)
2. Computing density distributions
3. Forward gravity calculation at observation points
4. Normalization and alignment to observations

sphinx_gallery_thumbnail_num = 5

"""
import sys

# %%
# Import Libraries
# ----------------

import numpy as np
from whoami import whoami
import multiprocessing as mp

import gempy as gp
import gempy_probability as gpp

import torch
import pyro
import pyro.distributions as dist

import arviz as az
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.model_definition.prob_model_factory_v2 import GemPyPyroConfig, GemPyPyroModelExtended

import dotenv

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
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multigravity_likelihood_diagonal, MultiGravityDiagonalLikelihood
from mineye.GeoModel.model_one.probabilistic_model_diagnostics import trace_pyro_model

dotenv.load_dotenv()

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
# Step 1: Load Gravity Observations
# ----------------------------------
#
# **What is gravity data?**
#
# Gravity measurements record variations in Earth's gravitational field at the surface.
# These variations arise from lateral density contrasts in the subsurface. Dense rocks
# (e.g., mafic intrusions) create positive anomalies, while less dense rocks (e.g.,
# sediments) create negative anomalies.
#
# **Units**: Gravity is measured in microGals (µGal), where 1 µGal = 10⁻⁸ m/s².
# Typical gravity anomalies from geological structures range from tens to thousands of µGal.
#
# The data includes:
#
# - **gravity_data**: Station locations (x, y, z coordinates)
# - **observed_gravity_ugal**: Measured gravity values in microGals

gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
print(f"\nGravity observations:")
print(f"  Number of measurements: {len(observed_gravity_ugal)}")
print(f"  Range: {observed_gravity_ugal.min():.1f} to {observed_gravity_ugal.max():.1f} µGal")
print(f"  Mean: {observed_gravity_ugal.mean():.1f} µGal")

# %%
# Step 2: Setup Geomodel with Gravity Configuration
# --------------------------------------------------
#
# **Geological Model Setup**
#
# We configure the GemPy geological model to compute gravity at observation locations.
# The model uses an implicit potential field method to interpolate geological structures
# from sparse input data (orientations and surface points).
#
# **Sigmoid slope parameter**: Controls the sharpness of lithological boundaries.
# Higher values (e.g., 100) create sharper transitions, approximating discontinuous
# interfaces between rock units. This affects gravity computation since density contrasts
# are strongest at sharp boundaries.

geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100

print(f"\n✓ Geomodel configured with {len(xy_ravel)} gravity measurement locations")

# %%
# Step 3: Compute Baseline Forward Model
# ---------------------------------------
#
# **The Forward Problem**
#
# Before inversion, we compute a forward model with initial (baseline) parameters.
# This serves multiple purposes:
#
# 1. Verify the forward model runs correctly
# 2. Assess initial fit to observations
# 3. Establish normalization parameters for the likelihood
#
# The forward gravity computation involves:
#
# 1. Interpolating geological structure from control points
# 2. Assigning densities to each rock unit
# 3. Computing gravitational attraction at each observation station
# 4. Accounting for station elevations and 3D geometry

baseline_fw_gravity_np = baseline(geo_model)
print(f"\nBaseline forward gravity:")
print(f"  Range: {baseline_fw_gravity_np.min():.1f} to {baseline_fw_gravity_np.max():.1f} µGal")
print(f"  Mean: {baseline_fw_gravity_np.mean():.1f} µGal")

# %%
# Step 4: Normalize Forward Model to Observations
# ------------------------------------------------
#
# **The Regional-Residual Problem**
#
# Gravity data contains both:
#
# - **Regional trends**: Long-wavelength signals from deep sources (e.g., Moho, crustal thickness)
# - **Residual anomalies**: Short-wavelength signals from shallow geological structures
#
# Since we're modeling shallow structures, we need to remove regional trends. The
# normalization step handles this by aligning the forward model to observations using
# a reference-based approach.
#
# **Alignment method**: "align_to_reference" finds an optimal offset that minimizes
# the difference between forward and observed gravity. The extrapolation buffer adds
# robustness by considering a range of possible offsets.
#
# This normalization is crucial because absolute gravity values depend on reference
# frames, but gravity *anomalies* (differences) are what matter for geological interpretation.

norm_params = normalize(
    baseline_fw_gravity_np,
    observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)
print(f"\nNormalization parameters:")
print(f"  Method: {norm_params['method']}")
print(f"  Parameters: {norm_params}")

# %%
# Step 5: Define Prior Distributions
# -----------------------------------
#
# **Prior Knowledge in Bayesian Inference**
#
# Priors encode our geological knowledge *before* seeing the gravity data. They serve
# multiple purposes:
#
# 1. **Regularization**: Prevent unphysical parameter values
# 2. **Information integration**: Incorporate knowledge from other data sources
# 3. **Computational efficiency**: Guide MCMC sampler to high-probability regions
#
# **Prior on Dips**
#
# We use a Normal distribution to represent our belief about layer orientations:
#
# .. math::
#
#     \\text{dip} \\sim \\mathcal{N}(10°, 10°)
#
# **Interpretation**:
#
# - **Mean (10°)**: Based on regional geology, we expect shallow-dipping layers
#   (nearly horizontal but with slight tilt)
# - **Std (10°)**: Allows for moderate uncertainty (±10° covers range of 0-30°)
# - This prior is *weakly informative*: it guides inference but can be overridden by data
#
# For a more complex model, we would also define priors on densities:
#
# .. code-block:: python
#
#     "density": dist.Normal(
#         loc=torch.tensor([2.9, 2.3]),  # plutonites, host rock
#         scale=torch.tensor(0.15)
#     ).to_event(1)
#
# Where typical crustal densities range from 2.2-3.0 g/cm³.

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
# Step 6: Define Deterministic Functions
# ---------------------------------------
#
# **Tracking Intermediate Quantities**
#
# Deterministic functions compute derived quantities during inference. These are not
# sampled but are deterministic transformations of sampled parameters and model outputs.
#
# We track:
#
# - **gravity_response_raw**: Unnormalized forward gravity (before alignment)
# - **gravity_response**: Normalized and aligned to observations (used in likelihood)
# - **mean_gravity / max_gravity**: Summary statistics for diagnostics
#
# These deterministics help with:
#
# 1. Debugging: Inspect intermediate values if inference fails
# 2. Visualization: Plot how gravity predictions evolve during sampling
# 3. Diagnostics: Check if normalized values are reasonable

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
# Step 7: Define Likelihood Function
# -----------------------------------
#
# **The Likelihood: Connecting Model to Data**
#
# The likelihood quantifies how probable the observed data is given model parameters.
# It defines the noise model:
#
# .. math::
#
#     y_i \\sim \\mathcal{N}(f_i(θ), σ^2)
#
# Where:
#
# - **y_i**: Observed gravity at station i
# - **f_i(θ)**: Forward model prediction at station i
# - **σ**: Measurement noise standard deviation
#
# **Diagonal vs Hierarchical Likelihood**
#
# We use a diagonal likelihood here (assumes independent noise across stations).
# For more robust inversion, a hierarchical likelihood can be used:
#
# .. code-block:: python
#
#     # Hierarchical: each station has its own noise level
#     likelihood_fn = generate_multigravity_likelihood_hierarchical_per_station(
#         norm_params=norm_params
#     )
#
# **Benefits of hierarchical modeling**:
#
# 1. **Automatic outlier detection**: Stations with high σ are downweighted
# 2. **Robustness**: Reduces influence of problematic measurements
# 3. **Realism**: Different stations have different measurement quality
#
# The hierarchical model estimates per-station noise:
#
# .. math::
#
#     σ_i \\sim \\text{HalfNormal}(τ) \\\\
#     y_i \\sim \\mathcal{N}(f_i(θ), σ_i^2)

likelihood_fn = MultiGravityDiagonalLikelihood(
    align_fn=align_forward_to_observed,
    norm_params=norm_params,
    sigma_value=5000.0
)
print("✓ Likelihood function created (diagonal covariance)")
print(f"  Assumed measurement noise: 5000.0 µGal")

# %%
# Step 8: Create Probabilistic Model
# -----------------------------------
#
# **Building the Full Bayesian Model**
#
# We now combine all components into a Pyro probabilistic model. This orchestrates
# the inference workflow:
#
# 1. **Sample from priors** (dips, densities)
# 2. **Update GeoModel** via `set_interp_input_fn`
# 3. **Run forward gravity computation**
# 4. **Compute deterministics** (tracked quantities)
# 5. **Evaluate likelihood** against observations
#
# **The `set_interp_input_fn` function**
#
# This function bridges Pyro's sampled parameters to GemPy's geological model:
#
# .. code-block:: python
#
#     def set_priors(samples: dict, geo_model: gp.data.GeoModel):
#         # Extract sampled dips
#         dips = samples["dips"]
#         # Update orientations in geo_model
#         geo_model.orientations_copy.dip = dips
#         # Extract sampled densities
#         densities = samples["density"]
#         # Update rock densities
#         geo_model.structural_frame.densities = densities
#         return geo_model
#
# This coupling allows MCMC to explore geological parameter space while
# computing physically-consistent gravity responses.

config = GemPyPyroConfig(
    priors=model_priors,
    set_interp_input_fn=create_orientation_modifier(key=r'dips'),
    likelihood_fn=likelihood_fn,
    pre_forward_deterministics=pre_forward_dets,
    post_forward_deterministics=post_forward_dets,
    obs_name="Gravity Measurement"
)

prob_model = GemPyPyroModelExtended(config)

print("✓ Probabilistic model created")

# %%
# Verify Model Structure
# -----------------------
#
# **Model Tracing**
#
# Before running expensive inference, we trace the model to verify its structure.
# This executes one forward pass and checks that all components work correctly.

gravity_observations_tensor = torch.tensor(observed_gravity_ugal, dtype=torch.float64)
trace = trace_pyro_model(prob_model, geo_model, gravity_observations_tensor)
print("✓ Model trace verified")

# %%
# Step 9: Prior Predictive Checks
# --------------------------------
#
# **Why Prior Predictive Checks?**
#
# Before running inference, we sample from the prior to answer critical questions:
#
# 1. **Range check**: Do prior predictions cover the observed values?
#    If not, the prior may be too restrictive or the model inadequate.
#
# 2. **Model adequacy**: Can *any* parameter combination explain the data?
#    If prior predictions are far from observations, we may need:
#    - More model complexity (additional parameters)
#    - Different physics (e.g., include magnetics, seismic)
#    - Revised priors (incorrect geological assumptions)
#
# 3. **Prior sensitivity**: How much do predictions vary under the prior?
#    High variability indicates the prior is uninformative; low variability
#    suggests the prior is too restrictive.
#
# **Expected behavior**:
#
# In this test case, we simulate 20 observations per iteration. Some forward models
# explain certain stations well but fail at others, suggesting:
#
# - The model may be oversimplified
# - Some stations could be outliers
# - Additional data types might not help without increasing model complexity
#
# Prior predictive sampling generates data *as if* we hadn't seen the observations yet.

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
# Step 10: Run MCMC Inference with NUTS
# --------------------------------------
#
# **MCMC and the NUTS Algorithm**
#
# Markov Chain Monte Carlo (MCMC) is a family of algorithms for sampling from
# probability distributions. We use the No-U-Turn Sampler (NUTS), a variant of
# Hamiltonian Monte Carlo (HMC) that:
#
# 1. **Uses gradient information**: Computes ∇log p(θ|y) to efficiently explore
#    high-dimensional parameter space
# 2. **Adapts automatically**: No manual tuning of trajectory length
# 3. **Provides diagnostics**: Divergences, tree depth, and acceptance rates
#    indicate sampling quality
#
# **How NUTS Works**
#
# NUTS simulates Hamiltonian dynamics to propose new parameter values:
#
# 1. Start at current parameter θ
# 2. Introduce momentum variable p ~ Normal(0, M)
# 3. Simulate Hamilton's equations using leapfrog integration
# 4. Build a binary tree of candidate states
# 5. Stop when trajectory makes a "U-turn" (starts returning)
# 6. Select new state from tree with probability proportional to posterior
#
# **Configuration Parameters**
#
# - **step_size (0.0001)**: Small steps for careful exploration of posterior.
#   Too large → many rejections; too small → slow mixing.
#
# - **adapt_step_size (True)**: Automatically tune step size during warmup
#   to achieve target acceptance probability.
#
# - **target_accept_prob (0.65)**: Aim for 65% acceptance rate. Lower values
#   (e.g., 0.6) → larger steps, faster exploration, more rejections. Higher values
#   (e.g., 0.9) → smaller steps, better for difficult posteriors.
#
# - **max_tree_depth (5)**: Computational budget per iteration. Each depth doubles
#   the number of gradient evaluations (2^5 = 32 max). Increase if you see
#   "reached max tree depth" warnings.
#
# - **init_strategy ('median')**: Start chains at median of prior. Alternatives:
#   'uniform' (random from prior), custom initial values.
#
# - **num_samples (200)**: Number of posterior samples after warmup. More samples
#   → better posterior approximation but slower.
#
# - **warmup_steps (200)**: Number of iterations to tune sampler. During warmup,
#   NUTS adapts step size and mass matrix. Samples are discarded.
#
# - **num_chains (1)**: Number of independent MCMC chains. Multiple chains help
#   diagnose convergence (R-hat statistic) but are computationally expensive.
#
# **Outputs**
#
# - **data.posterior**: Parameter samples (dips, densities)
# - **data.posterior_predictive**: Forward gravity predictions from posterior
# - **data.sample_stats**: Diagnostics (divergences, tree depth, acceptance rate)

print("\nRunning NUTS inference...")
print("  Warmup: 200 steps")
print("  Sampling: 200 samples")
print("  Chains: 1")

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
        num_chains=1
    ),
    plot_trace=True,
    run_posterior_predictive=True
)

print("✓ NUTS inference complete")

# %%
# Combine Prior and Posterior
# ----------------------------
#
# For comparison and diagnostics, we store both prior and posterior in the same
# ArviZ InferenceData object. This allows us to visualize how our beliefs changed
# after seeing the data.

data.extend(prior_inference_data)
print("✓ Prior and posterior combined")

# %%
# Analysis: Parameter Posterior Statistics
# -----------------------------------------
#
# **Interpreting the Posterior**
#
# The posterior distribution represents our updated beliefs about parameters after
# seeing the gravity data. Key questions:
#
# 1. **Did we learn from the data?**
#    Compare posterior std to prior std. Posterior should be narrower (more certain).
#
# 2. **Did beliefs shift?**
#    Compare posterior mean to prior mean. Shift indicates data provided information.
#
# 3. **Are parameters well-constrained?**
#    Small posterior std → data strongly constrains parameter.
#    Large posterior std → parameter remains uncertain (weak data signal).
#
# **Posterior concentration**
#
# The posterior concentrates around models that maximize the number of explained
# observations. This acts like robust regression: outliers (potentially interesting
# geological anomalies) may remain unexplained if they conflict with the majority
# of data.

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

# Quantify learning
uncertainty_reduction = (1 - posterior_dips.std() / prior_dips.std()) * 100
print(f"\nUncertainty reduction: {uncertainty_reduction:.1f}%")

# %%
# Analysis: Predictive Performance
# ---------------------------------
#
# **Posterior Predictive Checks**
#
# We compare posterior predictions to observations to assess model fit:
#
# **Residual analysis**:
#
# - **Mean residual ≈ 0**: Unbiased model (no systematic over/under-prediction)
# - **RMS residual**: Root-mean-square error quantifies overall misfit
# - **Residual patterns**: Spatial patterns indicate unmodeled geology
#
# **What if fit is poor?**
#
# If RMS residuals are large compared to expected measurement noise:
#
# 1. **Model inadequacy**: Need more complexity (additional layers, spatial
#    density variations, 3D structure)
# 2. **Wrong physics**: May need different geophysical methods
# 3. **Data issues**: Outliers, measurement errors, incorrect corrections

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
print(f"  Mean: {residuals.mean():.2f} µGal (bias)")
print(f"  RMS: {np.sqrt((residuals**2).mean()):.2f} µGal (fit quality)")
print(f"  Max absolute: {np.abs(residuals).max():.2f} µGal")

# %%
# Summary: The Inference Pipeline
# --------------------------------
#
# **Complete Workflow Recap**
#
# We've completed a full Bayesian gravity inversion:
#
# .. code-block:: none
#
#     Priors (dips, densities)
#         ↓
#     set_priors() → Update GeoModel
#         ↓
#     Forward gravity computation
#         ↓
#     Normalize & align to observations
#         ↓
#     Likelihood evaluation (diagonal)
#         ↓
#     MCMC sampling (NUTS)
#         ↓
#     Posterior (updated beliefs)
#         ↓
#     Posterior predictive (predictions with uncertainty)
#
# **Key Insights**
#
# 1. **Model limitations**: If prior predictive shows models that fit some stations
#    but not others, the forward model may be too simple. Adding more data types
#    (magnetics, seismics) won't help if the model structure is inadequate.
#
# 2. **Posterior concentration**: The posterior focuses on explaining the maximum
#    number of observations. Outliers may indicate interesting geological features
#    requiring model complexity increases.
#
# 3. **The forward model as a polynomial**: Just as y = θ₁x + θ₂ in linear regression,
#    here y = f(θ) where f is a complex geological simulation. The principles of
#    regression still apply, but f is much more sophisticated!
#
# **Next Steps**
#
# For production inversions, consider:
#
# 1. **Hierarchical likelihood**: Use per-station noise estimation for automatic
#    outlier detection and robustness
#
# 2. **Multiple chains**: Run 4+ chains to diagnose convergence with R-hat
#
# 3. **More parameters**: Invert for densities, layer thicknesses, additional
#    orientations, spatial density variations
#
# 4. **Joint inversion**: Combine gravity with magnetics, seismics, or other
#    geophysical data for better constraint
#
# 5. **Model comparison**: Use WAIC or LOO to compare alternative geological
#    hypotheses
#
# **Diagnostic checks** (not shown here but available):
#
# - Convergence diagnostics (R-hat, ESS)
# - Divergence analysis
# - Trace plots and autocorrelation
# - Prior-posterior comparison plots
# - Observed vs predicted scatter plots
# - Spatial residual maps
