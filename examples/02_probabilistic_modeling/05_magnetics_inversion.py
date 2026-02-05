"""
Bayesian Magnetic Inversion: TMI Inversion Workflow
===================================================
This tutorial demonstrates Bayesian inversion of Total Magnetic Intensity (TMI)
data using the GemPy-Probability framework. This example parallels
:ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py` but applies
the methodology to magnetic data.

**What's Different from Gravity?**

While the Bayesian framework remains the same (see the gravity tutorial for
theoretical details), magnetic inversions differ in:

1. **Forward physics**: We model magnetic susceptibility contrasts instead of
   density contrasts
2. **Earth's field dependence**: Magnetic responses depend on inclination,
   declination, and field intensity at the site
3. **Induced magnetization**: We assume rocks are magnetized by Earth's current
   magnetic field (no remanence)

**Key Parameters**

Instead of density (g/cm�), we invert for:

- **Susceptibility** (SI units): Dimensionless quantity ranging from ~0.001
  (sediments) to ~0.1 (mafic rocks)
- **Dips**: Layer orientations (same as gravity example)

**Data: Total Magnetic Intensity (TMI)**

TMI measures the magnitude of the magnetic field vector. Units are nanoTesla (nT).
Typical anomalies range from tens to thousands of nT depending on source depth
and susceptibility contrast.

For comprehensive theory on Bayesian inference, MCMC, and hierarchical modeling,
please refer to :ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py`.
"""

import os
import sys

# %%
# Import Libraries
# ----------------

import dotenv

dotenv.load_dotenv()

from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.model_one.visualization import gempy_viz, plot_many_observed_vs_forward, generate_gravity_uncertainty_plots

from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison

import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv

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
from mineye.GeoModel.model_one.model_setup import setup_geomodel
from mineye.GeoModel.model_one.probabilistic_model import normalize, set_magnetic_priors
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multimagnetic_likelihood_hierarchical_per_station
from mineye.GeoModel.model_one.visualization import probability_fields_for
import geopandas as gpd

# %%
# Define Model Extent and Resolution
# -----------------------------------

min_x = -707521
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
    project_name='magnetics_inversion',
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
# Step 1: Load Magnetic Observations
# -----------------------------------
#
# **Total Magnetic Intensity (TMI) Data**
#
# TMI measures variations in the magnitude of Earth's magnetic field caused by
# subsurface magnetic susceptibility contrasts. Magnetic rocks (e.g., plutons
# with magnetite) create positive anomalies; non-magnetic rocks (e.g.,
# sediments) create negative anomalies or lows.
#
# **Units**: nanoTesla (nT), where 1 nT = 10⁻⁹ Tesla. Earth's field is ~47,000 nT
# at this location; anomalies from geology are typically 10-1000 nT.

def read_magnetics(geophysical_dir):
    """Read magnetic data from geophysical directory."""
    magnetic_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson'))
    # Sample 20 points for this example
    sampled_magnetic = magnetic_data.sample(n=20, random_state=42)
    observed_magnetics = sampled_magnetic['MAG'].values  # in nT
    return sampled_magnetic, observed_magnetics


magnetic_data, observed_magnetics_nt = read_magnetics(geophysical_dir)
print(f"\nMagnetic observations:")
print(f"  Number of measurements: {len(observed_magnetics_nt)}")
print(f"  Range: {observed_magnetics_nt.min():.1f} to {observed_magnetics_nt.max():.1f} nT")
print(f"  Mean: {observed_magnetics_nt.mean():.1f} nT")


# %%
# Step 2: Setup Geomodel with Magnetic Configuration
# ---------------------------------------------------
#
# **Magnetic Forward Modeling**
#
# Magnetic modeling requires:
#
# 1. **Centered grid**: Measurement locations in 3D space
# 2. **Magnetic gradient tensor**: Pre-computed kernel for TMI calculation
# 3. **IGRF parameters**: International Geomagnetic Reference Field for location
#    (inclination, declination, intensity)
# 4. **Susceptibilities**: Initial values for each rock unit

def setup_magnetic_geomodel(magnetic_data, simple_geo_model):
    """Setup geomodel with magnetic measurement locations."""
    xy_ravel = np.column_stack([
            np.array(magnetic_data.geometry.x.values),
            np.array(magnetic_data.geometry.y.values),
            np.full(len(magnetic_data), 0)  # Surface elevation
    ])

    # Setup centered grid for magnetics
    gp.set_centered_grid(
        grid=simple_geo_model.grid,
        centers=xy_ravel,
        resolution=np.array([10, 10, 15]),
        radius=np.array([5000, 5000, 5000])
    )

    # Calculate magnetic gradient tensor
    from gempy_engine.modules.geophysics.magnetic_gradient import calculate_magnetic_gradient_tensor
    from gempy_engine.core.data.geophysics_input import MagneticsInput

    gradient_tensor_dict = calculate_magnetic_gradient_tensor(
        centered_grid=simple_geo_model.grid.centered_grid,
        igrf_params={
                "inclination": 60.79,  # Huelva, Spain (2025)
                "declination": 1.26,
                "intensity"  : 47258.4  # Earth's field strength in nT
        },
        compute_tmi=True,
        units_nT=True,
    )

    # Configure geophysics input
    simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
        magnetics_input=MagneticsInput(
            mag_kernel=gradient_tensor_dict['tmi_kernel'],
            susceptibilities=np.array([0.05, 0.001]),  # Initial: plutonites, host
            igrf_params={
                    "inclination": gradient_tensor_dict['inclination'],
                    "declination": gradient_tensor_dict['declination'],
                    "intensity"  : gradient_tensor_dict['intensity']
            }
        )
    )

    simple_geo_model.interpolation_options.mesh_extraction = False

    gp.set_active_grid(
        grid=simple_geo_model.grid,
        grid_type=[simple_geo_model.grid.GridTypes.CENTERED],
        reset=True
    )

    return simple_geo_model, xy_ravel


# %%
# Compute and Visualize Initial Model
# ------------------------------------

gp.compute_model(simple_geo_model)
gpv.plot_2d(simple_geo_model)

# %%
gpv.plot_3d(
    model=simple_geo_model,
    ve=5,
    image=False,
    kwargs_pyvista_bounds={
            'show_xlabels': False,
            'show_ylabels': False,
            'show_zlabels': False,
    }
)

# %%
geo_model, xy_ravel = setup_magnetic_geomodel(magnetic_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100

print(f"\n Geomodel configured with {len(xy_ravel)} magnetic measurement locations")


# %%
# Step 3: Compute Baseline Forward Model
# ---------------------------------------
#
# Compute forward magnetics with initial susceptibility values to assess
# baseline fit before inversion.

def baseline_magnetics(geo_model):
    """Compute baseline forward magnetics model."""
    sol = gp.compute_model(geo_model)
    return sol.magnetics.detach().cpu().numpy() if hasattr(sol.magnetics, 'detach') else sol.magnetics


baseline_fw_magnetics_np = baseline_magnetics(geo_model)
print(f"\nBaseline forward magnetics:")
print(f"  Range: {baseline_fw_magnetics_np.min():.1f} to {baseline_fw_magnetics_np.max():.1f} nT")
print(f"  Mean: {baseline_fw_magnetics_np.mean():.1f} nT")

# %%
# Step 4: Normalize Forward Model to Observations
# ------------------------------------------------
#
# Similar to gravity, we align the forward model to observations to handle
# regional trends and reference frame differences.

norm_params = normalize(
    baseline_fw_gravity_np=baseline_fw_magnetics_np,
    observed_gravity=observed_magnetics_nt,
    method="align_to_reference",
    extrapolation_buffer=0.3
)
print(f"\nNormalization parameters:")
print(f"  Method: {norm_params['method']}")

# %%
# Visualize baseline fit
plot_geophysics_comparison(
    forward_norm=align_forward_to_observed(baseline_fw_magnetics_np, norm_params),
    normalization_method="align_to_reference",
    observed_ugal=observed_magnetics_nt,
    xy_ravel=xy_ravel
)

# %%
# Step 5: Define Prior Distributions
# -----------------------------------
#
# **Priors for Magnetic Inversion**
#
# We define priors on:
#
# 1. **Dips**: Same as gravity example (N(10�, 10�))
# 2. **Susceptibilities**:
#    - Plutonites: 0.05 � 0.01 SI (typical for granitic rocks with magnetite)
#    - Host rock: 0.001 � 0.01 SI (low susceptibility sediments)

n_orientations = geo_model.orientations_copy.xyz.shape[0]
prior_key_dips = r'dips'
prior_key_susceptibility = r'susceptibility'

model_priors = {
        prior_key_dips          : dist.Normal(
            loc=(torch.ones(n_orientations) * 10.0),
            scale=torch.tensor(10.0, dtype=torch.float64),
            validate_args=True
        ).to_event(1),
        prior_key_susceptibility: dist.Normal(
            loc=(torch.tensor([
                    0.05,  # plutonites
                    0.001  # host rock
            ])),
            scale=torch.tensor(0.01),
        ).to_event(1)
}

print(f"\nPrior on susceptibilities:")
print(f"  Plutonites: 0.05 � 0.01 SI")
print(f"  Host rock: 0.001 � 0.01 SI")

# %%
# Step 6: Define Deterministic Functions
# ---------------------------------------

post_forward_dets = {
        "magnetic_response_raw": lambda samples, gm, sol: sol.magnetics,
        r'$\mu_{magnetics}$'    : lambda samples, gm, sol: align_forward_to_observed(
            sol.magnetics, norm_params
        ),
        "mean_magnetics"       : lambda samples, gm, sol: torch.mean(
            align_forward_to_observed(sol.magnetics, norm_params)
        ),
        "max_magnetics"        : lambda samples, gm, sol: torch.max(
            align_forward_to_observed(sol.magnetics, norm_params), 0
        ),
}


# %%
# Step 7: Define Likelihood Function
# -----------------------------------
#
# We use a hierarchical likelihood with per-station noise, analogous to the
# gravity example but with typical magnetic noise levels (~50 nT).

likelihood_fn = generate_multimagnetic_likelihood_hierarchical_per_station(
    norm_params=norm_params
)

print("Likelihood function created (hierarchical per-station)")
print(f"  Global mean noise prior: ~50.0 nT")


# %%
# Step 8: Set Interpolation Input Function
# -----------------------------------------
#
# This function updates GemPy model parameters from Pyro samples.

# %%
# Step 9: Create Probabilistic Model
# -----------------------------------

prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
    priors=model_priors,
    set_interp_input_fn=set_magnetic_priors,
    likelihood_fn=likelihood_fn,
    pre_forward_deterministics={},
    post_forward_deterministics=post_forward_dets,
    obs_name="Magnetic Measurement"
)

print("Probabilistic model created")
print("  Model components:")
print("    - Priors: dips (orientations), susceptibilities")
print("    - Forward model: GemPy geological interpolation + TMI")
print("    - Likelihood: Hierarchical per-station noise")

# %%
# Step 11: Prior Predictive Checks
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
#
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

print("\nRunning prior predictive sampling (100 samples)...")
prior_inference_data: az.InferenceData = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=torch.tensor(observed_magnetics_nt),
    n_samples=100,
    plot_trace=True
)

print("✓ Prior predictive sampling complete")

# %%
# Visualize Prior Geological Models
# ----------------------------------
#
# **Understanding gempy_viz**
#
# This function creates a 2D cross-section visualization showing:
#
# 1. **The geological model**: Interpolated lithological boundaries
# 2. **Prior uncertainty via KDE (Kernel Density Estimation)**:
#
#    - Background colored contours show probability density of boundary locations
#    - Darker/more saturated colors = higher probability
#    - Shows where the geological boundary is likely to be given our prior beliefs
#
# 3. **Representative realizations**: Individual model samples overlaid as contours
#
# **Why visualize the prior?**
#
# - Verify that prior predictions span a geologically reasonable range
# - Check if the model can produce diverse enough structures
# - Identify if priors are too restrictive (narrow KDE) or too vague (very wide KDE)
#
# **KDE interpretation**:
#
# - **Narrow, focused density**: Strong prior belief about structure location
# - **Wide, diffuse density**: High uncertainty in structure location
# - **Multiple modes**: Prior suggests multiple possible configurations

gempy_viz(
    geo_model=geo_model,
    prior_inference_data=prior_inference_data,
    n_samples=20
)

# %%
# Compare Multiple Prior Predictions to Observations
# ---------------------------------------------------
#
# **Understanding plot_many_observed_vs_forward**
#
# This diagnostic plot shows how well different prior samples fit the data.
# It helps answer: "Can ANY model from our prior explain the observations?"
#
# **Plot Structure**:
#
# - **X-axis**: Observed magnetics (sorted by value for clarity)
# - **Y-axis**: Forward-modeled magnetics from different prior samples
# - **Each colored line**: One realization from the prior distribution
# - **Red dashed line**: Perfect 1:1 agreement

plot_many_observed_vs_forward(
    forward_norm=(align_forward_to_observed(baseline_fw_magnetics_np, norm_params)),
    many_forward_norm=prior_inference_data.prior[r'$\mu_{magnetics}$'].values[0, -10:],
    observed_norm=observed_magnetics_nt
)

# %%
# Step 10: Load Pre-computed Results
# -----------------------------------
#
# To save computation time, we load results from a previous MCMC run.
# The simulation was performed with the same configuration as the gravity
# example (200 warmup steps, 200 samples, NUTS sampler).
#
# To run the inference yourself, set RUN_SIMULATION = True below.

RUN_SIMULATION = False

if RUN_SIMULATION:
    print("\nRunning NUTS inference...")

    magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)

    # Prior predictive checks
    prior_inference_data: az.InferenceData = gpp.run_predictive(
        prob_model=prob_model,
        geo_model=geo_model,
        y_obs_list=magnetic_observations_tensor,
        n_samples=100,
        plot_trace=True
    )

    # MCMC inference
    data = gpp.run_nuts_inference(
        prob_model=prob_model,
        geo_model=geo_model,
        y_obs_list=magnetic_observations_tensor,
        config=NUTSConfig(
            step_size=0.0001,
            adapt_step_size=True,
            target_accept_prob=0.65,
            max_tree_depth=5,
            init_strategy='median',
            num_samples=200,
            warmup_steps=200,
        ),
        plot_trace=True,
        run_posterior_predictive=True
    )

    data.extend(prior_inference_data)
    print(" NUTS inference complete")

else:
    from pathlib import Path
    import inspect

    # Get the directory of the current file
    current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
    data_path = current_dir / "arviz_data_magnetic_feb2026.nc"

    if not data_path.exists():
        raise FileNotFoundError(
            f"Data file not found at {data_path}. "
            f"Please run the simulation first with RUN_SIMULATION=True"
        )

    # Read the data file
    data = az.from_netcdf(str(data_path))
    print(f" Loaded inference results from {data_path}")

# %%
# **Interactive: Pre-computed Results**
#
# Running full MCMC chains can be computationally expensive. For convenience,
# pre-computed `.nc` files (ArviZ InferenceData) are available for download.
#
# To use them:
# 1. Download `arviz_data_05.nc` from the project repository.
# 2. Place it in the same directory as this script.
# 3. Set `RUN_SIMULATION = False` to load the results instantly.

# %%
# Analysis: Parameter Posterior Statistics
# -----------------------------------------

posterior_dips = data.posterior['dips'].values
print(f"\nPosterior dip statistics:")
print(f"  Shape: {posterior_dips.shape}")
print(f"  Mean: {posterior_dips.mean():.2f}�")
print(f"  Std: {posterior_dips.std():.2f}�")

posterior_suscept = data.posterior['susceptibility'].values
print(f"\nPosterior susceptibility statistics:")
print(f"  Plutonites - Mean: {posterior_suscept[:, :, 0].mean():.4f} SI")
print(f"  Plutonites - Std: {posterior_suscept[:, :, 0].std():.4f} SI")
print(f"  Host rock - Mean: {posterior_suscept[:, :, 1].mean():.4f} SI")
print(f"  Host rock - Std: {posterior_suscept[:, :, 1].std():.4f} SI")

if hasattr(data, 'prior'):
    prior_dips = data.prior['dips'].values
    uncertainty_reduction = (1 - posterior_dips.std() / prior_dips.std()) * 100
    print(f"\nUncertainty reduction in dips: {uncertainty_reduction:.1f}%")

# %%
# Analysis: Predictive Performance
# ---------------------------------

posterior_magnetics = data.posterior_predictive['Magnetic Measurement'].values
print(f"\nPosterior predictive magnetics:")
print(f"  Mean: {posterior_magnetics.mean():.1f} nT")
print(f"  Std: {posterior_magnetics.std():.1f} nT")

residuals = posterior_magnetics.mean(axis=(0, 1)) - observed_magnetics_nt
print(f"\nResiduals (posterior mean - observed):")
print(f"  Mean: {residuals.mean():.2f} nT (bias)")
print(f"  RMS: {np.sqrt((residuals ** 2).mean()):.2f} nT (fit quality)")

# %%
# Posterior Predictive: Model Fit Analysis
# -----------------------------------------
#
# **Understanding Posterior Convergence**
#
# In the plot_many_observed_vs_forward visualization, we observe that no single model
# can perfectly explain all observations simultaneously. This reveals an important aspect
# of Bayesian inference: the posterior distribution identifies models that best balance
# the fit across all measurement stations.
#
# **What the inference is doing**:
#
# The MCMC sampler finds parameter combinations that bring the maximum number of
# observations close to their measured values. Rather than perfectly fitting a subset
# of stations, the posterior concentrates on models that provide reasonable explanations
# across the entire dataset. This "compromise" solution is mathematically optimal under
# our likelihood model.
#
# **Geometric constraints from data**:
#
# Even with this distributed fit, the magnetic data significantly constrains the possible
# geological configurations. As visualized in gempy_viz, the posterior distribution shows
# much tighter bounds on structural geometry than the prior. The inference has successfully
# ruled out large regions of parameter space that are inconsistent with observations.
#
# **Implications**:
#
# - Stations that remain poorly fit may indicate localized geological complexity not
#   captured by our current model structure
# - The spread in posterior predictions quantifies irreducible uncertainty given model
#   assumptions
# - Future model refinements (additional layers, spatially-varying susceptibilities) could
#   improve station-by-station fit while maintaining these geometric constraints

plot_many_observed_vs_forward(
    forward_norm=(align_forward_to_observed(baseline_fw_magnetics_np, norm_params)),
    many_forward_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].values[0, -20:],
    observed_norm=observed_magnetics_nt
)

# %%
# Density plots comparing prior and posterior distributions

az.plot_density(
    data=[data, data.prior],
    var_names=["dips", "susceptibility"],
    filter_vars="like",
    hdi_prob=0.9999,
    shade=.2,
    data_labels=["Posterior", "Prior"],
    colors=[default_red, default_blue],
)

# %%
# Geological Model Uncertainty Visualization
# -------------------------------------------

gempy_viz(
    geo_model=geo_model,
    prior_inference_data=data,
    n_samples=20
)

# %%
# Spatial Comparison: Prior Predictions

plot_geophysics_comparison(
    forward_norm=data.prior[r'$\mu_{magnetics}$'].mean(axis=1),
    normalization_method='align_to_reference',
    observed_ugal=observed_magnetics_nt,
    xy_ravel=xy_ravel
)

# %%
# Spatial Comparison: Posterior Predictions

plot_geophysics_comparison(
    forward_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].mean(axis=1),
    normalization_method='align_to_reference',
    observed_ugal=observed_magnetics_nt,
    xy_ravel=xy_ravel
)

# %%
# Uncertainty Quantification: Prior

gravity_samples_norm, unit_label = generate_gravity_uncertainty_plots(
    gravity_samples_norm=data.prior[r'$\mu_{magnetics}$'].values[0, :],
    observed_gravity_ugal=observed_magnetics_nt,
    xy_ravel=xy_ravel
)

# %%
# Uncertainty Quantification: Posterior

if hasattr(data, 'posterior_predictive') and r'$\mu_{magnetics}$' in data.posterior_predictive:
    gravity_samples_norm, unit_label = generate_gravity_uncertainty_plots(
        gravity_samples_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].values[0, :],
        observed_gravity_ugal=observed_magnetics_nt,
        xy_ravel=xy_ravel
    )

# %%
# **Sigma Analysis: Outlier Detection**
#
# Hierarchical modeling allows us to identify stations with unusually high noise,
# which may indicate instrument problems, terrain correction errors, or localized
# geological complexity not captured by the model.

if "sigma_stations" in data.posterior_predictive:
    posterior_sigmas = data.posterior_predictive["sigma_stations"].values
    station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
    sigma_global_mean = station_noise_mean.mean()

    # Identify stations with noise > 2 standard deviations above mean
    problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
    print(f"\nPotential outlier stations identified: {problematic}")

    # Plot sigma distribution
    az.plot_density(
        data=[data, data.prior],
        var_names=["sigma_stations"],
        filter_vars="like",
        hdi_prob=0.9999,
        shade=.2,
        data_labels=["Posterior", "Prior"],
        colors=[default_red, default_blue],
    )
    plt.title("Per-Station Noise Distribution (Sigma)")
    plt.show()

# %%
# **Probability Density Fields and Information Entropy**
#
# To visualize the spatial uncertainty of the geological structure, we compute
# probability density fields and information entropy.


# Resetting the model
geo_model = gp.create_geomodel(
    project_name='gravity_inversion',
    extent=extent,
    refinement=refinement,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)
# Prior Probability Fields
print("\nComputing prior probability fields...")
topography_path = paths.get_topography_path()
probability_fields_for(
    geo_model=geo_model,
    inference_data=data.prior,
    topography_path=topography_path
)

# Posterior Probability Fields
if hasattr(data, 'posterior'):
    print("\nComputing posterior probability fields...")
    probability_fields_for(
        geo_model=geo_model,
        inference_data=data.posterior,
        topography_path=topography_path
    )

# %%
# **3D Entropy Visualization**
#
# We can also visualize uncertainty in 3D by injecting the entropy field back
# into the GemPy solutions object.

# Note: The 3D visualization is already handled inside probability_fields_for()
# by injecting the entropy field and calling gpv.plot_3d.

# %%
# Summary
# -------
#
# **Key Takeaways**
#
# 1. **Same Framework, Different Physics**: The Bayesian workflow for magnetic
#    inversion is identical to gravity (see example 04), but we model
#    susceptibility instead of density.
#
# 2. **TMI is Orientation-Dependent**: Magnetic anomalies depend on Earth's field
#    direction (inclination, declination), making interpretation more complex
#    than gravity.
#
# 3. **Hierarchical Noise Modeling**: Per-station noise inference automatically
#    handles outliers and varying data quality.
#
# 4. **Joint Interpretation**: For best results, consider joint gravity-magnetic
#    inversion to simultaneously constrain both density and susceptibility.
#
# For detailed explanations of:
# - Bayesian theory and MCMC
# - Hierarchical modeling
# - Prior predictive checks
# - Posterior diagnostics
#
# Please refer to :ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py`

# sphinx_gallery_thumbnail_number = 11
