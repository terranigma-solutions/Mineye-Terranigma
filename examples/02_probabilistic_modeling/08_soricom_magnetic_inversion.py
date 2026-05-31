"""
Bayesian Magnetic Inversion: Soricom Chromite Lens
===================================================
This tutorial demonstrates Bayesian inversion of Total Magnetic Intensity (TMI)
data using the Soricom fault structural model. It extends the magnetic inversion
methodology to a more complex geological setting with a fault and a chromite
lens hosted in ultramafic rocks.

**What Makes the Soricom Model Different?**

Unlike the simple Tharsis plutonite model (:ref:`sphx_glr_02_probabilistic_modeling_05_magnetics_inversion.py`),
the Soricom model features:

1. **A fault**: The Main_Fault truncates all formations (fault-first structural
   frame ordering)
2. **A chromite lens**: A thin, high-susceptibility target embedded in host rock
3. **Smaller extent**: ~500 m × 350 m area at UTM zone 34N coordinates
4. **Higher resolution**: Octree refinement level 5 on a much smaller domain

**Data and Preprocessing**

TMI data is extracted from a merged B1B2 raster at the Soricom prospect.
The raw measurements are in **nanoTesla (nT)** and include Earth's background
magnetic field (~47,500 nT at this location). Since the GemPy forward model
computes **TMI anomalies** (deviations from the background field), we subtract
the IGRF (International Geomagnetic Reference Field) intensity before inversion:

.. math::

    TMI_{\text{anomaly}} = TMI_{\text{measured}} - IGRF_{\text{intensity}}


For comprehensive theory on Bayesian inversion, MCMC, and hierarchical
likelihoods, see :ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py`.
"""

import os
import sys

# %%
# Import Libraries
# ----------------

import dotenv

dotenv.load_dotenv()

from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from gempy_probability.modules.plot.plot_gempy import plot_gempy
from mineye.GeoModel.model_one.visualization import (
    plot_many_observed_vs_forward,
    generate_gravity_uncertainty_plots,
    probability_fields_for
)
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison

import numpy as np
import matplotlib.pyplot as plt

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
# Import Soricom-specific helpers
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import (
    generate_multimagnetic_likelihood_fixed_std,
)
from gempy_engine.modules.geophysics.magnetic_gradient import calculate_magnetic_gradient_tensor
from gempy_engine.core.data.geophysics_input import MagneticsInput

# --- Helper: Z-position update function (used by gempy_viz and probability fields) ---
_original_z_coords = None


def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: np.ndarray, sample_idx: int):
    global _original_z_coords
    if _original_z_coords is None:
        _original_z_coords = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)

    scale_z = geo_model.input_transform.scale[2]
    if hasattr(sample_value, 'detach'):
        sample_value = sample_value.detach().cpu().numpy()
    shifts_m = sample_value / scale_z

    new_z = _original_z_coords.copy()
    # Point 0: Main_Fault (no change)
    # Points 1-12: host_rock -> shifted by shifts_m[0] (layer wide)
    new_z[1:13] = _original_z_coords[1:13] + shifts_m[0]
    # Points 13-21: chromite lense -> shifted by shifts_m[1:] (independent)
    new_z[13:22] = _original_z_coords[13:22] + shifts_m[1:]

    gp.modify_surface_points(
        geo_model=geo_model,
        Z=new_z
    )


def gempy_viz(geo_model: gp.data.GeoModel, prior_inference_data: az.InferenceData,
              n_samples=20, ve=3):
    gp.set_active_grid(
        grid=geo_model.grid,
        grid_type=[geo_model.grid.GridTypes.OCTREE],
        reset=True
    )
    geo_model.geophysics_input = None
    gp.compute_model(gempy_model=geo_model)

    p2d = gpv.plot_2d(
        model=geo_model,
        show_topography=False,
        legend=False,
        show_lith=False,
        show_data=False,
        show=False,
        ve=ve
    )

    original_z = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)
    global _original_z_coords
    _original_z_coords = None

    plot_gempy(
        geo_model=geo_model,
        n_samples=n_samples,
        samples=(prior_inference_data.prior['surface_points_z'].values[0, :]),
        update_model_fn=_update_model_for_plotting,
        gempy_plot=p2d,
    )

    if hasattr(prior_inference_data, 'posterior'):
        gp.modify_surface_points(geo_model=geo_model, Z=original_z)
        _original_z_coords = None
        n_surfaces = len(geo_model.structural_frame.elements_colors_contacts)
        plot_gempy(
            geo_model=geo_model,
            n_samples=n_samples,
            samples=(prior_inference_data.posterior['surface_points_z'].values[0, :]),
            update_model_fn=_update_model_for_plotting,
            gempy_plot=p2d,
            contour_colors=[default_red] * n_surfaces,
        )

    return p2d


# %%
# Step 1: Load Preprocessed Magnetic Observations
# ------------------------------------------------
#
# **Preprocessing Applied**
#
# The raw TMI raster values (~47,000 nT) include Earth's ambient magnetic field.
# To match the forward model (which computes anomalies), we subtract the IGRF
# intensity of 47,500 nT. The resulting anomalies range from -392 to -42 nT,
# consistent with the forward model output (~-550 to -475 nT).


def read_magnetics():
    """Read preprocessed magnetic data from extracted numpy arrays."""
    try:
        this_file = __file__
    except NameError:
        this_file = sys.argv[0]
    current_dir = os.path.dirname(os.path.abspath(this_file))
    xyz_path = os.path.join(current_dir, 'soricom_magnetic_xyz_adaptive.npy')
    mag_path = os.path.join(current_dir, 'soricom_magnetic_values_adaptive.npy')
    xyz = np.load(xyz_path)
    mag = np.load(mag_path)

    # Subsample 20 points for inversion (random subset, seed fixed)
    rng = np.random.default_rng(42)
    idx = rng.choice(len(xyz), size=min(20, len(xyz)), replace=False)
    sampled_xyz = xyz[idx]

    # IGRF: forward model outputs TMI anomalies, so subtract IGRF from raw TMI
    igrf_intensity_nT = 47_500.0
    observed_magnetics = mag[idx] - igrf_intensity_nT
    return sampled_xyz, observed_magnetics


magnetic_xyz, observed_magnetics_nt = read_magnetics()
print(f"\nMagnetic observations:")
print(f"  Number of measurements: {len(observed_magnetics_nt)}")
print(f"  Range: {observed_magnetics_nt.min():.1f} to {observed_magnetics_nt.max():.1f} nT")
print(f"  Mean: {observed_magnetics_nt.mean():.1f} nT")

# --- Prior keys ---
prior_key_surface_points_z = r'surface_points_z'
prior_key_susceptibility = r'susceptibility'

# %%
# **Z-Position Prior Setter**
#
# The ``_set_magnetic_priors`` function maps sampled parameter values onto the
# GemPy model. It handles two parameter groups:
#
# 1. ``surface_points_z``: applies vertical shifts to surface points (meter units)
# 2. ``susceptibility``: updates the magnetic susceptibility of each unit
#
# The host rock's 12 surface points share one shift (layer-wide movement),
# while the chromite lens's 9 points each have independent shifts.


def _set_magnetic_priors(samples: dict, geo_model: gp.data.GeoModel):
    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame
    interp_input = interpolation_input_from_structural_frame(geo_model)

    if prior_key_surface_points_z in samples:
        shifts = samples[prior_key_surface_points_z]
        coords = interp_input.surface_points.sp_coords.clone()
        # coords is on CPU (gempy_engine expects CPU tensors); bring shifts to CPU too
        shifts_cpu = shifts.to(coords.device)

        # Index 1:13 = 12 host_rock points → shifted by shifts[0] (layer wide)
        coords[1:13, 2] = coords[1:13, 2] + shifts_cpu[0]
        # Index 13:22 = 9 chromite lense points → shifted by shifts[1:] (independent)
        coords[13:22, 2] = coords[13:22, 2] + shifts_cpu[1:]

        interp_input.surface_points.sp_coords = coords

    if prior_key_susceptibility in samples:
        susceptibilities = samples[prior_key_susceptibility]
        if geo_model.geophysics_input and geo_model.geophysics_input.magnetics_input:
            geo_model.geophysics_input.magnetics_input.susceptibilities = susceptibilities

    return interp_input


# %%
# Step 2: Setup Geomodel with Soricom Fault Configuration
# --------------------------------------------------------

geo_model = gp.create_geomodel(
    project_name=SoricomModelConfig.PROJECT_NAME,
    extent=SoricomModelConfig.EXTENT,
    refinement=SoricomModelConfig.REFINEMENT,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=paths.get_soricom_orientations(),
        path_to_surface_points=paths.get_soricom_formation_points(),
    ),
)

geo_model.grid = geo_model.grid.init_octree_grid(
    extent=SoricomModelConfig.EXTENT,
    octree_levels=SoricomModelConfig.REFINEMENT,
)
geo_model.interpolation_options.number_octree_levels_surface = 2

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object=SoricomModelConfig.SURFACE_MAPPING,
)

geo_model.structural_frame.structural_groups[
    SoricomModelConfig.FAULT_GROUP_INDEX
].structural_relation = gp.data.StackRelationType.FAULT
geo_model.structural_frame.fault_relations = SoricomModelConfig.FAULT_RELATIONS_MATRIX

# %%
# Compute initial model to see geometry
gp.compute_model(geo_model)
gpv.plot_2d(geo_model)

# %%
gpv.plot_3d(
    model=geo_model,
    ve=5,
    image=False,
    kwargs_pyvista_bounds={
        'show_xlabels': False,
        'show_ylabels': False,
        'show_zlabels': False,
    },
)

# %%
# Step 3: Setup Magnetic Forward Model
# -------------------------------------
#
# **Magnetic Configuration**
#
# We configure centered grids at each measurement location and pre-compute the
# TMI gradient kernel (the linear operator that maps susceptibility to TMI).
#
# **Susceptibility values** (initial guess, in SI units):
#
# - ``Main_Fault``: 0.0001 (near-zero)
# - ``host_rock``: 0.0001 (ultramafic host, near-zero for inversion start)
# - ``chromite_lense``: 0.5 (high-susceptibility target)
# - ``basement``: 0.0001 (low-susceptibility country rock)
#
# **IGRF parameters** at Soricom (WGS 84 UTM zone 34N):
# - Inclination: 57.0°
# - Declination: 4.0°
# - Intensity: 47,500 nT

gp.set_centered_grid(
    grid=geo_model.grid,
    centers=magnetic_xyz,
    resolution=np.array([10, 10, 15]),
    radius=np.array([2000, 5000, 2000]),
)

gradient_tensor_dict = calculate_magnetic_gradient_tensor(
    centered_grid=geo_model.grid.centered_grid,
    igrf_params={
        "inclination": 57.0,
        "declination": 4.0,
        "intensity": 47500.0,
    },
    compute_tmi=True,
    units_nT=True,
)

geo_model.geophysics_input = gp.data.GeophysicsInput(
    magnetics_input=MagneticsInput(
        mag_kernel=gradient_tensor_dict['tmi_kernel'],
        susceptibilities=np.array([0.0001, 0.0001, 0.5, 0.0001]),
        igrf_params={
            "inclination": gradient_tensor_dict['inclination'],
            "declination": gradient_tensor_dict['declination'],
            "intensity": gradient_tensor_dict['intensity'],
        },
    ),
)

geo_model.interpolation_options.mesh_extraction = False

gp.set_active_grid(
    grid=geo_model.grid,
    grid_type=[geo_model.grid.GridTypes.CENTERED],
    reset=True,
)
gp.compute_model(geo_model)

geo_model.interpolation_options.sigmoid_slope = 100
print(f"\n✓ Geomodel configured with {len(magnetic_xyz)} magnetic measurement locations")

# %%
# Step 4: Compute Baseline Forward Model
# ---------------------------------------

sol = gp.compute_model(geo_model)
baseline_fw_magnetics_np = (
    sol.magnetics.detach().cpu().numpy()
    if hasattr(sol.magnetics, 'detach')
    else sol.magnetics
)
print(f"\nBaseline forward magnetics:")
print(f"  Range: {baseline_fw_magnetics_np.min():.1f} to {baseline_fw_magnetics_np.max():.1f} nT")
print(f"  Mean: {baseline_fw_magnetics_np.mean():.1f} nT")

# %%
# Step 5: Normalize Forward Model to Observations
# ------------------------------------------------
#
# We align the forward model anomalies to the observed anomaly scale using
# ``align_to_reference``, matching mean and standard deviation.

norm_params = normalize(
    baseline_fw_gravity_np=baseline_fw_magnetics_np,
    observed_gravity=observed_magnetics_nt,
    method="align_to_reference",
    extrapolation_buffer=0.3,
)
print(f"\nNormalization parameters:")
print(f"  Method: {norm_params['method']}")

# %%
# Visualize baseline fit
plot_geophysics_comparison(
    forward_norm=align_forward_to_observed(baseline_fw_magnetics_np, norm_params),
    normalization_method="align_to_reference",
    observed_ugal=observed_magnetics_nt,
    xy_ravel=magnetic_xyz,
    unit_label='nT',
)

# %%
# Step 6: Define Prior Distributions
# -----------------------------------
#
# **Two Parameter Groups**
#
# 1. **Surface point Z shifts** (:math:`\Delta z`):
#    - ``host_rock`` (12 surface points): one shared shift — layer-wide vertical movement
#    - ``chromite_lense`` (9 surface points): 9 independent shifts per point
#    - Main_Fault (1 point): fixed (no shift)
#    - Total: **10 parameters** (1 layer-wide + 9 per-point)
#    - Prior: :math:`\mathcal{N}(0, 15 \, \text{m})` in scaled coordinates
#
# 2. **Susceptibility** (SI units, LogNormal):
#    - :math:`\log(\kappa) \sim \mathcal{N}(\mu, \sigma)` for each unit
#    - Corresponds to ~0.0001 SI (fault/basement), ~0.0001 SI (host), ~0.5 SI (chromite)

model_priors = {
    prior_key_surface_points_z: dist.Normal(
        loc=torch.tensor([0.0] * 10, dtype=torch.float64),
        scale=torch.tensor([15.0] * 10, dtype=torch.float64) * geo_model.input_transform.scale[2],
        validate_args=True,
    ).to_event(1),
    prior_key_susceptibility: dist.LogNormal(
        loc=torch.tensor([-9.21, -9.21, -0.69, -9.21], dtype=torch.float64),
        scale=torch.tensor([0.1, 0.1, 0.2, 0.1], dtype=torch.float64),
    ).to_event(1),
}

# %%
# Step 7: Define Deterministic Functions
# ---------------------------------------

post_forward_dets = {
    "magnetic_response_raw": lambda samples, gm, sol: sol.magnetics,
    r'$\mu_{magnetics}$': lambda samples, gm, sol: align_forward_to_observed(
        sol.magnetics, norm_params,
    ),
    "mean_magnetics": lambda samples, gm, sol: torch.mean(
        align_forward_to_observed(sol.magnetics, norm_params),
    ),
    "max_magnetics": lambda samples, gm, sol: torch.max(
        align_forward_to_observed(sol.magnetics, norm_params), 0,
    ),
}

# %%
# Step 8: Define Likelihood Function
# -----------------------------------
#
# Fixed-standard-deviation likelihood with :math:`\sigma = 150` nT. This is
# approximately 2× the baseline residual standard deviation, giving NUTS room
# to explore parameter space while still penalizing poor fits.

likelihood_fn = generate_multimagnetic_likelihood_fixed_std(
    norm_params=norm_params,
    sigma_value=150.0,
)

print("Likelihood function created (fixed std = 150 nT)")

# %%
# Step 9: Create Probabilistic Model
# -----------------------------------

prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model(
    priors=model_priors,
    set_interp_input_fn=_set_magnetic_priors,
    likelihood_fn=likelihood_fn,
    obs_name="Magnetic Measurement (Soricom)",
)

print("Probabilistic model created")
print("  Model components:")
print("    - Priors: surface points Z shifts, susceptibility")
print("    - Forward model: GemPy geological interpolation + TMI")
print("    - Likelihood: Fixed standard deviation (150 nT)")
print("    - Deterministics: magnetic_response, mean_magnetics, max_magnetics")

# %%
# Step 10: Prior Predictive Checks
# --------------------------------

print("\nRunning prior predictive sampling (100 samples)...")
prior_inference_data: az.InferenceData = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=torch.tensor(observed_magnetics_nt),
    n_samples=100,
    plot_trace=True,
)

print("✓ Prior predictive sampling complete")

# %%
# Visualize Prior Geological Models
# ----------------------------------

gempy_viz(
    geo_model=geo_model,
    prior_inference_data=prior_inference_data,
    n_samples=20,
)

# %%
# Compare Multiple Prior Predictions to Observations

plot_many_observed_vs_forward(
    forward_norm=align_forward_to_observed(baseline_fw_magnetics_np, norm_params),
    many_forward_norm=prior_inference_data.prior[r'$\mu_{magnetics}$'].values[0, -10:],
    observed_norm=observed_magnetics_nt,
    unit_label='nT',
)

# %%
# Step 11: Load Pre-computed Results
# ------------------------------------
#
# Full MCMC inference takes hours. Load pre-computed results from a previous
# run with 200 warmup steps, 200 samples, and NUTS sampler.
#
# To run inference yourself, set ``RUN_SIMULATION = True`` below.

RUN_SIMULATION = False

if RUN_SIMULATION:
    print("\nRunning NUTS inference...")

    magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)

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
        run_posterior_predictive=True,
    )

    data.extend(prior_inference_data)
    print("✓ NUTS inference complete")

else:
    from pathlib import Path
    import inspect

    current_dir = Path(inspect.getfile(inspect.currentframe())).parent.resolve()
    data_path = current_dir / "arviz_data_magnetic_soricom_z.nc"

    if not data_path.exists():
        raise FileNotFoundError(
            f"Data file not found at {data_path}. "
            f"Please run the simulation first with RUN_SIMULATION=True"
        )

    data = az.from_netcdf(str(data_path))
    print(f"✓ Loaded inference results from {data_path}")

# %%
# Analysis: Parameter Posterior Statistics
# -----------------------------------------

posterior_z = data.posterior['surface_points_z'].values
print(f"\nPosterior Z-shift statistics (host rock layer-wide):")
print(f"  Mean: {posterior_z[:, :, 0].mean():.2f} m")
print(f"  Std:  {posterior_z[:, :, 0].std():.2f} m")
print(f"\nPosterior Z-shift statistics (chromite per-point, first 3):")
for i in range(1, min(4, posterior_z.shape[-1])):
    print(f"  Point {i}: mean={posterior_z[:, :, i].mean():.2f} m, std={posterior_z[:, :, i].std():.2f} m")

posterior_suscept = data.posterior['susceptibility'].values
print(f"\nPosterior susceptibility statistics:")
print(f"  Fault       - Mean: {posterior_suscept[:, :, 0].mean():.4f} SI")
print(f"  Host rock   - Mean: {posterior_suscept[:, :, 1].mean():.4f} SI")
print(f"  Chromite    - Mean: {posterior_suscept[:, :, 2].mean():.4f} SI")
print(f"  Basement    - Mean: {posterior_suscept[:, :, 3].mean():.4f} SI")

if hasattr(data, 'prior'):
    prior_z = data.prior['surface_points_z'].values
    uncertainty_reduction = (1 - posterior_z.std() / prior_z.std()) * 100
    print(f"\nUncertainty reduction in Z-positions: {uncertainty_reduction:.1f}%")
    prior_suscept = data.prior['susceptibility'].values
    suscept_uncertainty_reduction = (
        1 - posterior_suscept.std() / prior_suscept.std()
    ) * 100
    print(f"Uncertainty reduction in susceptibility: {suscept_uncertainty_reduction:.1f}%")

# %%
# Analysis: Predictive Performance
# ---------------------------------

posterior_magnetics = data.posterior_predictive['Magnetic Measurement (Soricom)'].values
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

plot_many_observed_vs_forward(
    forward_norm=align_forward_to_observed(baseline_fw_magnetics_np, norm_params),
    many_forward_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].values[0, -20:],
    observed_norm=observed_magnetics_nt,
    unit_label='nT',
)

# %%
# Density plots comparing prior and posterior distributions

az.plot_density(
    data=[data, data.prior],
    var_names=["surface_points_z", "susceptibility"],
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
    n_samples=20,
)

# %%
# Spatial Comparison: Prior Predictions

plot_geophysics_comparison(
    forward_norm=data.prior[r'$\mu_{magnetics}$'].mean(axis=1),
    normalization_method='align_to_reference',
    observed_ugal=observed_magnetics_nt,
    xy_ravel=magnetic_xyz,
    unit_label='nT',
)

# %%
# Spatial Comparison: Posterior Predictions

plot_geophysics_comparison(
    forward_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].mean(axis=1),
    normalization_method='align_to_reference',
    observed_ugal=observed_magnetics_nt,
    xy_ravel=magnetic_xyz,
    unit_label='nT',
)

# %%
# Uncertainty Quantification: Prior
#
# The ``generate_gravity_uncertainty_plots`` function works with any geophysical
# data type — we pass ``unit_label`` and ``response_label`` to override the
# default gravity labels.

magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
    gravity_samples_norm=data.prior[r'$\mu_{magnetics}$'].values[0, :],
    observed_gravity_ugal=observed_magnetics_nt,
    xy_ravel=magnetic_xyz,
    unit_label='nT',
    response_label='Magnetics',
    title_prefix='Magnetic',
)

# %%
# Uncertainty Quantification: Posterior

if hasattr(data, 'posterior_predictive') and r'$\mu_{magnetics}$' in data.posterior_predictive:
    magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
        gravity_samples_norm=data.posterior_predictive[r'$\mu_{magnetics}$'].values[0, :],
        observed_gravity_ugal=observed_magnetics_nt,
        xy_ravel=magnetic_xyz,
        unit_label='nT',
        response_label='Magnetics',
        title_prefix='Magnetic',
    )

# %%
# **Residual Analysis**
#
# With a fixed-standard-deviation likelihood, we examine the residual
# distribution directly. A well-calibrated model should produce residuals
# that are roughly normally distributed with mean near zero.

posterior_predictive_values = data.posterior_predictive[
    'Magnetic Measurement (Soricom)'
].values
residuals_all = posterior_predictive_values - observed_magnetics_nt[np.newaxis, np.newaxis, :]
residuals_flat = residuals_all.reshape(-1)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
axes[0].hist(residuals_flat, bins=50, density=True, alpha=0.7, color=default_blue)
axes[0].set_xlabel('Residual (nT)')
axes[0].set_ylabel('Density')
axes[0].set_title(f'Residual Distribution\nMean = {residuals_flat.mean():.1f} nT, Std = {residuals_flat.std():.1f} nT')

station_rmse = np.sqrt((posterior_predictive_values.mean(axis=(0, 1)) - observed_magnetics_nt) ** 2)
axes[1].bar(range(len(station_rmse)), station_rmse, color=default_red, alpha=0.7)
axes[1].set_xlabel('Station index')
axes[1].set_ylabel('RMSE (nT)')
axes[1].set_title('Per-Station RMSE')
axes[1].axhline(y=150, color='gray', linestyle='--', label='Likelihood σ = 150 nT')
axes[1].legend()
plt.tight_layout()
plt.show()

# %%
# **Probability Density Fields and Information Entropy**
#
# To visualize the spatial uncertainty of the geological structure, we compute
# probability density fields and information entropy.
#
# * **Probability Density Field**: Shows the probability of each geological unit
#   existing at any given location.
# * **Information Entropy**: Quantifies the total uncertainty. High entropy means
#   high uncertainty about which unit is present.
#
# Note: The 3D visualization is already handled inside probability_fields_for()
# by generating a multi-panel PyVista paper-quality figure.

# Recreate the model to reset grid state
geo_model = gp.create_geomodel(
    project_name=SoricomModelConfig.PROJECT_NAME,
    extent=SoricomModelConfig.EXTENT,
    refinement=SoricomModelConfig.REFINEMENT,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=paths.get_soricom_orientations(),
        path_to_surface_points=paths.get_soricom_formation_points(),
    ),
)

geo_model.grid = geo_model.grid.init_octree_grid(
    extent=SoricomModelConfig.EXTENT,
    octree_levels=SoricomModelConfig.REFINEMENT,
)
geo_model.interpolation_options.number_octree_levels_surface = 2

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object=SoricomModelConfig.SURFACE_MAPPING,
)

geo_model.structural_frame.structural_groups[
    SoricomModelConfig.FAULT_GROUP_INDEX
].structural_relation = gp.data.StackRelationType.FAULT
geo_model.structural_frame.fault_relations = SoricomModelConfig.FAULT_RELATIONS_MATRIX

# Prior Probability Fields
print("\nComputing prior probability fields...")
topography_path = paths.get_soricom_dem_path()

# Cache Z for restore between prior and posterior
original_z_prior = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)
_original_z_coords = None

probability_fields_for(
    geo_model=geo_model,
    inference_data=data.prior,
    topography_path=topography_path,
    var_name=prior_key_surface_points_z,
    update_model_fn=_update_model_for_plotting,
    ve=1,
)

# Posterior Probability Fields
if hasattr(data, 'posterior'):
    # Restore model Z so posterior OnlineProbability init uses baseline lithology
    gp.modify_surface_points(geo_model=geo_model, Z=original_z_prior)
    _original_z_coords = None
    print("\nComputing posterior probability fields...")
    probability_fields_for(
        geo_model=geo_model,
        inference_data=data.posterior,
        topography_path=topography_path,
        var_name=prior_key_surface_points_z,
        update_model_fn=_update_model_for_plotting,
        ve=1,
    )

# %%
# Summary
# -------
#
# **Key Takeaways**
#
# 1. **IGRF subtraction is essential**: Raw TMI values include Earth's background
#    field (~47,500 nT). The forward model computes anomalies, so the data must
#    match. Without this step, the likelihood is numerically flat and NUTS cannot
#    move.
#
# 2. **Z-position priors are geometrically intuitive**: Unlike dip priors (which
#    rotate orientation vectors), shifting surface point Z-coordinates directly
#    controls the depth of geological boundaries. This matches how geologists
#    think about uncertainty — "the contact could be ±15 m deeper."
#
# 3. **Layer-wide vs per-point shifting**: The host rock shares one Z-shift
#    across all 12 surface points (rigid-body movement). The chromite lens uses
#    independent per-point shifts (allowing lens shape variation).
#
# 4. **Fixed-std likelihood works for well-characterized noise**: A single
#    :math:`\sigma = 150` nT (approximately 2× the baseline residual std)
#    lets NUTS explore parameters freely while producing good data fit.
#
# 5. **LogNormal susceptibility keeps values positive**: Susceptibility must be
#    strictly positive (SI convention). The LogNormal prior enforces this
#    naturally, unlike a Normal prior that could sample negative values.
#
# **Comparison with Dip-Prior Approach** (see :ref:`gravity_inversion`)
#
# - **Z-priors**: Shift surface point coordinates → directly controls boundary
#   geometry → simpler parameterization with clear physical units (meters)
# - **Dip priors**: Modify orientation pole vectors → affects gradient field
#   globally → useful when orientation uncertainty dominates
# - Both approaches can be combined, but Z-priors are preferred for magnetic
#   inversion where magnetic contrast is the primary signal and exact interface
#   positions determine the anomaly shape.
#
# For theoretical background, see:
# - :ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py`
# - :ref:`sphx_glr_02_probabilistic_modeling_05_magnetics_inversion.py`

# sphinx_gallery_thumbnail_number = 4
