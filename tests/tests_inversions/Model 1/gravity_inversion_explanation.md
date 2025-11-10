# Gravity Inversion Tutorial: A Complete Workflow

This tutorial walks through a Bayesian gravity inversion workflow using the `test_gravity_inversion.py` test case. We'll explore how to invert gravity observations to estimate geological parameters like rock densities and layer orientations.

## Table of Contents
1. [Overview](#overview)
2. [Workflow Steps](#workflow-steps)
3. [The Complete Pipeline](#the-complete-pipeline)
4. [Analysis and Diagnostics](#analysis-and-diagnostics)
5. [Key Insights](#key-insights)

---

## Overview

**What are we doing?**
We're performing Bayesian inference on a geological model using gravity measurements. The goal is to estimate:
- **Dips**: Orientation angles of geological layers (degrees)
- **Densities**: Rock densities for different lithological units (g/cm³)

**Why Bayesian?**
Bayesian inference gives us:
- Posterior probability distributions (not just point estimates)
- Uncertainty quantification
- The ability to incorporate prior knowledge
- Detection of outliers and model misspecification

---

## Workflow Steps

### Step 1: Read Gravity Data
```python
gravity_data, observed_gravity_ugal = read_gravity(geophysical_dir)
```
**Location**: `test_gravity_inversion.py:96`

This loads gravity measurements from field instruments. The data includes:
- **gravity_data**: Station locations (x, y, z coordinates)
- **observed_gravity_ugal**: Measured gravity values in microGals (µGal)

---

### Step 2: Setup Initial GeoModel
```python
geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
geo_model.interpolation_options.sigmoid_slope = 100
```
**Location**: `test_gravity_inversion.py:99-100`

This creates the baseline geological model with:
- Gravity device locations from the observation data
- Initial structural configuration
- Interpolation parameters (sigmoid_slope controls interface sharpness)

**Normalization**: The forward gravity must be aligned to observed gravity:
```python
norm_params = normalize(
    baseline_fw_gravity_np=baseline(geo_model),
    observed_gravity=observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)
```
**Location**: `test_gravity_inversion.py:101-106`

This handles the **regional-residual separation problem** in gravity data—removing long-wavelength regional trends to focus on local geological features.

---

### Step 3: Define Priors

Priors encode our geological knowledge before seeing the data:

```python
model_priors = {
    "dips": dist.Normal(
        loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),
        scale=torch.tensor(10, dtype=torch.float64),
        validate_args=True
    ),
    "density": dist.Normal(
        loc=torch.tensor([
            2.9,  # plutonites
            2.3   # host rock
        ]),
        scale=torch.tensor(0.15)
    ).to_event(1)
}
```
**Location**: `test_gravity_inversion.py:109-122`

**Interpretation**:
- **Dips prior**: Expect ~10° dip angle, ±10° uncertainty (shallow-dipping layers)
- **Density prior**:
  - Plutonites: 2.9 g/cm³ (typical for granitic intrusions)
  - Host rock: 2.3 g/cm³ (lighter sedimentary/volcanic rocks)
  - Uncertainty: ±0.15 g/cm³ for both

The `.to_event(1)` treats the density vector as a multivariate parameter.

---

### Step 4: Define Deterministic Functions

These are computed quantities we want to track during inference:

```python
post_forward_dets = {
    "gravity_response_raw": lambda samples, gm, sol: sol.gravity,
    "gravity_response": lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),
    "mean_gravity": lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
    "max_gravity": lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
}
```
**Location**: `test_gravity_inversion.py:126-131`

- **gravity_response_raw**: Unnormalized forward gravity
- **gravity_response**: Normalized and aligned to observations (this is what we compare)
- **mean_gravity** / **max_gravity**: Summary statistics

---

### Step 5: Define Likelihood Function

The likelihood connects our model predictions to observations:

```python
likelihood_fn = generate_multigravity_likelihood_hierarchical_per_station(
    norm_params=norm_params
)
```
**Location**: `test_gravity_inversion.py:138-140`

**Why hierarchical per-station?**

Instead of assuming uniform noise across all stations:
```python
# Simple approach (not used):
# y ~ Normal(forward_model, global_sigma)
```

We use a **hierarchical model** that estimates different noise levels for each station:
```python
# Hierarchical approach:
# sigma_i ~ HalfNormal(tau)  # Each station has its own noise
# y_i ~ Normal(forward_model_i, sigma_i)
```

**Benefits**:
1. **Outlier detection**: Stations with high σ indicate problematic measurements or model misfit
2. **Robustness**: Reduces influence of noisy stations
3. **Realism**: Different stations have different measurement quality

See `probabilistic_model_likelihoods.py:generate_multigravity_likelihood_hierarchical_per_station` for implementation.

---

### Step 6: Build the Probabilistic Model

Now we combine everything into a Pyro probabilistic model:

```python
prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
    priors=model_priors,
    set_interp_input_fn=set_priors,
    likelihood_fn=likelihood_fn,
    pre_forward_deterministics={},
    post_forward_deterministics=post_forward_dets,
    obs_name="Gravity Measurement"
)
```
**Location**: `test_gravity_inversion.py:143-151`

**What does this do?**

The `GemPyPyroModel` orchestrates the inference workflow:
1. Sample from priors (dips, densities)
2. Update the GeoModel via `set_priors` function
3. Run forward gravity computation
4. Compute deterministics
5. Evaluate likelihood against observations

The `set_priors` function (from `probabilistic_model.py`) modifies the GemPy model:
```python
def set_priors(samples: dict, geo_model: gp.data.GeoModel):
    """Apply sampled parameters to the geological model"""
    # Update orientations with sampled dips
    # Update rock densities
    # Return modified model
```

---

### Step 7: Prior Predictive Checks

Before running inference, check if the priors are reasonable:

```python
prior_inference_data = gpp.run_predictive(
    prob_model=prob_model,
    geo_model=geo_model,
    y_obs_list=gravity_observations_tensor,
    n_samples=100,
    plot_trace=True
)
```
**Location**: `test_gravity_inversion.py:39-45, 58-64`

**What this tells us**:
1. **Range check**: Do prior predictions cover the observed values?
2. **Model adequacy**: Can *any* parameter combination explain the data?
3. **Prior sensitivity**: How much do predictions vary under the prior?

**Key insight from the test**:
Since we're simulating 20 observations per iteration, some forward models explain certain stations well but fail at others. This suggests:
- The model may be too simple (need more parameters or complexity)
- Some stations may be outliers
- Additional data (e.g., magnetics) might not help much with this model structure

---

### Step 8: Run MCMC Inference

Now we sample from the posterior using NUTS (No-U-Turn Sampler):

```python
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
```
**Location**: `test_gravity_inversion.py:68-86`

**NUTS Configuration Explained**:
- **step_size**: Small step (0.0001) for careful exploration
- **adapt_step_size**: Automatically tune during warmup
- **target_accept_prob**: 65% acceptance (balance exploration vs efficiency)
- **max_tree_depth**: Limit trajectory length (computational budget)
- **num_samples**: 200 posterior samples after warmup
- **warmup_steps**: 200 iterations to tune sampler
- **run_posterior_predictive**: Generate predictions from posterior

**Outputs**:
- `data.posterior`: Parameter samples (dips, densities)
- `data.posterior_predictive`: Forward gravity predictions
- `data.sample_stats`: Diagnostics (divergences, tree depth, etc.)

Save results:
```python
data.to_netcdf("arviz_data_Nov10_I_hierarchical.nc")
```
**Location**: `test_gravity_inversion.py:91`

---

## The Complete Pipeline

Here's the full workflow in `test_gravity_inversion()`:

```python
def test_gravity_inversion(self, simple_geo_model, geophysical_dir,
                           n_samples=50, arviz_data_filename="arviz_data.nc"):
    # 1-6: Create probabilistic model
    geo_model, observed_gravity_ugal, prob_model = self._create_probabilistic_model(
        geophysical_dir, simple_geo_model
    )

    # 7: Prior predictive checks
    prior_inference_data = gpp.run_predictive(...)

    # 8: MCMC inference
    data = gpp.run_nuts_inference(...)

    # Combine prior and posterior
    data.extend(prior_inference_data)

    # Save results
    data.to_netcdf(arviz_data_filename)
```
**Location**: `test_gravity_inversion.py:47-91`

---

## Analysis and Diagnostics

### Diagnostic Checks (`test_run_diagnostics`)

```python
data = az.from_netcdf("arviz_data.nc")
check_mcmc_quality(data, verbose=True, plot=True)
```
**Location**: `test_gravity_inversion.py:154-157`

This checks:
- **Convergence**: R-hat statistics (should be < 1.01)
- **Mixing**: Effective sample size (ESS)
- **Pathologies**: Divergences, maximum tree depth hits
- **Autocorrelation**: How independent are samples?

---

### Parameter Analysis (`test_run_analysis`)

```python
# Plot posterior distributions
az.plot_posterior(data, var_names=["dips", "density"])

# Compare prior vs posterior
az.plot_density(
    data=[data, data.prior],
    var_names=["dips", "density"],
    data_labels=["Posterior", "Prior"],
    colors=[default_red, default_blue],
)
```
**Location**: `test_gravity_inversion.py:237-261`

**What to look for**:
- **Posterior narrowing**: Did we learn from data? (posterior < prior width)
- **Posterior shift**: Did beliefs change from prior?
- **Multimodality**: Multiple plausible solutions?

---

### Predictive Analysis (`test_run_predictive_analysis`)

Compare observed vs predicted gravity:

```python
observed_norm = observed_gravity_ugal
forward_norm = data.prior['gravity_response'].mean(axis=1)
many_forward_norm = data.posterior_predictive['gravity_response'].values[0, -40:-20]

# Plot observed vs forward correlation
for fw in many_forward_norm:
    sorted_fw = fw[sorted_indices]
    ax.scatter(sorted_observed, sorted_fw, alpha=0.7)
    ax.plot(sorted_observed, sorted_fw, alpha=0.3)

# Add 1:1 reference line
ax.plot(lims, lims, 'r--', label='1:1 line')

# Compute correlation
correlation = np.corrcoef(observed_norm, forward_norm)[0, 1]
```
**Location**: `test_gravity_inversion.py:159-208`

**Interpretation**:
- Points near the 1:1 line → good fit
- Scatter away from line → model uncertainty or misfit
- Systematic deviation → model bias (need more complexity)

Use `plot_gravity_comparison()` to visualize spatially:
```python
plot_gravity_comparison(
    observed_ugal=observed_gravity_ugal,
    forward_norm=data.posterior_predictive['gravity_response'].mean(axis=1),
    xy_ravel=xy_ravel,
    normalization_method='align_to_reference'
)
```
**Location**: `test_gravity_inversion.py:278-283`

---

### Outlier Detection (`test_run_outlier_detection`)

The hierarchical likelihood estimates per-station noise:

```python
posterior_sigmas = data.posterior_predictive["sigma_stations"].values

# Identify problematic stations
station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
sigma_global_mean = station_noise_mean.mean()
problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
print(f"Check stations: {problematic}")

# Plot sigma distributions
az.plot_density(
    data=[data, data.prior],
    var_names=["sigma_stations"],
    data_labels=["Posterior", "Prior"],
)
```
**Location**: `test_gravity_inversion.py:210-234`

**Why is this useful?**

Stations with high σ indicate:
1. **Measurement errors**: Instrument malfunction, poor station setup
2. **Geological complexity**: Local features not captured by the model
3. **Outliers**: Data entry errors or extreme conditions

This is **automatic outlier detection**—no manual flagging needed!

---

### Geological Visualization (`gempy_viz`)

Finally, visualize the updated geological model:

```python
gempy_viz(geo_model, data)
```
**Location**: `test_gravity_inversion.py:301`

This shows:
- 3D geological cross-sections with posterior uncertainty
- Layer boundaries with probabilistic envelopes
- Density distributions overlaid on structures

---

## Key Insights

### 1. Model Limitations

**Observation**: Prior predictive shows models that fit some stations but not others.

**Implications**:
- The model is simplified—it may need:
  - More parameters (additional dips, layer thicknesses)
  - Spatial density variations within units
  - More complex structural geometries
- Adding more geophysical data (magnetics) won't help if the forward model is inadequate

---

### 2. Posterior Concentration

**Observation**: Posterior concentrates around models that maximize the number of explained observations.

**Implications**:
- The inversion acts like a **robust regression**
- Outliers (potentially interesting anomalies!) remain unexplained
- To capture outliers, we'd need a more complex model or mixture likelihoods

---

### 3. The Role of the Forward Model

**Key concept**: The forward model (gravity computation) acts as the "polynomial" in a regression:

```
y = f(θ) + ε
```

Where:
- **y**: Observed gravity
- **f(θ)**: Complex forward model (GemPy interpolation + gravity computation)
- **θ**: Parameters (dips, densities)
- **ε**: Noise (station-specific in hierarchical model)

Unlike linear regression where `f` is simple (e.g., `θ₁x + θ₂`), here `f` is a full geological simulation!

---

### 4. Inference Flow Summary

```
Priors (dips, densities)
    ↓
set_priors() → Update GeoModel
    ↓
Forward gravity computation
    ↓
Normalize & align to observations
    ↓
Likelihood evaluation (hierarchical)
    ↓
MCMC sampling (NUTS)
    ↓
Posterior (updated beliefs)
    ↓
Posterior predictive (predictions with uncertainty)
```

---

## Running the Tests

### Full inversion:
```bash
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_gravity_inversion
```

### Prior predictive only:
```bash
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_run_predictive
```

### Analysis (requires saved results):
```bash
# Diagnostics
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_run_diagnostics

# Parameter analysis
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_run_analysis

# Predictive checks
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_run_predictive_analysis

# Outlier detection
pytest test_gravity_inversion.py::TestProbabilisticInversion::test_run_outlier_detection
```

---

## Further Reading

**Key files to explore**:
- `model_setup.py`: Initial GeoModel construction
- `probabilistic_model.py`: Prior setting functions (`set_priors`, `normalize`)
- `probabilistic_model_likelihoods.py`: Likelihood functions (hierarchical model)
- `geophysics.py`: Forward gravity computation and alignment
- `_pyro_runner.py`: MCMC inference engine
- `probabilistic_analysis.py`: Visualization utilities

**Concepts to study**:
- Bayesian inference and MCMC
- Hamiltonian Monte Carlo (HMC) and NUTS
- Hierarchical models for robust inference
- Regional-residual separation in potential fields
- ArviZ for Bayesian diagnostics

---

## Troubleshooting

### Low acceptance probability
- Decrease `step_size` or increase `target_accept_prob`

### Divergences
- Increase `target_accept_prob` (e.g., 0.9)
- Reparameterize priors (use centered vs non-centered)
- Check for numerical issues in forward model

### Slow mixing
- Increase `warmup_steps`
- Run more chains to diagnose
- Consider reparameterization

### Poor fit
- Check prior predictive—is data within prior range?
- Increase model complexity
- Investigate outliers with hierarchical diagnostics

---

**Happy inverting!** 🎉
