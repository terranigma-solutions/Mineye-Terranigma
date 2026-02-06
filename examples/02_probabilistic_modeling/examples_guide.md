# Probabilistic Inversion Examples Guide

This guide provides an exhaustive description of how Bayesian inversion examples (Sphinx Gallery notebooks) should be structured. The gravity inversion example (`04_gravity_inversion.py`) serves as the gold standard template.

---

## Overview

Each inversion example should be a complete, self-contained tutorial that guides users through the entire Bayesian inversion workflow. The example should be educational, explaining not just *what* to do but *why* each step matters.

**Target audience**: Geoscientists with basic Python knowledge who want to learn probabilistic geological modeling.

**Key principles**:
1. **Comprehensive documentation**: Every code block should have explanatory text
2. **Mathematical foundations**: Include relevant equations with LaTeX formatting
3. **Practical guidance**: Explain parameter choices and their implications
4. **Diagnostic focus**: Show how to verify each step before proceeding
5. **Visualization-rich**: Include plots at every major step

---

## Required Sections

### 1. Title and Introduction (Lines 1-75 in gravity example)

**Purpose**: Set the context and explain what the tutorial accomplishes.

**Required elements**:
- **Title**: Descriptive title with "Bayesian" or "Probabilistic" keyword
- **What are we doing?**: Clear statement of the inversion goal
- **Why Bayesian Inference?**: Explain advantages over deterministic approaches
  - Posterior probability distributions
  - Uncertainty quantification
  - Prior knowledge incorporation
  - Outlier detection (if using hierarchical models)
- **The Bayesian Framework**: Bayes' theorem with mathematical notation
  ```
  p(θ|y) = p(y|θ) p(θ) / p(y)
  ```
- **The Forward Model**: Explain the relationship y = f(θ) + ε

**Example structure**:
```python
"""
Bayesian [Data Type] Inversion: Complete Workflow
==================================================

This tutorial provides a comprehensive walkthrough of Bayesian [data type] inversion...

**What are we doing?**
...

**Why Bayesian Inference?**
...

**The Bayesian Framework**
...

**The Forward Model**
...
"""
```

---

### 2. Import Libraries (Lines 75-115)

**Purpose**: Load all required packages with clear organization.

**Required elements**:
- Environment setup (dotenv, seeds)
- Core scientific packages (numpy, matplotlib)
- GemPy ecosystem (gempy, gempy_probability, gempy_viewer)
- Probabilistic programming (torch, pyro, pyro.distributions)
- Diagnostics (arviz)
- Project-specific helpers

**Best practices**:
- Set random seeds for reproducibility
- Group imports logically
- Add comments for non-obvious imports

```python
# Set random seeds for reproducibility
seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)
np.random.seed(1234)
```

---

### 3. Model Configuration (Lines 115-165)

**Purpose**: Define the geological model extent and resolution.

**Required elements**:
- Model extent (min/max x, y, z coordinates)
- Resolution or refinement level
- Data paths
- Initial GemPy model creation

**Documentation should explain**:
- Why these extent values were chosen
- Trade-offs between resolution and computation time
- How to adapt for different study areas

---

### 4. Data Loading (Lines 165-195)

**Purpose**: Load and describe the observational data.

**Required elements**:
- Data loading function call
- Print statements showing:
  - Number of measurements
  - Data range (min, max)
  - Mean value
  - Units

**Documentation should explain**:
- What the data represents physically
- Units and typical ranges
- Data quality considerations

**Example**:
```python
# %%
# Step 1: Load [Data Type] Observations
# --------------------------------------
#
# **What is [data type] data?**
#
# [Explanation of the physical measurement...]
#
# **Units**: [Units and typical ranges]

data, observed_values = read_data(data_dir)
print(f"\n[Data type] observations:")
print(f"  Number of measurements: {len(observed_values)}")
print(f"  Range: {observed_values.min():.1f} to {observed_values.max():.1f} [units]")
print(f"  Mean: {observed_values.mean():.1f} [units]")
```

---

### 5. Model Setup with Geophysics Configuration (Lines 195-250)

**Purpose**: Configure the GemPy model for forward geophysical computation.

**Required elements**:
- Compute initial model
- Visualize 2D and 3D
- Configure geophysics (measurement locations, parameters)
- Set sigmoid slope for boundary sharpness

**Documentation should explain**:
- What the sigmoid slope parameter controls
- How measurement locations are incorporated
- Any preprocessing steps

---

### 6. Baseline Forward Model (Lines 250-285)

**Purpose**: Compute initial forward model to establish reference.

**Required elements**:
- Baseline computation
- Print forward model statistics
- Explain the forward problem

**Documentation should explain**:
- Purpose of baseline computation
- What the forward model computes
- How to interpret initial results

---

### 7. Normalization/Alignment (Lines 285-325)

**Purpose**: Align forward model to observations.

**Required elements**:
- Normalization function call
- Print normalization parameters
- Visualization comparing observed vs forward

**Documentation should explain**:
- Why normalization is needed (regional-residual problem)
- Available normalization methods
- How to interpret the comparison plots

**Visualization documentation** (plot_geophysics_comparison):
- Panel 1: Observed data spatial distribution
- Panel 2: Forward model predictions
- Panel 3: Residuals (observed - forward)
- Panel 4: Correlation plot with R coefficient

---

### 8. Prior Distributions (Lines 325-385)

**Purpose**: Define prior beliefs about model parameters.

**Required elements**:
- Prior dictionary with all parameters
- Mathematical notation for each prior
- Print statements summarizing priors

**Documentation should explain**:
- What each parameter represents physically
- Why specific distributions were chosen
- How to interpret prior parameters (mean, std)
- Difference between informative and uninformative priors

**Example**:
```python
# %%
# Step 5: Define Prior Distributions
# -----------------------------------
#
# **Prior Knowledge in Bayesian Inference**
#
# Priors encode our geological knowledge *before* seeing the data...
#
# **Prior on [Parameter]**
#
# We use a [Distribution] distribution:
#
# .. math::
#
#     \text{parameter} \sim \mathcal{N}(\mu, \sigma)
#
# **Interpretation**:
# - **Mean (μ)**: [Why this value]
# - **Std (σ)**: [What uncertainty this represents]

model_priors = {
    'parameter_name': dist.Normal(
        loc=torch.tensor(...),
        scale=torch.tensor(...)
    ).to_event(1),
}
```

---

### 9. Likelihood Function (Lines 405-540) ⭐ CRITICAL SECTION

**Purpose**: Define how model predictions connect to observations.

**This is the most variable component across different inversions and requires the most detailed documentation.**

**Required elements**:
- Likelihood function creation
- Mathematical formulation
- Explanation of noise model

**Documentation MUST explain**:

#### 9.1 Basic Likelihood Concept
- What the likelihood represents: p(data|parameters)
- The noise model assumption
- Mathematical notation

#### 9.2 Available Likelihood Functions

Document each likelihood function available in `probabilistic_model_likelihoods.py`:

**Diagonal Likelihood** (`generate_multigravity_likelihood_diagonal`):
```
y_i ~ Normal(f_i(θ), σ²)
```
- Assumes independent noise across stations
- Fixed noise level
- Simplest option, good for initial testing

**Hierarchical Per-Station** (`generate_multigravity_likelihood_hierarchical_per_station`):
```
Level 1 (Hyperpriors):
  μ_log_σ ~ Normal(log(5000), 0.5)
  τ_log_σ ~ HalfNormal(0.5)

Level 2 (Station-specific noise):
  log(σ_i) ~ Normal(μ_log_σ, τ_log_σ)

Level 3 (Observations):
  y_i ~ Normal(f_i(θ), σ_i²)
```
- Each station has its own noise level
- Automatic outlier detection
- Partial pooling between stations

**Hierarchical with Spatial Correlation** (`generate_multigravity_likelihood_hierarchical`):
- Includes spatial covariance structure
- Uses Gaussian kernel for correlation
- Student-t distribution for robustness

#### 9.3 Showing Likelihood Function Source Code

**IMPORTANT**: To help users understand what the likelihood function does internally, include a cell that displays the function's source code:

```python
# %%
# **Inspecting the Likelihood Function**
#
# To understand exactly what the likelihood function does, we can inspect its source code:

import inspect
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import (
    generate_multigravity_likelihood_hierarchical_per_station
)

print("Likelihood function source code:")
print("=" * 60)
print(inspect.getsource(generate_multigravity_likelihood_hierarchical_per_station))
```

This allows users to see:
- How simulated data is aligned to observations
- What distributions are used
- What parameters are sampled vs fixed
- How the final distribution is constructed

#### 9.4 Parameter Interpretation

For hierarchical models, explain each parameter:
- **μ_log_σ**: Global mean noise level (log scale)
- **τ_log_σ**: Variability between stations
- **σ_i**: Per-station noise standard deviation

#### 9.5 Benefits and Trade-offs

Compare different likelihood choices:
| Likelihood Type | Pros | Cons | Use When |
|----------------|------|------|----------|
| Diagonal | Simple, fast | No outlier handling | Clean data |
| Hierarchical per-station | Robust, automatic outliers | More parameters | Real-world data |
| Spatial correlation | Captures correlations | Computationally expensive | Dense station networks |

---

### 10. Probabilistic Model Creation (Lines 540-715)

**Purpose**: Combine all components into a Pyro model.

**Required elements**:
- `make_gempy_pyro_model` or `make_gempy_pyro_model_extended` call
- Print statement confirming model creation

**Documentation should explain**:
- What each argument does:
  - `priors`: Prior distributions dictionary
  - `set_interp_input_fn`: Bridge between Pyro and GemPy
  - `likelihood_fn`: Connects predictions to observations
  - `obs_name`: Label for ArviZ output
- The internal workflow (sample → update model → forward → likelihood)
- Why this architecture enables automatic differentiation

**Include workflow diagram**:
```
Sample θ from priors
       ↓
set_interp_input_fn(θ) → Update GeoModel
       ↓
Forward model computation
       ↓
likelihood_fn(Solutions) → p(data|θ)
       ↓
MCMC uses gradients to propose next θ
```

---

### 11. Prior Predictive Checks (Lines 715-845)

**Purpose**: Verify model adequacy before expensive inference.

**Required elements**:
- `gpp.run_predictive()` call
- Geological visualization (`gempy_viz`)
- Forward model comparison (`plot_many_observed_vs_forward`)

**Documentation should explain**:
- Why prior predictive checks matter
- What to look for:
  - Do prior predictions cover observed values?
  - Can ANY parameter combination explain the data?
  - How variable are predictions under the prior?
- How to interpret each visualization

**Key questions to answer**:
1. Range check: Do predictions span observations?
2. Model adequacy: Can the model explain the data?
3. Prior sensitivity: How much do predictions vary?

---

### 12. MCMC Inference (Lines 845-960)

**Purpose**: Run the actual Bayesian inference.

**Required elements**:
- `gpp.run_nuts_inference()` call with `NUTSConfig`
- Option to load pre-computed results
- Print statements showing progress

**Documentation should explain**:
- How NUTS works (Hamiltonian Monte Carlo variant)
- Each configuration parameter:
  - `step_size`: Leapfrog integration step
  - `adapt_step_size`: Automatic tuning during warmup
  - `target_accept_prob`: Target acceptance rate
  - `max_tree_depth`: Computational budget per iteration
  - `init_strategy`: Starting point for chains
  - `num_samples`: Posterior samples after warmup
  - `warmup_steps`: Adaptation iterations
  - `num_chains`: Independent chains for diagnostics

**Include pre-computed results option**:
```python
RUN_SIMULATION = False
if RUN_SIMULATION:
    data = gpp.run_nuts_inference(...)
    data.to_netcdf("results.nc")
else:
    data = az.from_netcdf("results.nc")
```

---

### 13. Posterior Analysis (Lines 960-1140)

**Purpose**: Analyze and interpret inference results.

**Required elements**:
- Parameter statistics (mean, std, median)
- Comparison with prior (uncertainty reduction)
- Predictive performance (residuals, RMS)
- Visualizations:
  - `az.plot_density()` for prior vs posterior
  - `plot_many_observed_vs_forward()` for posterior predictions
  - `gempy_viz()` for geological uncertainty

**Documentation should explain**:
- How to interpret posterior statistics
- What uncertainty reduction means
- How to assess model fit quality
- What residual patterns indicate

**Key metrics**:
- Uncertainty reduction: `(1 - posterior_std / prior_std) * 100%`
- Mean residual: Should be ≈ 0 (unbiased)
- RMS residual: Overall misfit measure

---

### 14. Advanced Diagnostics (Lines 1140-1220)

**Purpose**: Deeper analysis for hierarchical models and uncertainty.

**Required elements** (if applicable):
- Sigma analysis for outlier detection
- Probability density fields
- Information entropy visualization
- 3D uncertainty visualization

**Documentation should explain**:
- How to identify problematic stations
- What entropy represents
- How to interpret probability fields

---

### 15. Summary and Conclusions (Lines 1220-1368)

**Purpose**: Recap the workflow and provide guidance for extensions.

**Required elements**:
- Workflow diagram (ASCII art)
- Key functions and their roles
- Key insights from the tutorial
- Next steps for production inversions
- Available diagnostic checks

**Workflow diagram example**:
```
┌─────────────────────────────────────────────────────────────────┐
│ 1. PRIOR DEFINITION                                             │
│    Define prior distributions for parameters                    │
└─────────────────────────────────────────────────────────────────┘
                               ↓
┌─────────────────────────────────────────────────────────────────┐
│ 2. PROBABILISTIC MODEL CONSTRUCTION                             │
│    make_gempy_pyro_model() orchestrates all components          │
└─────────────────────────────────────────────────────────────────┘
                               ↓
... [continue for all steps]
```

**Key insights to include**:
1. Hierarchical modeling benefits
2. Importance of prior predictive checks
3. Residual pattern interpretation
4. Connection to regression problems
5. Value of multiple visualization types

**Next steps suggestions**:
- Multiple chains for convergence diagnostics
- More samples for publication quality
- Additional parameters to invert
- Joint inversion with multiple data types
- Model comparison techniques

---

## Likelihood Functions Reference

The likelihood functions in `probabilistic_model_likelihoods.py` are the most variable component across inversions. Here's a complete reference:

### Gravity Likelihoods

| Function | Description | Use Case |
|----------|-------------|----------|
| `MultiGravityDiagonalLikelihood` | Dataclass with fixed diagonal noise | Simple, fast inference |
| `generate_multigravity_likelihood_diagonal` | Function returning diagonal likelihood | Testing, clean data |
| `generate_multigravity_likelihood_hierarchical_per_station` | Per-station noise inference | **Recommended for real data** |
| `generate_multigravity_likelihood_per_station_stable` | Bounded per-station noise | VI stability |
| `generate_multigravity_likelihood_hierarchical` | Spatial correlation + Student-t | Dense networks |

### Magnetic Likelihoods

| Function | Description | Use Case |
|----------|-------------|----------|
| `generate_multimagnetic_likelihood_hierarchical_per_station` | Per-station noise for magnetics | Real magnetic data |

### Remote Sensing Likelihoods

| Function | Description | Use Case |
|----------|-------------|----------|
| `enmap_likelihood_fn` | Ordinal classification likelihood | EnMAP spectral data |

### How to Display Likelihood Source in Examples

Always include a cell showing the likelihood function source code:

```python
# %%
# **Understanding the Likelihood Function**
#
# The likelihood function is crucial for connecting model predictions to observations.
# Let's inspect what it does internally:

import inspect
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import (
    generate_multigravity_likelihood_hierarchical_per_station  # or your chosen likelihood
)

# Display the function signature and docstring
help(generate_multigravity_likelihood_hierarchical_per_station)

# %%
# **Full Source Code**
#
# For complete transparency, here's the full implementation:

print(inspect.getsource(generate_multigravity_likelihood_hierarchical_per_station))
```

This helps users:
1. Understand exactly what distributions are used
2. See what parameters are sampled vs fixed
3. Learn how to create custom likelihoods
4. Debug issues with their own inversions

---

## Code Style Guidelines

### Sphinx Gallery Formatting

- Use `# %%` for cell separators
- Use `# %%` followed by title for section headers
- Use `#` comments for documentation (converted to RST)
- Use `.. math::` for LaTeX equations
- Use `.. code-block:: python` for code examples in docs

### Documentation Density

- Every code cell should have preceding documentation
- Explain *why*, not just *what*
- Include interpretation guidance for outputs
- Add print statements showing key values

### Visualization Standards

- Always label axes with units
- Use consistent colormaps (viridis_r for data, RdBu_r for residuals)
- Include colorbars
- Add titles describing what the plot shows

---

## Checklist for New Inversion Examples

Before submitting a new inversion example, verify:

- [ ] Title clearly states "Bayesian [Data Type] Inversion"
- [ ] Introduction explains Bayesian framework and forward model
- [ ] All imports are organized and seeds are set
- [ ] Data loading includes statistics printout
- [ ] Priors have mathematical notation and interpretation
- [ ] **Likelihood function is documented with source code display**
- [ ] Prior predictive checks are performed and explained
- [ ] NUTS configuration parameters are documented
- [ ] Pre-computed results option is available
- [ ] Posterior analysis includes uncertainty reduction metrics
- [ ] Summary includes workflow diagram
- [ ] All visualizations have explanatory text
- [ ] Code follows Sphinx Gallery formatting

---

## Example Length Guidelines

| Example Type | Target Lines | Sections |
|--------------|--------------|----------|
| Full tutorial (like gravity) | 1200-1500 | All 15 sections |
| Standard inversion | 600-900 | Sections 1-13, 15 |
| Quick example | 250-400 | Abbreviated sections |

The gravity inversion example (1368 lines) represents the gold standard for comprehensive documentation. Other examples can be shorter but should maintain the same structure and documentation quality for the sections they include.
