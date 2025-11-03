# Gravity Analysis Visualization Guide

This document explains the different visualization plots used for analyzing gravity data, comparing forward models with observations, and understanding uncertainty propagation in probabilistic geophysical inversions.

---

## Table of Contents

1. [Forward Model Comparison Plots](#1-forward-model-comparison-plots)
2. [Uncertainty Visualization Plots](#2-uncertainty-visualization-plots)
3. [Profile Analysis Plots](#3-profile-analysis-plots)
4. [Interpolated Uncertainty Maps](#4-interpolated-uncertainty-maps)
5. [Normalization Methods](#5-normalization-methods)
6. [Interpretation Guidelines](#6-interpretation-guidelines)

---

## 1. Forward Model Comparison Plots

**Function:** `_plot_comparison()`

Creates a 2×2 grid comparing observed gravity with forward model predictions.

### Plot 1: Observed Gravity (Top Left)

**What it shows:**
- Spatial distribution of measured/observed gravity data
- Each point represents an actual measurement location
- Color indicates gravity magnitude

**What to expect:**
- Smooth spatial patterns reflecting geological structures
- Higher values (warm colors) over dense geological bodies
- Lower values (cool colors) over less dense formations
- Possible noise or measurement artifacts

**Interpretation tips:**
- Look for spatial trends (gradients, anomalies)
- Identify regions of high/low gravity
- Note any irregular patterns that might indicate data quality issues

### Plot 2: Forward Model Gravity (Top Right)

**What it shows:**
- Predicted gravity from the geological model
- Same spatial locations as observed data
- Uses same colorbar scale as observed (if normalization is applied)

**What to expect:**
- Smoother patterns than observations (models are inherently smoothed)
- Similar overall spatial trends if model is good
- May miss small-scale features not represented in the geological model

**Interpretation tips:**
- Compare visual patterns with observed data
- Look for areas where model captures main features
- Identify regions where model diverges from observations

### Plot 3: Residuals (Bottom Left)

**What it shows:**
- Difference between observed and predicted gravity: `Residual = Observed - Forward`
- Uses symmetric colorbar (blue = negative, red = positive)

**What to expect:**
- **Near-zero values (white)**: Model fits well
- **Positive values (red)**: Model underpredicts (observed > predicted)
- **Negative values (blue)**: Model overpredicts (predicted > observed)
- Random scatter: Good fit with random errors
- Systematic patterns: Model bias or missing geological features

**Interpretation tips:**
- ✅ **Good model**: Random, small-magnitude residuals
- ⚠️ **Poor model**: Large systematic patterns in residuals
- 🔍 **Investigate**: Clusters of large residuals may indicate:
  - Missing geological features
  - Incorrect density contrasts
  - Model structural errors
  - Data quality issues in specific areas

### Plot 4: Correlation Plot (Bottom Right)

**What it shows:**
- Scatter plot: Observed (x-axis) vs. Forward Model (y-axis)
- Each point represents one measurement location
- Red dashed line: Perfect 1:1 correlation

**What to expect:**
- **Perfect model**: All points on 1:1 line
- **Good model**: Points cluster around 1:1 line
- **Poor model**: Wide scatter, points far from 1:1 line
- **Correlation coefficient (R)**: Displayed in text box (range: -1 to +1)

**Interpretation tips:**
- **R ≈ 1.0**: Excellent correlation (model captures observations well)
- **R = 0.7-0.9**: Good correlation (model captures main trends)
- **R < 0.7**: Poor correlation (model needs improvement)
- **Systematic deviation from 1:1 line**: Model bias (offset or scale issue)
- **Points above line**: Model underpredicts
- **Points below line**: Model overpredicts

---

## 2. Uncertainty Visualization Plots

**Function:** `plot_gravity_with_uncertainty()`

Creates a 2×2 grid showing uncertainty quantification from probabilistic sampling (MCMC, prior predictive).

### Plot 1: Mean Gravity with Error Bars (Top Left)

**What it shows:**
- Mean predicted gravity across all samples
- Error bars scaled by standard deviation (for visibility)
- Color indicates mean gravity magnitude

**What to expect:**
- Mean represents "best estimate" from all samples
- Error bars show relative uncertainty at each location
- Larger error bars = higher uncertainty

**Interpretation tips:**
- Focus on mean values for "best estimate"
- Large error bars indicate model parameter uncertainty propagates strongly to that location
- Small error bars indicate robust predictions regardless of parameter uncertainty

### Plot 2: Prediction Uncertainty (Top Right)

**What it shows:**
- Standard deviation of predictions at each location
- Direct measure of uncertainty magnitude
- Yellow/Red colormap: brighter = more uncertain

**What to expect:**
- **Low uncertainty (dark)**: Model predictions are consistent across parameter samples
- **High uncertainty (bright)**: Model predictions vary widely across parameter samples
- Spatial patterns in uncertainty reveal where model is most/least sensitive

**Interpretation tips:**
- 🎯 **Low uncertainty regions**: Reliable predictions, model structure well-constrained
- ⚠️ **High uncertainty regions**: 
  - Parameter uncertainty strongly affects predictions
  - May need more data in these areas
  - Consider fixing certain parameters if uncertainty is unacceptable
- Compare with data coverage: sparse data often → high uncertainty

### Plot 3: Coefficient of Variation (Bottom Left)

**What it shows:**
- Relative uncertainty: `CV = (std / |mean|) × 100%`
- Accounts for magnitude: 10 μGal uncertainty on 1000 μGal signal is low, but high on 50 μGal signal

**What to expect:**
- **Low CV (< 10%)**: Uncertainty small relative to signal
- **Moderate CV (10-30%)**: Noticeable but manageable uncertainty
- **High CV (> 30%)**: Uncertainty dominates, predictions unreliable

**Interpretation tips:**
- Better metric than absolute uncertainty when signal varies spatially
- High CV in low-gravity regions is common (small signal → large relative uncertainty)
- Use CV to identify where uncertainty matters most for your application

### Plot 4: Observed vs. Predicted with Uncertainty (Bottom Right)

**What it shows:**
- Scatter plot with error bars representing confidence intervals
- Gray error bars: Uncertainty bounds (e.g., 95% CI)
- Red line: Perfect 1:1 correlation

**If observed data provided:**
- **R (correlation coefficient)**: How well mean predictions match observations
- **RMSE (Root Mean Square Error)**: Average prediction error magnitude
- Error bars: Uncertainty ranges around each prediction

**What to expect:**
- Points should cluster around 1:1 line
- Error bars should overlap 1:1 line if model is consistent with data
- Some scatter due to random errors

**Interpretation tips:**
- ✅ **Observations fall within error bars**: Model uncertainty properly captures data variability
- ⚠️ **Observations systematically outside error bars**: 
  - Model uncertainty underestimated
  - Systematic model error not captured by parameter uncertainty
  - Data errors not accounted for
- **RMSE vs. Uncertainty**: If RMSE >> mean uncertainty → systematic model error

---

## 3. Profile Analysis Plots

**Function:** `plot_gravity_uncertainty_profiles()`

Creates 2×2 grid with profiles along different directions, showing gravity variation with uncertainty bands.

### The Four Profiles

1. **X Profile (Top Left)**: West-East transect
2. **Y Profile (Top Right)**: South-North transect
3. **Diagonal 1 (Bottom Left)**: SW-NE transect (X + Y)
4. **Diagonal 2 (Bottom Right)**: NW-SE transect (X - Y)

### What each profile shows:

- **Blue line**: Mean prediction across samples
- **Light blue band**: Confidence interval (e.g., 95%)
- **Gray lines**: Individual samples (subset shown for clarity)
- **Red points**: Observed data (if provided)

**What to expect:**
- Smooth mean prediction line
- Wider confidence bands where uncertainty is higher
- Observations should fall within confidence bands
- Individual samples show variability around mean

**Interpretation tips:**
- **Narrow bands**: Well-constrained predictions
- **Wide bands**: High uncertainty, parameters strongly affect predictions
- **Systematic offset from observations**: Model bias
- **Observations outside bands**: Model inadequacy or underestimated uncertainty
- Look for spatial patterns in uncertainty (e.g., wider at edges)

### Why multiple profiles?

- Gravity anomalies may be elongated (e.g., dikes, faults)
- Different profiles reveal different aspects of 3D structures
- X/Y profiles: Standard cardinal directions
- Diagonal profiles: Capture oblique trends

---

## 4. Interpolated Uncertainty Maps

**Function:** `plot_gravity_uncertainty_map_interpolated()`

Creates smooth, interpolated maps for publication-quality visualization.

### Plot 1: Mean Gravity (Interpolated)

**What it shows:**
- Smooth contour map of mean predicted gravity
- Interpolated between measurement points using cubic interpolation
- Red crosses: Actual measurement locations
- White dashed contours: Observed data (if provided)

**What to expect:**
- Smooth, continuous field (no gaps between points)
- Colors represent gravity magnitude
- Contours show iso-gravity lines

**Interpretation tips:**
- Compare solid (predicted) vs dashed (observed) contours
- Good match → model captures observations
- Offset contours → systematic model error
- Note interpolation artifacts far from measurement points

### Plot 2: Uncertainty Map (Interpolated)

**What it shows:**
- Smooth contour map of prediction standard deviation
- Interpolated uncertainty field
- Blue crosses: Measurement locations

**What to expect:**
- Often lowest uncertainty near measurement clusters
- Higher uncertainty in data gaps
- Smooth gradients in uncertainty

**Interpretation tips:**
- **Low uncertainty (cool colors)**: Reliable predictions
- **High uncertainty (warm colors)**: Need more data or better constraints
- Edge effects: Uncertainty typically increases at boundaries
- Use to plan future surveys: Target high-uncertainty regions

---

## 5. Normalization Methods

### Why normalize?

When comparing observed and forward model gravity:
- May have different units (mGal vs μGal)
- May have systematic offsets (datum issues)
- May have different scales (density errors)

Normalization ensures fair comparison and interpretable colorbars.

### Available Methods

#### `align_to_reference` ⭐ **RECOMMENDED**

**What it does:**