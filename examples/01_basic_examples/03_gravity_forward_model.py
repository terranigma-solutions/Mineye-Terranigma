"""
Forward Gravity Modeling - Tharsis Region
==========================================

This example demonstrates forward gravity modeling using a 3D geological model.
We compute the gravity response of the geological structure at real measurement locations.

**Overview**

Forward gravity modeling predicts the gravitational field anomaly caused by subsurface
density contrasts. This is a crucial tool in geophysical exploration for:

* Validating geological interpretations against field measurements
* Identifying density anomalies associated with ore deposits or intrusions
* Constraining geological models before inversion
* Understanding the relationship between geology and geophysical signatures

**Geophysical Background**

Gravity surveys measure tiny variations in Earth's gravitational field caused by lateral
density variations in the subsurface. The gravity anomaly (Δg) at a point depends on:

* The **density contrast** between geological units (Δρ)
* The **volume** and **geometry** of anomalous bodies
* The **distance** from the observation point to the density anomaly

Typical units:

* mGal (milligal) = 10⁻⁵ m/s² = 10⁻³ cm/s²
* μGal (microgal) = 10⁻⁸ m/s² = 10⁻⁶ cm/s²

Modern gravimeters can measure to ~10 μGal precision.

**The Tharsis Case Study**

The Tharsis plutonite intrusion (density ~2.9 g/cm³) intruded into sedimentary host rocks
(density ~2.3 g/cm³), creating a positive gravity anomaly. This density contrast of
~0.6 g/cm³ produces a measurable signal that we can model and compare with observations.

**Workflow Summary**

1. Create geological model (similar to Example 01)
2. Load field gravity measurements
3. Set up computation grids at measurement locations
4. Assign density values to geological units
5. Compute forward gravity response
6. Compare with observations and analyze residuals
"""

# %%
# Import Libraries
# ----------------
#
# **Core Libraries**:
#
# * **NumPy**: Numerical array operations
# * **PyTorch**: Tensor operations with automatic differentiation (useful for inversions)
# * **Matplotlib**: Plotting and visualization
#
# **Geological Modeling**:
#
# * **GemPy**: 3D implicit geological modeling framework
# * **GeoPandas**: Geospatial data handling (for loading gravity measurements with coordinates)
#
# **Mineye Utilities**:
#
# * ``align_forward_to_observed``: Aligns modeled gravity to observed data statistics
# * ``normalize``: Data normalization and preprocessing functions

import dotenv
dotenv.load_dotenv()

from mineye.GeoModel.helper_plotter import plot_model_and_gravity_sensors
import numpy as np
import torch
from matplotlib import pyplot as plt

import gempy as gp
import geopandas as gpd

from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize

# Set random seed for reproducibility
np.random.seed(1234)

# %%
# Import paths configuration
from mineye.config import paths

# %%
# Define Model Extent and Resolution
# -----------------------------------
# The model covers the Tharsis mining district

min_x = -707_521  # * Cropping the corrupted area of the geotiff 
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
extent = [min_x, max_x, min_y, max_y, model_depth, max_z]

# Model resolution: use octree with refinement level 5
resolution = None  # Using octrees instead of regular grid
refinement = 5

# %%
# Get Data Paths
# --------------
# Load paths to structural and geophysical data

mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
gravity_data_path = paths.get_gravity_data_path()

print(f"Orientations: {mod_or_path}")
print(f"Points: {mod_pts_path}")
print(f"Gravity data: {gravity_data_path}")

# %%
# Create GemPy Geological Model
# ------------------------------
# Build the model structure with imported structural data

simple_geo_model = gp.create_geomodel(
    project_name='gravity_model',
    extent=extent,
    refinement=refinement,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

# %%
# Map Geological Units
# --------------------
# Define the stratigraphic stack

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={
            "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }
)

# %%
# Load Gravity Observations
# --------------------------
# Read actual gravity measurements from the field
#
# **Gravity Survey Data**:
#
# Gravity surveys involve systematic measurements of gravitational acceleration across
# an area of interest. The data typically includes:
#
# * **Spatial coordinates** (X, Y, Z): Location of each measurement station
# * **Bouguer anomaly**: Gravity after correcting for latitude, elevation, terrain, and tidal effects
# * **Free-air anomaly**: Simpler correction (elevation only)
# * **Measurement metadata**: Date, instrument ID, drift corrections, etc.
#
# The Bouguer anomaly (``VALU_BOU267``) represents the gravity signal attributed purely to
# subsurface density variations, making it ideal for geological interpretation.
#
# **Data Quality Considerations**:
#
# * Survey spacing and coverage affect resolution
# * Measurement precision (typically 10-50 μGal for ground surveys)
# * Systematic errors from incomplete corrections
# * Regional trends vs. local anomalies

gravity_data = gpd.read_file(gravity_data_path)
observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

print(f"Number of gravity observations: {len(observed_gravity)}")
print(f"Gravity range: {observed_gravity.min():.2f} to {observed_gravity.max():.2f} mGal")

# %%
# Set Up Gravity Measurement Locations
# -------------------------------------
# Create measurement grid at actual device locations
#
# **Centered Grid Approach**:
#
# For efficient gravity computation, GemPy uses "centered grids" - small 3D volumes
# around each measurement point. This approach:
#
# * Focuses computation only where observations exist (not the entire model volume)
# * Maintains high resolution near observation points
# * Significantly reduces computational cost for large numbers of stations
#
# Each measurement location becomes the center of a local grid where the density
# distribution is evaluated in detail.

xy_ravel = np.column_stack([
        np.array(gravity_data.geometry.x.values),
        np.array(gravity_data.geometry.y.values),
        np.full(len(gravity_data), 505)  # Set Z to surface elevation
])

print(f"Using {len(xy_ravel)} actual measurement points")

# %%
# Configure Centered Grid for Gravity Computation
# ------------------------------------------------
# Set up centered grid around each measurement point
#
# **Grid Parameters**:
#
# * **centers**: Coordinates of measurement stations
# * **resolution**: Number of cells in each direction [X, Y, Z]
#   - [10, 10, 15] creates 10×10×15 = 1,500 cells per station
# * **radius**: Extent of local grid in meters [X, Y, Z]
#   - [5000, 5000, 5000] = 5 km radius in all directions
#
# The grid resolution balances accuracy vs. computation time. Higher resolution captures
# finer details of the density distribution but increases computation cost.
#
# The 5 km radius is chosen to capture the significant gravity contributions from the
# plutonite body while keeping grids manageable. Distant mass has diminishing influence
# due to the 1/r² decay of gravitational attraction.

gp.set_centered_grid(
    grid=simple_geo_model.grid,
    centers=xy_ravel,
    resolution=np.array([10, 10, 15]),
    radius=np.array([2000, 5000, 2000])
)

# %%
# Calculate Gravity Gradient
# ---------------------------
# Compute the gravity gradient (tz component) for the centered grid
#
# **Gravity Tensor Components**:
#
# The gravity field can be represented as a gradient tensor with 9 components.
# For vertical gravity (what gravimeters measure), we primarily use:
#
# * **g_z** (or **t_z**): Vertical component of gravitational acceleration
# * This is what ground-based gravimeters measure
#
# The gradient calculation determines the geometric contribution of each grid cell
# to the total gravity at the observation point, based on:
#
# * Distance from cell to observation point
# * Cell volume
# * Orientation (vertical direction)
#
# These gradients are later multiplied by density values to get the final gravity response.

gravity_gradient = gp.calculate_gravity_gradient(simple_geo_model.grid.centered_grid)
print(f"Gravity gradient tensor shape: {gravity_gradient.shape}")

# %%
# Set Density Values
# ------------------
# Define densities for geological units
#
# **Rock Densities**:
#
# Density contrasts drive gravity anomalies. Typical rock densities:
#
# * **Plutonic rocks** (granite, diorite): 2.60-2.95 g/cm³
# * **Sedimentary rocks** (shale, sandstone): 2.20-2.70 g/cm³
# * **Volcanic rocks** (basalt): 2.70-3.10 g/cm³
# * **Ore deposits** (massive sulfides): 3.50-5.00 g/cm³
#
# **Tharsis Densities**:
#
# * **Tournaisian Plutonites**: 2.9 g/cm³ (dense intrusive rocks)
# * **Devonian Sedimentary Host**: 2.3 g/cm³ (lower density host rocks)
# * **Density contrast**: +0.6 g/cm³ (creates positive gravity anomaly)
#
# These values are based on:
#
# * Laboratory measurements of rock samples
# * Literature values for similar rock types
# * Calibration with gravity data (in inverse modeling)
#
# .. note::
#    Density uncertainties significantly impact modeling results. Typical uncertainties
#    are ±0.05-0.10 g/cm³, which should be considered in probabilistic inversions.

density_plutonites = 2.9 - 2.67  # g/cm³
density_sedimentary_host = 2.3 - 2.67 # g/cm³

simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
    tz=gravity_gradient,
    densities=np.array([density_plutonites, density_sedimentary_host])
)

print(f"Plutonite density: {density_plutonites} g/cm³")
print(f"Sedimentary host density: {density_sedimentary_host} g/cm³")
print(f"Density contrast: {density_plutonites - density_sedimentary_host:.1f} g/cm³")

# %%
# Compute Forward Gravity Model
# ------------------------------
# Run the interpolation and gravity computation
#
# **Forward Modeling Process**:
#
# 1. **Geological interpolation**: Determine which geological unit occupies each grid cell
# 2. **Density assignment**: Assign the appropriate density to each cell based on lithology
# 3. **Gravity summation**: Sum the gravitational contribution of all cells using:
#
#    .. math::
#       g_z = G \\sum_{i} \\frac{\\Delta\\rho_i \\cdot V_i \\cdot z_i}{r_i^3}
#
#    where:
#
#    * G = gravitational constant (6.674×10⁻¹¹ m³/(kg·s²))
#    * Δρᵢ = density of cell i
#    * Vᵢ = volume of cell i
#    * zᵢ = vertical distance from cell to observation point
#    * rᵢ = total distance from cell to observation point
#
# 4. **Unit conversion**: Convert to mGal or μGal
#
# The ``mesh_extraction=False`` flag skips 3D mesh generation (not needed for gravity).

simple_geo_model.interpolation_options.mesh_extraction = True

# topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')
gp.set_topography_from_file(
    grid=simple_geo_model.grid,
    filepath=paths.get_topography_path(),
    crop_to_extent=[
            simple_geo_model.grid.extent[0],
            simple_geo_model.grid.extent[2],
            simple_geo_model.grid.extent[1],
            simple_geo_model.grid.extent[3]
    ]
)

sol = gp.compute_model(simple_geo_model)

print("✓ Forward gravity model computed successfully")

# %%
# Extract Gravity Results
# -----------------------
# Get the computed gravity response
#
# **Data Alignment**:
#
# The ``align_forward_to_observed`` function adjusts the forward model to match the
# statistical properties of observations. This accounts for:
#
# * **Regional trends**: Long-wavelength variations not captured by the local model
# * **Reference level differences**: Arbitrary datums in gravity data
# * **Systematic biases**: Uncorrected effects or model simplifications
#
# The normalization uses the ``align_to_reference`` method with a 30% extrapolation
# buffer to handle edge effects.
#
# .. note::
#    This alignment is a practical necessity for forward modeling but should be
#    carefully considered in inverse problems where we're trying to fit the data.

observed_gravity_ugal = observed_gravity * 1000  # Convert mGal to μGal
norm_params = normalize(
    sol.gravity,
    observed_gravity_ugal,
    method="align_to_reference",
    extrapolation_buffer=0.3
)
grav = align_forward_to_observed(sol.gravity, norm_params)

print(f"\nComputed gravity values:")
print(f"  Shape: {grav.shape}")
print(f"  Range: {grav.min():.2f} to {grav.max():.2f} μGal")
print(f"  Mean: {grav.mean():.2f} μGal")

# %%
# Compare with Observations
# -------------------------
# Calculate residuals between observed and computed gravity
#
# **Residual Analysis**:
#
# Residuals (observed - modeled) reveal how well the geological model explains the data:
#
# * **Mean residual**: Should be close to zero after alignment; systematic bias indicates
#   regional trends or model errors
# * **Standard deviation**: Measures scatter; affected by model accuracy and data noise
# * **RMSE** (Root Mean Square Error): Overall measure of fit quality
#
#   .. math::
#      RMSE = \\sqrt{\\frac{1}{N} \\sum_{i=1}^{N} (g_{obs,i} - g_{mod,i})^2}
#
# **Interpreting Residuals**:
#
# * **Systematic patterns**: Suggest missing geological features or incorrect geometry
# * **Random scatter**: Dominated by measurement noise or small-scale geology
# * **Large residuals**: May indicate:
#   - Incorrect density values
#   - Unmodeled geological structures
#   - Data quality issues
#   - 3D effects not captured by 2D modeling

if isinstance(grav, torch.Tensor):
    grav = grav.detach().numpy()

residuals = observed_gravity_ugal - grav

print(f"\nGravity residuals:")
print(f"  Mean: {residuals.mean():.2f} μGal")
print(f"  Std: {residuals.std():.2f} μGal")
print(f"  RMS: {np.sqrt(np.mean(residuals ** 2)):.2f} μGal")

# %%
# Visualize the Model
# -------------------
# Create 2D and 3D visualizations of the geological model
#
# **Model Visualization**:
#
# Before analyzing geophysical results, it's important to visualize the geological
# model itself to understand:
#
# * Geometry of geological units
# * Position of observation points relative to the modeled structures
# * Topographic influences
# * Overall model plausibility

import gempy_viewer as gpv
import pyvista as pv

# Configure PyVista to use smaller fonts for axes
pv.global_theme.font.size = 10
pv.global_theme.font.label_size = 10
pv.global_theme.font.title_size = 12

gpv.plot_2d(simple_geo_model)




plot_model_and_gravity_sensors(simple_geo_model, xy_ravel, grav)

# %%
# Visualize Gravity Results
# --------------------------
# Plot forward gravity and comparison with observations
#
# **Geophysical Data Visualization**:
#
# Spatial plots help identify:
#
# * Anomaly patterns and their relationship to geology
# * Data coverage and survey design
# * Edge effects and boundary artifacts
# * Correlation between observed and modeled signatures

from mineye.GeoModel.plotting.probabilistic_analysis import plot_fw_geophysics, plot_comparison

plot_fw_geophysics(
    fw_values=grav,
    observed_data=gravity_data,
    xy_ravel=xy_ravel,
    label=r'Forward Model Gravity ($\mu$gal)',
    title='Forward Model Gravity Results'
)
plot_comparison(observed_gravity, grav, xy_ravel, unit_label=r'$\mu$Gal')

# %%
# Visualize Forward Model Results
# --------------------------------
#
# **Detailed Gravity Map**:
#
# This map shows the predicted gravity field from our geological model.
# Key features to observe:
#
# * **Positive anomalies**: Occur over the dense plutonite intrusion
# * **Negative anomalies**: May occur at the edges or over less dense rocks
# * **Spatial trends**: Reflect the geometry and extent of density contrasts
# * **Amplitude**: Related to the magnitude of density contrast and body volume
#
# The colormap (viridis_r) is reversed so higher values (positive anomalies) appear
# warmer, following common geophysical convention.

grav_forward = grav
print(f"✓ Forward model computed")
print(f"  Gravity range: [{grav_forward.min():.2f}, {grav_forward.max():.2f}] µGal")

fig, ax = plt.subplots(figsize=(10, 8))

scatter = ax.scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=grav_forward, s=50,
    cmap='viridis_r', alpha=0.8,
    edgecolors='black', linewidth=0.5
)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label(r'Forward Model Gravity (µGal)', fontsize=12)

ax.set_title('Forward Gravity Model', fontsize=14, fontweight='bold')
ax.set_xlabel('X (m)')
ax.set_ylabel('Y (m)')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.show()

# %%
# Compare with Observed Gravity
# ------------------------------
#
# **Three-Panel Comparison**:
#
# This comprehensive view shows:
#
# 1. **Observed Gravity** (Left):
#    - Actual field measurements
#    - Represents the "truth" we're trying to match
#    - Contains geological signal + noise + regional trends
#
# 2. **Forward Model** (Center):
#    - Predicted gravity from our geological interpretation
#    - Shows what the model "thinks" the gravity should be
#    - Ideally should match observed patterns
#
# 3. **Residuals** (Right):
#    - Difference between observed and modeled (Obs - Model)
#    - **Red**: Model under-predicts (observed > modeled)
#    - **Blue**: Model over-predicts (observed < modeled)
#    - **White/near-zero**: Good fit
#
# **What to look for**:
#
# * Systematic patterns in residuals suggest model improvements needed
# * Random residuals indicate we're at the noise level
# * Large residuals in specific areas point to local geological features not captured

# Convert observed from mGal to µGal
observed_ugal = observed_gravity_ugal

# Create comparison plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

# Observed
sc1 = axes[0].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=observed_ugal, s=50, cmap='viridis_r',
    edgecolors='black', linewidth=0.5
)
axes[0].set_title('Observed Gravity')
axes[0].set_xlabel('X (m)')
axes[0].set_ylabel('Y (m)')
plt.colorbar(sc1, ax=axes[0], label='µGal')

# Forward model
sc2 = axes[1].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=grav_forward, s=50, cmap='viridis_r',
    edgecolors='black', linewidth=0.5
)
axes[1].set_title('Forward Model')
axes[1].set_xlabel('X (m)')
plt.colorbar(sc2, ax=axes[1], label='µGal')

# Residuals
sc3 = axes[2].scatter(
    xy_ravel[:, 0], xy_ravel[:, 1],
    c=residuals, s=50, cmap='RdBu_r',
    edgecolors='black', linewidth=0.5
)
axes[2].set_title('Residuals (Obs - Model)')
axes[2].set_xlabel('X (m)')
plt.colorbar(sc3, ax=axes[2], label='µGal')

plt.tight_layout()
plt.show()

# %%
# Correlation Analysis
# --------------------
#
# **Cross-Plot Analysis**:
#
# The scatter plot of observed vs. modeled gravity is a powerful diagnostic tool:
#
# * **1:1 line** (red dashed): Perfect agreement would have all points on this line
# * **Correlation coefficient (R)**: Measures linear relationship strength (-1 to +1)
#   - R > 0.9: Excellent correlation
#   - R = 0.7-0.9: Good correlation
#   - R < 0.7: Poor correlation
# * **R²**: Fraction of variance explained by the model (0 to 1)
# * **RMSE**: Average prediction error in µGal
#
# **Interpreting the Plot**:
#
# * Points above the 1:1 line: Model under-predicts
# * Points below the 1:1 line: Model over-predicts
# * Scatter around the line: Combination of model error and measurement noise
# * Systematic deviations: Indicate model bias or missing physics
#
# This analysis helps determine if the geological model is a reasonable representation
# of the subsurface structure.

fig, ax = plt.subplots(figsize=(8, 8))

ax.scatter(observed_ugal, grav_forward, alpha=0.6, s=60,
           edgecolors='black', linewidth=0.5)

# 1:1 line
lims = [
        min(ax.get_xlim()[0], ax.get_ylim()[0]),
        max(ax.get_xlim()[1], ax.get_ylim()[1])
]
ax.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')

# Calculate statistics
correlation = np.corrcoef(observed_ugal, grav_forward)[0, 1]
rmse = np.sqrt(np.mean(residuals ** 2))

ax.set_xlabel('Observed (µGal)', fontsize=12)
ax.set_ylabel('Forward Model (µGal)', fontsize=12)
ax.set_title('Observed vs Modeled Gravity', fontsize=14, fontweight='bold')
ax.grid(True, alpha=0.3)
ax.legend()

# Add statistics text box
textstr = f'R = {correlation:.3f}\nR² = {correlation ** 2:.3f}\nRMSE = {rmse:.2f} µGal'
ax.text(0.05, 0.95, textstr, transform=ax.transAxes,
        verticalalignment='top',
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
        fontsize=11)

plt.tight_layout()
plt.show()

# %%
# Print Summary Statistics
# -------------------------

print("\n" + "=" * 50)
print("GRAVITY COMPARISON STATISTICS")
print("=" * 50)
print(f"Number of measurements: {len(observed_ugal)}")
print(f"\nObserved gravity:")
print(f"  Mean: {observed_ugal.mean():.2f} µGal")
print(f"  Std:  {observed_ugal.std():.2f} µGal")
print(f"  Range: [{observed_ugal.min():.2f}, {observed_ugal.max():.2f}] µGal")
print(f"\nForward model:")
print(f"  Mean: {grav_forward.mean():.2f} µGal")
print(f"  Std:  {grav_forward.std():.2f} µGal")
print(f"  Range: [{grav_forward.min():.2f}, {grav_forward.max():.2f}] µGal")
print(f"\nResiduals:")
print(f"  Mean: {residuals.mean():.2f} µGal")
print(f"  Std:  {residuals.std():.2f} µGal")
print(f"  RMSE: {rmse:.2f} µGal")
print(f"  MAE:  {np.abs(residuals).mean():.2f} µGal")
print(f"\nCorrelation: {correlation:.4f} (R² = {correlation ** 2:.4f})")
print("=" * 50)

# %%
# Summary and Interpretation
# ---------------------------

# **Model Limitations**:
#
# Forward models are simplifications of reality and have inherent limitations:
#
# * **Density assumptions**: We used uniform densities, but real rocks are heterogeneous
# * **Geological simplification**: Simple geometry may not capture all structural complexity
# * **Data coverage**: Limited measurement points constrain resolution
# * **Regional trends**: Long-wavelength features may not be fully captured
# * **3D effects**: Some anomalies may be caused by out-of-plane structures
#
# **Next Steps and Applications**:
#
# This forward modeling foundation enables:
#
# * **Probabilistic Inversion** (Example 04): Use Bayesian methods to infer densities
#   and geological parameters from the gravity data
# * **Uncertainty Quantification**: Propagate uncertainties in geology to geophysics
# * **Joint Inversion**: Combine gravity with magnetics, seismic, or other data
# * **Model Refinement**: Use residuals to iteratively improve geological interpretation
# * **Resource Estimation**: Better constrain ore body geometry and properties
#
# **Further Reading**:
#
# * Blakely, R. J. (1996). Potential Theory in Gravity and Magnetic Applications.
#   Cambridge University Press.
# * Li, Y., & Oldenburg, D. W. (1998). 3-D inversion of gravity data. Geophysics, 63(1), 109-119.
# * de la Varga et al. (2019). GemPy 1.0: open-source stochastic geological modeling
#   and inversion. Geoscientific Model Development, 12(1), 1-32.
#
# .. seealso::
#    * :doc:`../examples_probabilistic/04_gravity_inversion` - Bayesian gravity inversion
#    * :doc:`../examples_probabilistic/05_magnetics_inversion` - Joint geophysical inversion
