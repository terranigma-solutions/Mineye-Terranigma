"""
EnMap Likelihood and Model Comparison
=====================================

This example demonstrates how to compare a 3D geological model's predictions
with lithological labels extracted from EnMap hyperspectral data.

**Overview**

Once we have a 3D geological model and surface lithological information (from EnMap),
we can evaluate how well the model honors the surface observations. This comparison
is essential for:

1. **Model Validation**: Quantifying the accuracy of the geological interpretation at the surface.
2. **Likelihood Definition**: Defining a misfit function for probabilistic inversions.
3. **Residual Analysis**: Identifying areas where the geological model fails to explain surface data.

**Workflow**

1. Load the EnMap extracted points (see Example 02).
2. Set these points as a `custom_grid` in the GemPy model.
3. Compute the model to get predicted labels at these locations.
4. Map EnMap class IDs to GemPy lithology IDs.
5. Calculate accuracy and visualize residuals.
"""

# %%
# Import Libraries
# ----------------

import os
import numpy as np
import matplotlib.pyplot as plt
import gempy as gp
import gempy_viewer as gpv
from mineye.config import paths

# Set random seed for reproducibility
np.random.seed(1234)

# %%
# Load Model and Data
# -------------------
# We use the Tharsis geological model and the EnMap points extracted in the previous step.

# 1. Define Model Extent
extent = [-707521, -675558, 4526832, 4551949, -500, 505]

# 2. Get Data Paths
mod_or_path = paths.get_orientations_path()
mod_pts_path = paths.get_points_path()
topo_path = paths.get_topography_path()

# 3. Create GemPy Model
simple_geo_model = gp.create_geomodel(
    project_name='enmap_comparison',
    extent=extent,
    refinement=5,
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

# Set topography
gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topo_path)

# %%
# Load EnMap Extracted Data
# -------------------------
# For this example, we assume central points have been extracted.
# We define a helper function to extract these points from the EnMap results.

import rasterio
from rasterio.windows import from_bounds

def extract_points_central_reduced(raster_path, extent, min_distance=25, topo_path=None):
    """Extract points from the center of bodies using distance transform."""
    from skimage.segmentation import find_boundaries
    from skimage.feature import peak_local_max
    from scipy import ndimage

    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        
        data_mapped = data.copy()
        mask_nan = np.isnan(data)
        data_mapped[data_mapped == 3] = 0
        
        data_temp = data_mapped.copy()
        data_temp[mask_nan] = 255
        boundaries = find_boundaries(data_temp, mode='thick')
        
        dist_mask = ~boundaries & ~mask_nan
        dist_transform = ndimage.distance_transform_edt(dist_mask)
        
        unique_labels = np.unique(data_mapped)
        unique_labels = unique_labels[~np.isnan(unique_labels) & (unique_labels != 1)]
        
        all_ii, all_jj, all_labels = [], [], []
        for label_val in unique_labels:
            mask = (data_mapped == label_val)
            peaks = peak_local_max(dist_transform, min_distance=min_distance, labels=mask)
            if len(peaks) > 0:
                all_ii.extend(peaks[:, 0]); all_jj.extend(peaks[:, 1])
                all_labels.extend([label_val] * len(peaks))
        
        ii, jj = np.array(all_ii), np.array(all_jj)
        xs, ys = rasterio.transform.xy(transform, ii.tolist(), jj.tolist())
        
        if topo_path:
            with rasterio.open(topo_path) as topo_src:
                zs = np.array([val[0] for val in topo_src.sample(zip(xs, ys))])
        else:
            zs = np.full_like(xs, extent[5])
            
        return np.column_stack((xs, ys, zs)), np.array(all_labels)

base_dir = paths.get_base_dir()
enmap_path = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap', 'EPSG3857_EnMap_result_n4_betajump0.1.tif')

print("Extracting EnMap points for comparison...")
xyz_central, labels_enmap = extract_points_central_reduced(enmap_path, extent, min_distance=50, topo_path=topo_path)

print(f"Loaded {len(xyz_central)} points from EnMap extraction.")

# %%
# Compute Model on Custom Grid
# ----------------------------
# We set the EnMap point locations as a custom grid to evaluate the model exactly at those points.

# 1. Set custom grid
gp.set_custom_grid(simple_geo_model.grid, xyz_central)

# 2. Compute model
gp.compute_model(simple_geo_model)

# 3. Get GemPy predicted labels at custom grid points
# These are stored in solutions.raw_arrays.custom
labels_gempy = simple_geo_model.solutions.raw_arrays.custom.astype(int)

# %%
# Label Mapping and Accuracy
# --------------------------
# EnMap and GemPy use different ID systems. We must map them to compare results.
#
# **Automated Best Mapping**
#
# Since the class IDs from unsupervised segmentation don't necessarily match GemPy's
# lithology IDs, we use an automated "best mapping" logic that finds the permutation
# of IDs that maximizes the agreement between the two datasets.

from itertools import permutations

def find_best_mapping(observed_labels, predicted_labels):
    obs_unique = np.unique(observed_labels)
    pred_unique = np.unique(predicted_labels)
    
    best_acc = -1
    best_map = {}
    
    # Try all permutations of mapping observed labels to predicted labels
    for p in permutations(pred_unique, len(obs_unique)):
        mapping = dict(zip(obs_unique, p))
        mapped = np.vectorize(mapping.get)(observed_labels)
        acc = np.mean(mapped == predicted_labels)
        if acc > best_acc:
            best_acc = acc
            best_map = mapping
            
    return best_map, best_acc

best_mapping, best_accuracy = find_best_mapping(labels_enmap, labels_gempy)

mapped_enmap_labels = np.vectorize(best_mapping.get)(labels_enmap)

# Calculate residuals (where labels don't match)
residuals = (mapped_enmap_labels != labels_gempy)

print(f"Best automated mapping found: {best_mapping}")
print(f"Overall Accuracy: {best_accuracy:.2%}")

# %%
# **Performance Metrics: Confusion Matrix**
#
# A confusion matrix provides a detailed breakdown of which geological units are
# being misclassified.

from sklearn.metrics import confusion_matrix, classification_report
import seaborn as sns

cm = confusion_matrix(mapped_enmap_labels, labels_gempy)
plt.figure(figsize=(8, 6))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
plt.xlabel('GemPy (Predicted)')
plt.ylabel('EnMap (Observed)')
plt.title('Confusion Matrix: EnMap vs. GemPy')
plt.show()

print("\nClassification Report:")
print(classification_report(mapped_enmap_labels, labels_gempy))

# %%
# Visualization
# -------------

x, y = xyz_central[:, 0], xyz_central[:, 1]
fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True, sharey=True)

# Plot EnMap Labels (Mapped)
sc0 = axes[0].scatter(x, y, c=mapped_enmap_labels, cmap='viridis', s=20, edgecolors='k', linewidth=0.5)
axes[0].set_title('EnMap Labels (Mapped)')
plt.colorbar(sc0, ax=axes[0])

# Plot GemPy Labels
sc1 = axes[1].scatter(x, y, c=labels_gempy, cmap='viridis', s=20, edgecolors='k', linewidth=0.5)
axes[1].set_title('GemPy Predicted Labels')
plt.colorbar(sc1, ax=axes[1])

# Plot Residuals
sc2 = axes[2].scatter(x, y, c=residuals, cmap='Reds', s=20, edgecolors='k', linewidth=0.5)
axes[2].set_title(f'Residuals (Mismatches)\nAccuracy: {best_accuracy:.2%}')
plt.colorbar(sc2, ax=axes[2], label='1 = Mismatch')

for ax in axes:
    ax.set_aspect('equal')
    ax.set_xlabel('X (m)')
axes[0].set_ylabel('Y (m)')

plt.tight_layout()
plt.show()
