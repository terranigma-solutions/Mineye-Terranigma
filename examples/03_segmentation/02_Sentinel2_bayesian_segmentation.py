"""
Bayesian Segmentation of Satellite Imagery
===========================================

This example demonstrates Bayesian segmentation of Sentinel-2 satellite imagery
using the BaySeg algorithm for geological mapping.

**Key concepts:**

* Multi-spectral satellite image analysis
* Bayesian inference for image segmentation
* Markov Random Fields (MRF) for spatial coherence
* Lithological unit classification

**Workflow:**

1. Load Sentinel-2 bands and crop to ROI
2. Initialize BaySeg segmenter
3. Run Bayesian segmentation
4. Visualize results and diagnostics
5. Analyze class distributions

.. note::
   This example requires Sentinel-2 data and the bayseg package.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import time
import os

try:
    from bayseg import BaySeg
    BAYSEG_AVAILABLE = True
except ImportError:
    BAYSEG_AVAILABLE = False
    print("⚠ BaySeg package not installed")
    print("Install with: pip install bayseg")

try:
    from mineye.BayesianSegmentation.full_workflow import run_workflow
    WORKFLOW_AVAILABLE = True
except ImportError:
    WORKFLOW_AVAILABLE = False
    print("⚠ mineye package not correctly installed or full_workflow not found")

# %%
# Configuration
# -------------
# Set segmentation parameters and data paths

# Robust path handling: determine project root
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd()

# If we are in 'examples/03_segmentation', project root is 2 levels up
if 'examples/03_segmentation' in script_dir:
    project_root = os.path.abspath(os.path.join(script_dir, "../../"))
else:
    project_root = os.getcwd()

# Data paths for Sentinel-2 bands
data_subdir = "examples/Data/Segmentation_Input_Data/Sentinel2/Tharsis_all_Sentinel2 /combined/20m"
bands = {
    "B4": os.path.join(project_root, data_subdir, "merged_B04_20m.jp2"),
    "B6": os.path.join(project_root, data_subdir, "merged_B06_20m.jp2"),
    "B7": os.path.join(project_root, data_subdir, "merged_B07_20m.jp2"),
    "B8A": os.path.join(project_root, data_subdir, "merged_B8A_20m.jp2"),
    "B11": os.path.join(project_root, data_subdir, "merged_B11_20m.jp2"),
    "B12": os.path.join(project_root, data_subdir, "merged_B12_20m.jp2"),
    "TCI": os.path.join(project_root, data_subdir, "merged_TCI_20m.jp2"),
    "SCL": os.path.join(project_root, data_subdir, "merged_SCL_20m.jp2"),
}

# ROI and segmentation parameters
bounds = (204498, 4170995, 227828, 4187076)
n_classes = 6
beta_init = 30.0
n_iterations = 10 #using only a few steps for demonstration. in a real scenario 400 should be enough
beta_jump_length = 0.1

output_dir = os.path.join(project_root, "examples/Data/Segmentation_Output_Data")
os.makedirs(output_dir, exist_ok=True)

print(f"Segmentation configuration:")
print(f"  Classes: {n_classes}")
print(f"  Beta (spatial smoothness): {beta_init}")
print(f"  Iterations: {n_iterations}")

# %%
# Run Full Workflow
# -----------------
# This will load, crop, mask, and segment the data

if WORKFLOW_AVAILABLE and BAYSEG_AVAILABLE:
    print("\nRunning Bayesian Segmentation workflow...")

    # We use run_workflow which handles data loading, cropping, soil masking, and running BaySeg
    # It also saves the results and generates diagnostic plots.
    run_workflow(
        bands=bands,
        bounds=bounds,
        n_classes=n_classes,
        beta_init=beta_init,
        n_iterations=n_iterations,
        beta_jump_length=beta_jump_length,
        use_soil_mask=True,
        save_npy=True,
        plot_tci=True,
        ref_band="B4",
        output_prefix=os.path.join(output_dir, "Sentinel2_BaySeg"),
    )

    # Load the results for further visualization in this script if desired
    # though run_workflow already generates some plots and saves results.
    # We construct the filename the same way run_workflow does.
    result_filename = f"Sentinel2_BaySeg_labels_n{n_classes}_betajump{beta_jump_length}.npy"
    result_path = os.path.join(output_dir, result_filename)

    try:
        labels = np.load(result_path)
        DATA_AVAILABLE = True
        print(f"✓ Loaded segmentation results from {result_path}")
    except FileNotFoundError:
        DATA_AVAILABLE = False
        print(f"⚠ Segmentation results not found at {result_path}")
else:
    DATA_AVAILABLE = False
    if not WORKFLOW_AVAILABLE:
        print("⚠ Workflow not available.")
    if not BAYSEG_AVAILABLE:
        print("⚠ BaySeg not available.")

# %%
# Visualize Segmentation Results
# -------------------------------

if DATA_AVAILABLE and BAYSEG_AVAILABLE:
    fig, ax = plt.subplots(figsize=(12, 10))

    # Use a distinct colormap for classes
    cmap = ListedColormap(plt.cm.tab10(range(n_classes)))

    im = ax.imshow(labels, cmap=cmap)
    ax.set_title(f"Bayesian Segmentation Result ({n_classes} classes)",
                 fontsize=14, fontweight='bold')
    ax.axis('off')

    # Add colorbar with class labels
    cbar = plt.colorbar(im, ax=ax, ticks=range(n_classes),
                        fraction=0.046, pad=0.04)
    cbar.set_label("Lithological Class", fontsize=12)

    plt.tight_layout()
    plt.show()

# %%
# Diagnostic Plots
# ----------------
#
# .. note::
#    Diagnostic plots (log-likelihood evolution and acceptance ratios)
#    are generated automatically by the ``run_workflow`` function.

# %%
# Class Distribution Analysis
# ---------------------------

if DATA_AVAILABLE and BAYSEG_AVAILABLE:
    unique, counts = np.unique(labels, return_counts=True)
    total_pixels = labels.size

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Bar plot
    colors = plt.cm.tab10(range(n_classes))
    ax1.bar(unique, counts, color=colors[:len(unique)])
    ax1.set_xlabel('Class', fontsize=12)
    ax1.set_ylabel('Pixel Count', fontsize=12)
    ax1.set_title('Class Distribution', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')

    # Pie chart
    ax2.pie(counts, labels=[f'Class {i}' for i in unique],
            colors=colors[:len(unique)],
            autopct='%1.1f%%', startangle=90)
    ax2.set_title('Class Proportions', fontsize=13, fontweight='bold')

    plt.tight_layout()
    plt.show()

    # Print statistics
    print("\nClass Distribution:")
    print("="*40)
    for class_id, count in zip(unique, counts):
        pct = 100 * count / total_pixels
        print(f"  Class {class_id}: {count:,} pixels ({pct:.2f}%)")
    print("="*40)
    print(f"  Total: {total_pixels:,} pixels")

# %%
# Summary
# -------
#
# This example demonstrated:
#
# * Bayesian segmentation of satellite imagery
# * MCMC-based inference with spatial priors
# * Visualization of segmentation results
# * Diagnostic analysis of convergence
#
# **Key advantages of BaySeg:**
#
# * Accounts for spatial coherence (MRF prior)
# * Provides uncertainty estimates
# * No need for labeled training data
# * Robust to noise and missing data
#
# **Next steps:**
#
# * Validate results with ground truth data
# * Integrate with geological models
# * Use for resource exploration
#
# See full workflow at: ``examples/Segmentation/``

if not DATA_AVAILABLE:
    print("\n" + "="*60)
    print("To run this example:")
    print("1. Install bayseg: pip install bayseg")
    print("2. Prepare data: run examples/Segmentation/prepare_data.py")
    print("3. Or see examples/Segmentation/run_segmentation.py")
    print("="*60)
elif not BAYSEG_AVAILABLE:
    print("\n" + "="*60)
    print("Install bayseg package:")
    print("  pip install bayseg")
    print("="*60)
