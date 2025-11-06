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

1. Load preprocessed Sentinel-2 bands
2. Initialize BaySeg segmenter
3. Run Bayesian segmentation
4. Visualize results and diagnostics
5. Analyze class distributions

.. note::
   This example requires preprocessed Sentinel-2 data and the bayseg package.
"""

# %%
# Import Libraries
# ----------------

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import time

try:
    from bayseg import BaySeg
    BAYSEG_AVAILABLE = True
except ImportError:
    BAYSEG_AVAILABLE = False
    print("⚠ BaySeg package not installed")
    print("Install with: pip install bayseg")

# %%
# Load Preprocessed Data
# ----------------------
# Load the stacked Sentinel-2 bands

try:
    layers = np.load("sentinel2_bayseg_input.npy")
    DATA_AVAILABLE = True
    print(f"✓ Loaded satellite data")
    print(f"  Shape: {layers.shape}")
    print(f"  (rows, cols, bands): {layers.shape}")
except FileNotFoundError:
    DATA_AVAILABLE = False
    print("⚠ Data file not found: sentinel2_bayseg_input.npy")
    print("  Run prepare_data.py first to create input data")

# %%
# Configuration
# -------------
# Set segmentation parameters

n_classes = 6  # Number of lithological classes
beta_init = 10.0  # Initial spatial smoothness parameter
n_iterations = 200  # Number of MCMC iterations
beta_jump_length = 0.1  # Step size for beta updates

print(f"Segmentation configuration:")
print(f"  Classes: {n_classes}")
print(f"  Beta (spatial smoothness): {beta_init}")
print(f"  Iterations: {n_iterations}")

# %%
# Initialize Bayesian Segmenter
# -----------------------------

if DATA_AVAILABLE and BAYSEG_AVAILABLE:
    print("\nInitializing BaySeg...")
    start_time = time.time()

    seg = BaySeg(
        data=layers,
        n_labels=n_classes,
        beta_init=beta_init
    )

    print(f"✓ Segmenter initialized ({time.time() - start_time:.2f}s)")

# %%
# Run Segmentation
# ----------------
# Perform Bayesian segmentation using MCMC

if DATA_AVAILABLE and BAYSEG_AVAILABLE:
    print(f"\nRunning segmentation ({n_iterations} iterations)...")
    start_time = time.time()

    labels = seg.fit(n=n_iterations, beta_jump_length=beta_jump_length)

    elapsed = time.time() - start_time
    print(f"✓ Segmentation complete ({elapsed:.2f}s)")
    print(f"  Speed: {elapsed/n_iterations*1000:.1f} ms/iteration")

    # Save results
    np.save("bayseg_lithology_labels.npy", labels)
    print("✓ Results saved to bayseg_lithology_labels.npy")

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
# Show MCMC convergence diagnostics

if DATA_AVAILABLE and BAYSEG_AVAILABLE:
    # Plot log-likelihood evolution
    seg.diagnostics_plot(transpose=False)

    # Plot acceptance ratios for MCMC proposals
    seg.plot_acc_ratios()

    print("\n✓ Diagnostic plots generated")

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
