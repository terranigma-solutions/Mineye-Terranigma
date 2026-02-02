"""
EnMap Lithological Segmentation
==============================

This example demonstrates the workflow for lithological segmentation using EnMap hyperspectral data.
It covers feature extraction, data preparation, and the Bayesian segmentation process.

**Overview**

EnMap (Environmental Mapping and Analysis Program) is a German hyperspectral satellite mission.
This example shows how to process EnMap L2A products to identify different lithological units.

**Key Steps**

1.  **Feature Extraction**: Extracting meaningful features from the hyperspectral cube (MNF, absorption depths, derivative PCs).
2.  **Preprocessing**: Detrending and background-field removal to account for sensor artifacts.
3.  **Bayesian Segmentation**: Running the core segmentation engine to classify pixels and estimate uncertainty.
4.  **Export**: Saving the results as georeferenced GeoTIFFs.

.. note::
   This example requires an EnMap L2A product folder. Update the `enmap_folder` path to point to your local data.
"""

# %%
# Import Libraries
# ----------------
#
# We import the necessary modules for EnMap feature extraction and the core segmentation workflow.

import os
from mineye.BayesianSegmentation.EnMap_feature_extraction import (
    enmap_to_feature_stack,
    plot_feature_quicklooks,
    save_detrended_swir_column_median_plot,
    features_to_memory_datasets,
    write_features_to_files,
)
from mineye.BayesianSegmentation.full_workflow import run_workflow

# %%
# Parameters
# ----------
#
# Define the input paths and hyperparameters for the segmentation.
# The `enmap_folder` should point to a valid EnMAP L2A directory containing `SPECTRAL_IMAGE.TIF` and metadata.

# Attempt to get paths from Sphinx config if running in gallery
import os
try:
    from sphinx.util import logging
    # This is a bit of a hack to check if we are in a sphinx build
    # In a real gallery build, we might use environment variables
    enmap_root = os.getenv('ENMAP_DATA_ROOT', 'examples/Data/Segmentation_Input_Data/Enmap')
    output_root = os.getenv('SEGMENTATION_OUTPUT_ROOT', 'examples/Data/Segmentation_Output_Data')
except ImportError:
    enmap_root = 'examples/Data/Segmentation_Input_Data/Enmap'
    output_root = 'examples/Data/Segmentation_Output_Data'

# Paths
enmap_folder: str = os.path.join(enmap_root, "ENMAP01-____L2A-DT0000026661_20230712T114038Z_001_V010402_20240818T134118Z")
output_prefix: str = os.path.join(output_root, "EnMap")

# If the path is relative, make it relative to the project root
if not os.path.isabs(enmap_folder):
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # __file__ is not defined when running in some environments (like Sphinx Gallery exec)
        # In Sphinx Gallery, the CWD is typically the directory where the script is located
        script_dir = os.getcwd()
    
    # Sphinx Gallery execution check: if we are in 'docs/source/examples_segmentation'
    # or similar, the project root is 3 levels up. Otherwise, assume 2 levels up.
    if 'docs/source' in script_dir:
        project_root = os.path.abspath(os.path.join(script_dir, "../../../"))
    else:
        project_root = os.path.abspath(os.path.join(script_dir, "../../"))
    
    enmap_folder = os.path.join(project_root, enmap_folder)
    output_prefix = os.path.join(project_root, output_prefix)

# Segmentation hyperparameters
n_classes: int = 8
iterations: int = 5
beta_init: float = 30.0
beta_jump: float = 0.1
save_npy: bool = True

# Detrending / background-field removal settings
detrend: bool = True
detrend_sigma: float = 120.0
detrend_downsample: int | None = 4
detrend_buffer_px: int = 20
detrend_den_thresh: float = 0.2

# Additional exclusions and handling
edge_px: int = 40
mnf_buffer_px: int | None = None
destripe: bool = False

# Workflow options
use_soil_mask: bool = False
plot_tci: bool = False

# %%
# Feature Extraction
# ------------------
#
# We extract a stack of descriptive features from the EnMap cube.
# This includes Minimum Noise Fraction (MNF) transforms, specific absorption depth maps, and derivative-based Principal Components.

if not os.path.isdir(enmap_folder):
    print(f"EnMAP folder not found at {enmap_folder}. Skipping extraction for demonstration.")
    features, meta = None, None
else:
    features, meta = enmap_to_feature_stack(
        enmap_folder,
        detrend=detrend,
        detrend_sigma=detrend_sigma,
        detrend_downsample=detrend_downsample,
        detrend_buffer_px=detrend_buffer_px,
        detrend_den_thresh=detrend_den_thresh,
        mnf_buffer_px=mnf_buffer_px,
        edge_px=edge_px,
        destripe=destripe,
    )
    print(f"Extracted {len(features)} feature layers")

# %%
# QA Plotting
# -----------
#
# Visualize the extracted features and diagnostic plots to ensure the preprocessing (like detrending) worked as expected.

if features:
    output_prefix_abs = os.path.abspath(output_prefix)
    try:
        # Pass show=True to ensure Sphinx Gallery captures the plots
        print("[Workflow] Plotting feature quicklooks...")
        plot_feature_quicklooks(features, output_prefix_abs, show=True)
    except Exception as e:
        print(f"[WARN] Could not create feature quicklooks: {e}")

    try:
        diag_path = save_detrended_swir_column_median_plot(
            enmap_folder=enmap_folder,
            out_prefix=output_prefix_abs,
            target_nm=2200.0,
            use_qa_mask=True,
            edge_px=edge_px,
            detrend_sigma=detrend_sigma,
            detrend_downsample=detrend_downsample,
            detrend_buffer_px=detrend_buffer_px,
            detrend_den_thresh=detrend_den_thresh,
        )
        print(f"Saved diagnostic plot to {diag_path}")
    except Exception as e:
        print(f"Could not create diagnostic plot: {e}")

# %%
# Data Preparation
# ----------------
#
# Wrap the extracted features into in-memory datasets and persist them to disk for the segmentation engine.

if features and meta:
    bands_dict, open_datasets, memfiles = features_to_memory_datasets(features, meta)
    
    ref_key = next(iter(bands_dict.keys()))
    ref_ds = bands_dict[ref_key]
    ref_bounds = ref_ds.bounds
    bounds = (
        float(ref_bounds.left), float(ref_bounds.bottom), float(ref_bounds.right), float(ref_bounds.top)
    )
    
    paths_dict = write_features_to_files(bands_dict, meta, output_prefix)

    # %%
    # Bayesian Segmentation
    # ---------------------
    #
    # Run the Markov Chain Monte Carlo (MCMC) based Bayesian segmentation.
    # This process iteratively updates the class labels and calculates the entropy (uncertainty).

    try:
        run_workflow(
            bands=paths_dict,
            bounds=bounds,
            n_classes=n_classes,
            beta_init=beta_init,
            n_iterations=iterations,
            beta_jump_length=beta_jump,
            use_soil_mask=use_soil_mask,
            save_npy=save_npy,
            plot_tci=plot_tci,
            ref_band=ref_key,
            output_prefix=output_prefix,
        )
    finally:
        # Cleanup open file handles
        for ds in open_datasets:
            try:
                ds.close()
            except:
                pass
        for mf in memfiles:
            try:
                mf.close()
            except:
                pass
else:
    print("Skipping segmentation because features were not extracted.")
