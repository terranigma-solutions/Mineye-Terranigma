"""
EnMap Lithological Segmentation Example
===================================================

This example demonstrates how to:

- Load an EnMAP L2A product
- Extract a stack of descriptive features (MNF, absorption depths, derivative PCs)
- Wrap each feature as an in-memory, georeferenced raster
- Run the Bayesian segmentation workflow on the feature stack
- Export georeferenced label and entropy GeoTIFFs

How to use
----------
- Set `enmap_folder` to your EnMAP L2A product directory (must contain SPECTRAL_IMAGE.TIF and QA layers)
- Wavelengths are parsed automatically from the EnMAP metadata XML; no external wavelengths file is needed.
"""

# %%
# Import Libraries
# ----------------

from mineye.BayesianSegmentation.EnMap_feature_extraction import (
    enmap_to_feature_stack,
    plot_feature_quicklooks,
    plot_feature_histograms,
    save_detrended_swir_column_median_plot,
    features_to_memory_datasets,
    write_features_to_files,
)
from mineye.BayesianSegmentation.full_workflow import run_workflow
import os

# %%
# Parameters
# ----------
# Set these paths to your local data. Wavelengths are parsed from the ENMAP metadata XML automatically.

# Examples: enmap_folder = "/path/to/ENMAP_L2A_TILE/"

enmap_folder: str = "/Users/simonvirgo/Downloads/Enmap_data/ENMAP01-____L2A-DT0000026661_20230712T114038Z_001_V010402_20240818T134118Z"  # REQUIRED: EnMAP L2A folder path

# Segmentation hyperparameters
n_classes: int = 8
iterations: int = 400
beta_init: float = 30.0
beta_jump: float = 0.1
output_prefix: str = "examples/Data/Segmentation_Input_Data/Segmentation_Output_Data/EnMap"  # output prefix for results
save_npy: bool = True  # save intermediate arrays from full_workflow

# Detrending / background-field removal (recommended for the current EnMAP dataset)
# The dominant artifact here is typically a smooth 2D cross-track bias (not narrow detector stripes).
detrend: bool = True
detrend_sigma: float = 120.0  # pixels; larger sigma removes broader cross-track fields
detrend_downsample: int | None = 4  # 2/4 for speed; set None for full-res
detrend_buffer_px: int = 20  # erode footprint to avoid rim-driven background estimates
detrend_den_thresh: float = 0.2  # feather correction where G(V) is small (fraction of den.max)

# Additional hard edge exclusion (array border), independent of QA and footprint erosion.
edge_px: int = 40

# MNF rim handling: fit MNF/noise on an eroded interior footprint.
# Default: reuse the detrending buffer.
mnf_buffer_px: int | None = None

# Destriping (keep available for datasets that truly have detector striping)
destripe: bool = False
destripe_frac: float = 1.0

# Residual-based destriping settings (only used if destripe=True)
destripe_poly: int | None = None
destripe_reference_kernel: int | None = 51
destripe_smooth_cols: int | None = 21
destripe_reference_downsample: int | None = 4

# Derived features are not raw spectral bands, so we disable soil mask and TCI plotting
use_soil_mask: bool = False
plot_tci: bool = False

# %%
# Run Workflow
# ------------

if __name__ == "__main__":
    if not enmap_folder or not os.path.isdir(enmap_folder):
        raise FileNotFoundError(
            f"Please set 'enmap_folder' to a valid EnMAP L2A product directory. Got: {enmap_folder}"
        )

    print(f"[Workflow] ENMAP folder: {enmap_folder}")

    # 1. Feature Extraction
    # ---------------------
    # Extracts MNF, absorption depths, and derivative PCs from the EnMAP cube.
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
        destripe_frac=destripe_frac,
        destripe_poly=destripe_poly,
        destripe_reference_kernel=destripe_reference_kernel,
        destripe_reference_downsample=destripe_reference_downsample,
        destripe_smooth_cols=destripe_smooth_cols,
    )
    print(f"[Workflow] Extracted {len(features)} feature layers")

    # 2. QA plotting
    # --------------
    # Visualize the extracted features and diagnostic plots.
    try:
        output_prefix_abs = os.path.abspath(output_prefix)
        plot_feature_quicklooks(features, output_prefix_abs)
        print(f"[Workflow] Saved feature quicklooks to {output_prefix_abs}_*.png")
    except Exception as e:
        print(f"[WARN] Could not create feature quicklooks: {e}")

    try:
        output_prefix_abs = os.path.abspath(output_prefix)
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
        print(f"[Workflow] Saved detrended SWIR column-median diagnostic to {diag_path}")
    except Exception as e:
        print(f"[WARN] Could not create SWIR column-median diagnostic: {e}")

    # 3. Wrap and save features
    # -------------------------
    # Prepare in-memory datasets for the segmentation engine and persist them to disk.
    bands_dict, open_datasets, memfiles = features_to_memory_datasets(features, meta)
    try:
        # Choose a reference layer and compute full bounds (xmin, ymin, xmax, ymax)
        ref_key = next(iter(bands_dict.keys()))
        ref_ds = bands_dict[ref_key]
        ref_bounds = ref_ds.bounds
        bounds = (
            float(ref_bounds.left), float(ref_bounds.bottom), float(ref_bounds.right), float(ref_bounds.top)
        )

        print(f"[Workflow] Reference layer: {ref_key}")
        print(f"[Workflow] Bounds -> xmin={bounds[0]:.3f}, ymin={bounds[1]:.3f}, xmax={bounds[2]:.3f}, ymax={bounds[3]:.3f}")

        # Persist features to disk
        paths_dict = write_features_to_files(bands_dict, meta, output_prefix)

        # 4. Run Bayesian Segmentation
        # ----------------------------
        # Run the actual segmentation engine.
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
        # 5. Cleanup
        # ----------
        # Ensure all open file handles are closed.
        for ds in open_datasets:
            try:
                ds.close()
            except Exception:
                pass
        for mf in memfiles:
            try:
                mf.close()
            except Exception:
                pass

