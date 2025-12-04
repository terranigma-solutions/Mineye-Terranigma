"""
EnMap Lithological Segmentation (cell-based example)
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

This example follows a linear, cell-based structure similar to the gravity forward model example.
"""

# %%
# Import Libraries
# ----------------

from __future__ import annotations

import os
from typing import Dict, List, Tuple

import numpy as np
import rasterio
from rasterio.io import MemoryFile

from mineye.BayesianSegmentation.EnMap_feature_extraction import enmap_to_feature_stack
from mineye.BayesianSegmentation.full_workflow import run_workflow


# %%
# Helper Functions
# ----------------

import matplotlib.pyplot as plt

def _percentile_scale(arr, p_lo=2.0, p_hi=98.0):
    a = np.asarray(arr, dtype=float)
    vmin, vmax = np.nanpercentile(a, [p_lo, p_hi])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = np.nanmin(a), np.nanmax(a)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = 0.0, 1.0
    a = (a - vmin) / (vmax - vmin)
    return np.clip(a, 0, 1)


def plot_feature_quicklooks(features: Dict[str, np.ndarray], out_prefix: str, max_panels: int = 9):
    """Create quicklook plots of feature bands for QA.

    - If MNF_01..MNF_03 exist, an RGB composite is saved.
    - Up to `max_panels` single-band quicklooks are tiled and saved.
    """
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    # 1) RGB composite from first three MNF components (if present)
    mnf_keys = [k for k in features.keys() if k.startswith("MNF_")]
    if len(mnf_keys) >= 3:
        # Sort to ensure MNF_01, MNF_02, MNF_03 order
        mnf_keys_sorted = sorted(mnf_keys)[:3]
        r = _percentile_scale(features[mnf_keys_sorted[0]])
        g = _percentile_scale(features[mnf_keys_sorted[1]])
        b = _percentile_scale(features[mnf_keys_sorted[2]])
        rgb = np.dstack([r, g, b])
        plt.figure(figsize=(6, 6))
        plt.imshow(rgb)
        plt.title(f"RGB composite: {mnf_keys_sorted[0]}, {mnf_keys_sorted[1]}, {mnf_keys_sorted[2]}")
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_MNF_RGB.png", dpi=200)
        plt.close()

    # 2) Tile up to max_panels single-band quicklooks
    keys = list(features.keys())
    keys.sort()
    sel = keys[:max_panels]
    if sel:
        n = len(sel)
        cols = min(3, n)
        rows = int(np.ceil(n / cols))
        plt.figure(figsize=(4 * cols, 4 * rows))
        for i, k in enumerate(sel, 1):
            plt.subplot(rows, cols, i)
            img = _percentile_scale(features[k])
            plt.imshow(img, cmap='viridis')
            plt.title(k)
            plt.axis('off')
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_feature_quicklooks.png", dpi=200)
        plt.close()


def _features_to_memory_datasets(features: Dict[str, np.ndarray], meta: dict):
    """Create in-memory rasterio datasets for each 2D feature layer.

    Returns:
        bands_dict: dict[str, rasterio DatasetReader]
        datasets: list of open dataset handles (to close later)
        memfiles: list of MemoryFile objects (to close later)
    """
    transform = meta.get("transform")
    crs = meta.get("crs")
    height = meta.get("height")
    width = meta.get("width")

    bands_dict: Dict[str, rasterio.io.DatasetReader] = {}
    datasets: List[rasterio.io.DatasetWriter] = []
    memfiles: List[MemoryFile] = []

    for name, arr in features.items():
        if arr.ndim != 2:
            raise ValueError(f"Feature '{name}' must be a 2D array, got shape {arr.shape}")
        if arr.shape != (height, width):
            if arr.shape == (width, height):
                arr = arr.T
            else:
                raise ValueError(
                    f"Feature '{name}' has shape {arr.shape}, expected ({height}, {width})"
                )

        data = np.asarray(arr, dtype=np.float32)
        nodata = -9999.0
        data = np.where(np.isfinite(data), data, nodata)

        mem = MemoryFile()
        profile = {
            "driver": "GTiff",
            "dtype": "float32",
            "count": 1,
            "height": height,
            "width": width,
            "transform": transform,
            "crs": crs,
            "nodata": nodata,
        }
        ds = mem.open(**profile)
        ds.write(data, 1)

        bands_dict[name] = ds
        datasets.append(ds)
        memfiles.append(mem)

    return bands_dict, datasets, memfiles


# %%
# Parameters
# ----------
# Set these paths to your local data. Wavelengths are parsed from the ENMAP metadata XML automatically.

# Example: enmap_folder = "/path/to/ENMAP_L2A_TILE/"

enmap_folder: str = "/Users/simonvirgo/Downloads/Enmap_data/ENMAP01-____L2A-DT0000026661_20230712T114038Z_001_V010402_20240818T134118Z"  # REQUIRED: EnMAP L2A folder path

# Segmentation hyperparameters
n_classes: int = 6
iterations: int = 400
beta_init: float = 30.0
beta_jump: float = 0.1
output_prefix: str = "examples/Segmentation/EnMap"  # output prefix for results
save_npy: bool = True  # save intermediate arrays from full_workflow

# Derived features are not raw spectral bands, so we disable soil mask and TCI plotting
use_soil_mask: bool = False
plot_tci: bool = False


# %%
# Validate Inputs
# ---------------

if not enmap_folder or not os.path.isdir(enmap_folder):
    raise FileNotFoundError(
        "Please set 'enmap_folder' to a valid EnMAP L2A product directory containing SPECTRAL_IMAGE.TIF."
    )

print(f"ENMAP folder: {enmap_folder}")


# %%
# Feature Extraction
# ------------------
# Compute a stack of features from the EnMAP spectral cube.

features, meta = enmap_to_feature_stack(enmap_folder, n_mnf=12, n_deriv_pcs=3)
print(f"Extracted {len(features)} feature layers")

# QA quicklooks of feature bands (saves PNGs alongside output_prefix)
try:
    output_prefix_abs = os.path.abspath(output_prefix)
    plot_feature_quicklooks(features, output_prefix_abs)
    print(f"Saved feature quicklooks to {output_prefix_abs}_*.png")
except Exception as e:
    print(f"[WARN] Could not create feature quicklooks: {e}")


# %%
# Wrap Features as In-Memory Rasters
# ----------------------------------
# Each feature becomes a georeferenced in-memory rasterio dataset.

bands_dict, open_datasets, memfiles = _features_to_memory_datasets(features, meta)

# Choose a reference layer and compute full bounds (xmin, ymin, xmax, ymax)
ref_key = next(iter(bands_dict.keys()))
ref_ds = bands_dict[ref_key]
ref_bounds = ref_ds.bounds  # rasterio.coords.BoundingBox
bounds: Tuple[float, float, float, float] = (
    float(ref_bounds.left), float(ref_bounds.bottom), float(ref_bounds.right), float(ref_bounds.top)
)

print(f"Reference layer: {ref_key}")
print(
    f"Bounds -> xmin={bounds[0]:.3f}, ymin={bounds[1]:.3f}, xmax={bounds[2]:.3f}, ymax={bounds[3]:.3f}"
)


# %%
# Run Bayesian Segmentation
# -------------------------
# Feed the feature rasters into the general segmentation workflow.

try:
    run_workflow(
        bands=bands_dict,  # dict of open rasterio datasets
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
    # Ensure all in-memory datasets are closed to free resources
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

# %%
# Summary
# -------
# - This example extracted features from EnMAP reflectance, wrapped them into in-memory rasters, and
#   passed them to the generic Bayesian segmentation workflow.
# - The outputs include georeferenced GeoTIFFs for labels and entropy with the provided output prefix.
