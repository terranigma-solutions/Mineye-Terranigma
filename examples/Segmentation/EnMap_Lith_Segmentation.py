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

enmap_folder: str = ""  # REQUIRED: EnMAP L2A folder path

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
