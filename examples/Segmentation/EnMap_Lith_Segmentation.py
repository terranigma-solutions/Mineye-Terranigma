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

from mineye.BayesianSegmentation.EnMap_feature_extraction import enmap_to_feature_stack, create_rasterio_layer
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

    Accepts feature values that are either 2D numpy arrays or rasterio DatasetReader objects.
    - If MNF_01..MNF_03 exist, an RGB composite is saved.
    - Up to `max_panels` single-band quicklooks are tiled and saved.
    """
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    def as_array(val):
        # Rasterio dataset -> read first band
        if hasattr(val, "read") and hasattr(val, "profile"):
            try:
                arr = val.read(1)
                # Respect nodata on quicklooks (our pipeline uses -9999.0 sentinel)
                nodata = getattr(val, "nodata", None)
                if nodata is None:
                    try:
                        nodata = val.profile.get("nodata")
                    except Exception:
                        nodata = None
                if nodata is not None:
                    arr = arr.astype(np.float32, copy=False)
                    arr = np.where(arr == nodata, np.nan, arr)
                return arr
            except Exception:
                pass
        return np.asarray(val)

    # 1) RGB composite from first three MNF components (if present)
    mnf_keys = [k for k in features.keys() if k.startswith("MNF_")]
    if len(mnf_keys) >= 3:
        # Sort to ensure MNF_01, MNF_02, MNF_03 order
        mnf_keys_sorted = sorted(mnf_keys)[:3]
        r = _percentile_scale(as_array(features[mnf_keys_sorted[0]]))
        g = _percentile_scale(as_array(features[mnf_keys_sorted[1]]))
        b = _percentile_scale(as_array(features[mnf_keys_sorted[2]]))
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
            img = _percentile_scale(as_array(features[k]))
            plt.imshow(img, cmap='viridis')
            plt.title(k)
            plt.axis('off')
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_feature_quicklooks.png", dpi=200)
        plt.close()


def _features_to_memory_datasets(features: Dict[str, np.ndarray], meta: dict):
    """Create in-memory rasterio datasets for each 2D feature layer.

    Accepts either 2D numpy arrays or already-wrapped rasterio DatasetReader objects.

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
    datasets: List[rasterio.io.DatasetReader] = []
    memfiles: List[MemoryFile] = []

    for name, val in features.items():
        # Case 1: already a rasterio dataset
        if hasattr(val, "read") and hasattr(val, "profile"):
            ds = val
            # Basic sanity checks on shape
            try:
                h, w = ds.height, ds.width
            except Exception:
                h, w = None, None
            if h is not None and w is not None:
                if (height is not None and h != height) or (width is not None and w != width):
                    raise ValueError(f"Feature '{name}' dataset size ({h},{w}) does not match meta ({height},{width}).")
            # Keep underlying MemoryFile alive if attached
            mf = getattr(ds, "_memfile", None)
            if isinstance(mf, MemoryFile):
                memfiles.append(mf)
            bands_dict[name] = ds
            datasets.append(ds)
            continue

        # Case 2: numpy array -> wrap into rasterio MemoryFile
        arr = np.asarray(val)
        if arr.ndim != 2:
            raise ValueError(f"Feature '{name}' must be a 2D array, got shape {arr.shape}")
        if arr.shape != (height, width):
            if arr.shape == (width, height):
                arr = arr.T
            else:
                raise ValueError(
                    f"Feature '{name}' has shape {arr.shape}, expected ({height}, {width})"
                )
        # Use existing helper to wrap and keep MemoryFile attached
        ds = create_rasterio_layer(arr.astype(np.float32, copy=False), meta)
        bands_dict[name] = ds
        datasets.append(ds)
        mf = getattr(ds, "_memfile", None)
        if isinstance(mf, MemoryFile):
            memfiles.append(mf)

    return bands_dict, datasets, memfiles


def _sanitize_feature_name(name: str) -> str:
    import re
    s = re.sub(r"[^A-Za-z0-9_]+", "_", str(name)).strip("_")
    return s or "layer"


def _write_features_to_files(features: Dict[str, np.ndarray], meta: dict, out_prefix: str,
                             preferred_driver: str = "JP2OpenJPEG") -> Dict[str, str]:
    """Save each feature layer to a single-band raster file and return dict of name -> path.

    Ensures all output rasters have the same (height, width) as provided in meta.
    If an input array is accidentally transposed (width, height), it will be
    auto-transposed to (height, width). Any other mismatch raises a ValueError.

    Tries to write JPEG2000 (JP2OpenJPEG). If the driver is unavailable or writing fails,
    falls back to GeoTIFF. Preserves CRS/transform and writes float32.
    """
    os.makedirs(os.path.dirname(out_prefix), exist_ok=True)

    files: Dict[str, str] = {}

    # Expected size from meta
    expected_h = int(meta.get("height"))
    expected_w = int(meta.get("width"))
    if not expected_h or not expected_w:
        raise ValueError("meta must include valid 'height' and 'width' for writing features")

    # Base profile from meta (fixed size for all outputs)
    sentinel = -9999.0
    base_profile = {
        "dtype": "float32",
        "count": 1,
        "crs": meta.get("crs"),
        "transform": meta.get("transform"),
        "height": expected_h,
        "width": expected_w,
        "nodata": sentinel,
    }

    for name, val in features.items():
        safe = _sanitize_feature_name(name)
        # Prepare data as 2D float32 array
        if hasattr(val, "read") and hasattr(val, "profile"):
            data = val.read(1).astype(np.float32, copy=False)
        else:
            arr = np.asarray(val)
            if arr.ndim == 3 and arr.shape[0] == 1:
                arr = arr[0]
            if arr.ndim != 2:
                raise ValueError(f"Feature '{name}' must be 2D for writing, got shape {arr.shape}")
            data = arr.astype(np.float32, copy=False)

        # Enforce expected shape (H, W)
        if data.shape == (expected_w, expected_h):
            # likely transposed; fix
            data = data.T
            print(f"[EnMap] Note: transposed feature '{name}' from (W,H) to (H,W) to match meta.")
        elif data.shape != (expected_h, expected_w):
            raise ValueError(
                f"Feature '{name}' has shape {data.shape}, expected ({expected_h}, {expected_w})."
            )

        # Attempt JP2 first
        out_path_jp2 = f"{out_prefix}_{safe}.jp2"
        out_path_tif = f"{out_prefix}_{safe}.tif"

        profile = base_profile.copy()
        profile.update({"driver": preferred_driver})

        wrote_path = None
        try:
            with rasterio.open(out_path_jp2, "w", **profile) as dst:
                dst.write(data, 1)
            wrote_path = out_path_jp2
        except Exception as e:
            print(f"[WARN] JP2 write failed for '{name}' with driver '{preferred_driver}': {e}. Falling back to GeoTIFF.")
            # Fallback to GeoTIFF
            profile_gtiff = base_profile.copy()
            profile_gtiff.update({"driver": "GTiff"})
            with rasterio.open(out_path_tif, "w", **profile_gtiff) as dst:
                dst.write(data, 1)
            wrote_path = out_path_tif

        # Post-write validation: reopen and confirm shape
        try:
            with rasterio.open(wrote_path) as chk:
                if (chk.height, chk.width) != (expected_h, expected_w):
                    raise ValueError(
                        f"Wrote feature '{name}' to {wrote_path} with wrong size "
                        f"({chk.height},{chk.width}), expected ({expected_h},{expected_w})."
                    )
        except Exception as e:
            # If validation fails, ensure we surface a clear error
            raise

        files[name] = wrote_path

    print(f"[EnMap] Saved {len(files)} feature layers to files with prefix '{out_prefix}_*.(jp2|tif)', shape=({expected_h},{expected_w}).")
    return files


# %%
# Parameters
# ----------
# Set these paths to your local data. Wavelengths are parsed from the ENMAP metadata XML automatically.

# Example: enmap_folder = "/path/to/ENMAP_L2A_TILE/"

enmap_folder: str = "/Users/simonvirgo/Downloads/Enmap_data/ENMAP01-____L2A-DT0000026661_20230712T114038Z_001_V010402_20240818T134118Z"  # REQUIRED: EnMAP L2A folder path

# Segmentation hyperparameters
n_classes: int = 8
iterations: int = 400
beta_init: float = 30.0
beta_jump: float = 0.1
output_prefix: str = "examples/Segmentation/EnMap"  # output prefix for results
save_npy: bool = True  # save intermediate arrays from full_workflow

# Detrending / background-field removal (recommended for the current EnMAP dataset)
# The dominant artifact here is typically a smooth 2D cross-track bias (not narrow detector stripes).
detrend: bool = True
detrend_sigma: float = 75.0  # pixels; try 50–100
detrend_downsample: int | None = 4  # 2/4 for speed; set None for full-res

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

features, meta = enmap_to_feature_stack(
    enmap_folder,
    n_mnf=12,
    n_deriv_pcs=3,
    detrend=detrend,
    detrend_sigma=detrend_sigma,
    detrend_downsample=detrend_downsample,
    destripe=destripe,
    destripe_frac=destripe_frac,
    destripe_poly=destripe_poly,
    destripe_reference_kernel=destripe_reference_kernel,
    destripe_reference_downsample=destripe_reference_downsample,
    destripe_smooth_cols=destripe_smooth_cols,
)
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

# Persist features to disk as JP2 (fallback to GeoTIFF) and feed paths into workflow
paths_dict = _write_features_to_files(bands_dict, meta, output_prefix)
print(f"Prepared {len(paths_dict)} feature file paths for workflow input.")

# %%
# Run Bayesian Segmentation
# -------------------------
# Feed the feature rasters into the general segmentation workflow.

try:
    run_workflow(
        bands=paths_dict,  # dict of file paths (JP2 or TIF)
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
