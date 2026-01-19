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

from mineye.BayesianSegmentation.EnMap_feature_extraction import (
    enmap_to_feature_stack,
    create_rasterio_layer,
    parse_enmap_wavelengths,
    build_mask,
    remove_background_field,
)
from mineye.BayesianSegmentation.full_workflow import run_workflow


# %%
# Helper Functions
# ----------------

import matplotlib.pyplot as plt


def save_detrended_swir_column_median_plot(
        enmap_folder: str,
        out_prefix: str,
        target_nm: float = 2200.0,
        use_qa_mask: bool = True,
        edge_px: int = 0,
        detrend_sigma: float = 80.0,
        detrend_downsample: int | None = 4,
        detrend_buffer_px: int = 20,
        detrend_den_thresh: float = 0.2,
) -> str:
    """Save a diagnostic plot of the per-column median after detrending for a SWIR band.

    This is a lightweight check to see whether any systematic cross-track structure remains.
    The plot is saved as a PNG (no interactive window required).

    Returns:
        Path to the written PNG.
    """
    files = os.listdir(enmap_folder)

    cube_name = next(f for f in files if "SPECTRAL_IMAGE" in f and f.upper().endswith((".TIF", ".TIFF")))
    xml_name = next(f for f in files if "METADATA.XML" in f.upper())
    cube_path = os.path.join(enmap_folder, cube_name)
    xml_path = os.path.join(enmap_folder, xml_name)

    wls = parse_enmap_wavelengths(xml_path)
    if wls.size == 0:
        raise ValueError("No wavelengths parsed from metadata XML; cannot select SWIR band.")

    band0 = int(np.nanargmin(np.abs(wls - float(target_nm))))
    wl = float(wls[band0])

    mask_bad = None
    if use_qa_mask:
        qa_masks = {}
        for f in files:
            fu = f.upper()
            if ("QL_PIXELMASK" in fu or "QL_QUALITY" in fu) and fu.endswith((".TIF", ".TIFF")):
                fp = os.path.join(enmap_folder, f)
                with rasterio.open(fp) as src:
                    qa_masks[f] = src.read(1)
        if qa_masks:
            mask_bad = build_mask(qa_masks)

    with rasterio.open(cube_path) as src:
        img = src.read(band0 + 1).astype(np.float32, copy=False)

    # Mirror the scaling heuristic used in `load_enmap_dataset`.
    mx = float(np.nanmax(img))
    mn = float(np.nanmin(img))
    if 2.0 < mx <= 10000.0 and mn >= 0.0:
        img = img / 10000.0

    # Optional: exclude a fixed-width rim at array borders (helps prevent edge artifacts from
    # dominating this diagnostic on rotated/partial footprints).
    if edge_px is not None and int(edge_px) > 0:
        ep = int(edge_px)
        rows, cols = img.shape
        y = np.minimum(np.arange(rows), np.arange(rows)[::-1])
        x = np.minimum(np.arange(cols), np.arange(cols)[::-1])
        Y, X = np.meshgrid(y, x, indexing="ij")
        edge_valid = (np.minimum(Y, X) >= ep)
        edge_bad = ~edge_valid
        if mask_bad is None:
            mask_bad = edge_bad
        else:
            mask_bad = (mask_bad | edge_bad)

    img_dt = remove_background_field(
        img[None, :, :],
        mask_bad=mask_bad,
        sigma=float(detrend_sigma),
        downsample=detrend_downsample,
        buffer_px=int(detrend_buffer_px),
        den_thresh=float(detrend_den_thresh),
        progress_every=None,
    )[0]

    if mask_bad is not None:
        img_dt = np.where(mask_bad, np.nan, img_dt)

    col = np.nanmedian(img_dt, axis=0)

    out_path = f"{out_prefix}_column_median_after_detrend_wl{wl:.1f}nm_band{band0 + 1:03d}.png"
    plt.figure(figsize=(10, 3))
    plt.plot(col)
    plt.title(f"column median (after detrending) — {wl:.1f} nm (band {band0 + 1})")
    plt.xlabel("column")
    plt.ylabel("nanmedian")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    return out_path

def _percentile_scale(arr, p_lo=2.0, p_hi=98.0):
    a = np.asarray(arr, dtype=float)
    vmin, vmax = np.nanpercentile(a, [p_lo, p_hi])
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = np.nanmin(a), np.nanmax(a)
        if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
            vmin, vmax = 0.0, 1.0
    a = (a - vmin) / (vmax - vmin)
    return np.clip(a, 0, 1)


def plot_feature_quicklooks(
        features: Dict[str, np.ndarray],
        out_prefix: str,
        max_panels: int = 9,
        panels_per_page: int = 12,
):
    """Create quicklook plots of feature layers for QA.

    Accepts feature values that are either 2D numpy arrays or rasterio DatasetReader objects.

    Outputs:
    - If `MNF_01..MNF_03` exist, an RGB composite is saved.
    - Grouped quicklook pages for *all* `MNF_*`, `Depth_*`, and `Deriv_*` layers.
    - A small mixed overview (up to `max_panels` layers across all keys) for convenience.
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

    def _plot_tiles(keys: List[str], out_path: str, title: str | None = None):
        if not keys:
            return
        n = len(keys)
        cols = min(3, n)
        rows = int(np.ceil(n / cols))
        plt.figure(figsize=(4 * cols, 4 * rows))
        for i, k in enumerate(keys, 1):
            plt.subplot(rows, cols, i)
            arr = as_array(features[k])

            if k.startswith("Depth_"):
                img = _percentile_scale(arr, 10.0, 99.5)  # tighter stretch for QA
            else:
                img = _percentile_scale(arr, 2.0, 98.0)
            plt.imshow(img, cmap="viridis")
            plt.title(k)
            plt.axis("off")
        if title:
            plt.suptitle(title)
        plt.tight_layout()
        plt.savefig(out_path, dpi=200)
        plt.close()

    def _plot_group_all(prefix: str, label: str):
        group_keys = sorted([k for k in features.keys() if k.startswith(prefix)])
        if not group_keys:
            return
        pp = int(panels_per_page)
        if pp <= 0:
            pp = len(group_keys)
        pages = int(np.ceil(len(group_keys) / pp))
        for p in range(pages):
            start = p * pp
            end = min((p + 1) * pp, len(group_keys))
            chunk = group_keys[start:end]
            suffix = f"_p{p + 1:02d}" if pages > 1 else ""
            _plot_tiles(
                chunk,
                out_path=f"{out_prefix}_{label}_quicklooks{suffix}.png",
                title=f"{label} quicklooks ({start + 1}–{end} / {len(group_keys)})",
            )

    # 1) RGB composite from first three MNF components (if present)
    mnf_keys = sorted([k for k in features.keys() if k.startswith("MNF_")])
    if all(k in features for k in ("MNF_01", "MNF_02", "MNF_03")):
        r = _percentile_scale(as_array(features["MNF_01"]))
        g = _percentile_scale(as_array(features["MNF_02"]))
        b = _percentile_scale(as_array(features["MNF_03"]))
        rgb = np.dstack([r, g, b])
        plt.figure(figsize=(6, 6))
        plt.imshow(rgb)
        plt.title("RGB composite: MNF_01, MNF_02, MNF_03")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_MNF_RGB.png", dpi=200)
        plt.close()
    elif len(mnf_keys) >= 3:
        # Fallback: use the first 3 MNF_* keys if naming differs.
        r = _percentile_scale(as_array(features[mnf_keys[0]]))
        g = _percentile_scale(as_array(features[mnf_keys[1]]))
        b = _percentile_scale(as_array(features[mnf_keys[2]]))
        rgb = np.dstack([r, g, b])
        plt.figure(figsize=(6, 6))
        plt.imshow(rgb)
        plt.title(f"RGB composite: {mnf_keys[0]}, {mnf_keys[1]}, {mnf_keys[2]}")
        plt.axis("off")
        plt.tight_layout()
        plt.savefig(f"{out_prefix}_MNF_RGB.png", dpi=200)
        plt.close()

    # 2) Grouped quicklooks for all layers of interest
    _plot_group_all("Depth_", "Depth")
    _plot_group_all("Deriv_", "Deriv")
    _plot_group_all("MNF_", "MNF")

    # 3) Small mixed overview (legacy behavior) for a quick glance
    keys = sorted(list(features.keys()))
    sel = keys[: int(max_panels) if max_panels is not None else len(keys)]
    if sel:
        _plot_tiles(sel, out_path=f"{out_prefix}_feature_quicklooks.png", title="Feature overview")


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
    # Conservative initial feature set:
    # - MNF_01 only
    # - Derivative PCs 1–2 only
    n_mnf=1,
    n_deriv_pcs=2,
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
print(f"Extracted {len(features)} feature layers")

# QA quicklooks of feature bands (saves PNGs alongside output_prefix)
try:
    output_prefix_abs = os.path.abspath(output_prefix)
    plot_feature_quicklooks(features, output_prefix_abs)
    print(f"Saved feature quicklooks to {output_prefix_abs}_*.png")
except Exception as e:
    print(f"[WARN] Could not create feature quicklooks: {e}")

# Column-median diagnostic on a SWIR band after detrending (saves PNG alongside output_prefix)
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
    print(f"Saved detrended SWIR column-median diagnostic to {diag_path}")
except Exception as e:
    print(f"[WARN] Could not create SWIR column-median diagnostic: {e}")


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
