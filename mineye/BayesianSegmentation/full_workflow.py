
import rasterio
from rasterio.enums import Resampling
import numpy as np
import os
from rasterio.windows import from_bounds
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from bayseg.bayseg import BaySeg
from bayseg.bayseg import labels_map, compute_labels_prob, compute_ie
import time
from contextlib import contextmanager

@contextmanager
def _as_rio_dataset(obj):
    """
    Yield a rasterio dataset from either a file path or an already-open rasterio dataset-like object.
    If a path is provided, it will be opened and closed within this context.
    If a dataset object is provided, it will be yielded as-is without closing.
    """
    if isinstance(obj, str):
        with rasterio.open(obj) as ds:
            yield ds
    elif hasattr(obj, 'read') and hasattr(obj, 'transform') and hasattr(obj, 'crs'):
        # If the dataset is closed, raise a clear error
        if hasattr(obj, 'closed') and getattr(obj, 'closed'):
            raise ValueError("Provided rasterio dataset object is closed. Please pass an open dataset.")
        yield obj
    else:
        raise TypeError("Band value must be a file path (str) or an open rasterio dataset object.")

def crop_to_bounds(src, bounds, transform):
    """
    Crop an image to specified geographic bounds.

    Args:
        src: An open rasterio dataset (used to read and to access size/bounds).
        bounds: Tuple (xmin, ymin, xmax, ymax) in the same CRS as the image.
        transform: Affine transform of the raster (usually src.transform).

    Returns:
        data: np.ndarray for the cropped area (same return type as crop_by_rectangle).

    Raises:
        ValueError: If there is no overlap between the requested bounds and the raster.
    """
    # Local imports to avoid changing module-level imports
    from rasterio.windows import from_bounds as win_from_bounds, transform as win_transform, Window

    # Unpack and normalize user-provided bounds
    xmin, ymin, xmax, ymax = bounds
    if xmin > xmax:
        xmin, xmax = xmax, xmin
    if ymin > ymax:
        ymin, ymax = ymax, ymin

    # Intersect requested bounds with raster bounds to avoid out-of-range windows
    r_left, r_bottom, r_right, r_top = src.bounds
    ix_left = max(xmin, r_left)
    ix_right = min(xmax, r_right)
    ix_bottom = max(ymin, r_bottom)
    ix_top = min(ymax, r_top)

    if not (ix_left < ix_right and ix_bottom < ix_top):
        r_bounds_str = f"raster bounds: left={r_left}, bottom={r_bottom}, right={r_right}, top={r_top}"
        req_bounds_str = f"requested bounds: left={xmin}, bottom={ymin}, right={xmax}, top={ymax}"
        raise ValueError(f"Requested bounds do not overlap the raster. ({req_bounds_str}; {r_bounds_str})")

    # Compute a pixel window from the intersected bounds
    window = win_from_bounds(ix_left, ix_bottom, ix_right, ix_top, transform=transform)

    # Round to pixel boundaries: floor the offsets and ceil the shape to cover the full area
    # This avoids losing edge pixels due to float rounding.
    window = window.round_offsets(op='floor').round_shape(op='ceil')

    # Ensure the window is within the raster dimensions
    raster_window = Window(col_off=0, row_off=0, width=src.width, height=src.height)
    window = window.intersection(raster_window)

    if window.width <= 0 or window.height <= 0:
        raise ValueError("Computed crop window is empty after clipping to raster bounds.")

    # Enforce even window dimensions by expanding by one pixel if needed
    # Adjust width (columns)
    if int(round(window.width)) % 2 == 1:
        # Prefer expanding to the right if possible
        if window.col_off + window.width < src.width:
            window = Window(col_off=window.col_off, row_off=window.row_off,
                            width=window.width + 1, height=window.height)
        # Otherwise, expand to the left if possible
        elif window.col_off > 0:
            window = Window(col_off=window.col_off - 1, row_off=window.row_off,
                            width=window.width + 1, height=window.height)
        # Clip to raster bounds just in case
        window = window.intersection(raster_window)

    # Adjust height (rows)
    if int(round(window.height)) % 2 == 1:
        # Prefer expanding downward if possible
        if window.row_off + window.height < src.height:
            window = Window(col_off=window.col_off, row_off=window.row_off,
                            width=window.width, height=window.height + 1)
        # Otherwise, expand upward if possible
        elif window.row_off > 0:
            window = Window(col_off=window.col_off, row_off=window.row_off - 1,
                            width=window.width, height=window.height + 1)
        # Clip to raster bounds just in case
        window = window.intersection(raster_window)

    # Read the data using the window; keep as numpy array output to match crop_by_rectangle
    data = src.read(window=window)

    # Optional: helpful debug print
    print(f"Crop window -> rows {int(window.row_off)}:{int(window.row_off + window.height)}, "
          f"cols {int(window.col_off)}:{int(window.col_off + window.width)}")

    # Return only the cropped data to match crop_by_rectangle's return type (numpy array)
    return data

def crop_by_rectangle(data, row_start, row_end, col_start, col_end):
    """
    Crop a numpy array using row and column indices
    Args:
        data: numpy array to crop
        row_start: starting row index
        row_end: ending row index
        col_start: starting column index
        col_end: ending column index
    Returns:
        cropped data array
    """
    # Ensure indices are within bounds
    row_start = max(0, row_start)
    col_start = max(0, col_start)
    row_end = min(data.shape[0], row_end)
    col_end = min(data.shape[1], col_end)
    
    print(f"Crop window: rows {row_start}:{row_end}, cols {col_start}:{col_end}")
    
    # Crop the data
    cropped_data = data[row_start:row_end, col_start:col_end]
    print(f"Cropped data shape: {cropped_data.shape}")
    return cropped_data

def apply_soil_mask(img_stack, scl_source, bounds):
    """
    Apply soil mask to the image stack using SCL layer
    Args:
        img_stack: numpy array of shape (rows, cols, bands)
        scl_source: path to the SCL layer or an open rasterio dataset
        bounds: tuple of (xmin, ymin, xmax, ymax)
    Returns:
        masked image stack with -1 for non-soil pixels
    """
    # Load and process the SCL layer to create a soil mask
    with _as_rio_dataset(scl_source) as src:
        scl_data = crop_to_bounds(src, bounds, src.transform)
        # Ensure 2D SCL array
        if scl_data.ndim == 3 and scl_data.shape[0] == 1:
            scl_data = scl_data[0]
        # Create soil mask (SCL class 5 represents bare soil)
        soil_mask = (scl_data == 5)

    # Apply soil mask to the image stack
    soil_mask_3d = np.repeat(soil_mask[:, :, np.newaxis], img_stack.shape[2], axis=2)
    masked_stack = np.where(soil_mask_3d, img_stack, -9999)
    
    print(f"Number of soil pixels: {np.sum(soil_mask)}")
    print(f"Percentage of soil pixels: {np.sum(soil_mask) / soil_mask.size * 100:.2f}%")
    
    return masked_stack

import argparse


def run_workflow(bands: dict,
                 bounds: tuple,
                 n_classes: int = 4,
                 beta_init: float = 30.0,
                 n_iterations: int = 400,
                 beta_jump_length: float = 0.1,
                 use_soil_mask: bool = True,
                 save_npy: bool = True,
                 plot_tci: bool = True,
                 ref_band: str = None,
                 output_prefix: str = "segmentation",
                 stencil: str | None = "4p"):
    """
    Generalized segmentation workflow.

    Args:
        bands: dict mapping band names to file paths or open rasterio datasets. May include optional keys 'SCL' and 'TCI'.
        bounds: (xmin, ymin, xmax, ymax) in the image CRS.
        n_classes: number of segmentation classes (labels).
        beta_init: initial beta for BaySeg.
        n_iterations: number of MCMC iterations.
        beta_jump_length: jump length for beta tuning.
        use_soil_mask: whether to apply SCL-based soil mask (needs bands['SCL']).
        save_npy: whether to save intermediate input and results as .npy files.
        plot_tci: whether to plot TCI if available in bands.
        ref_band: band name to use as georeferencing reference; if None, use 'B4' if present, else first band.
        output_prefix: prefix for output files (GeoTIFF and NPY outputs).
        stencil: neighborhood stencil for BaySeg ('4p' or '8p'). Default '4p' avoids NumPy ragged-array issues in older bayseg.
    """
    # Choose a reference band for metadata/geotransform
    band_keys = [k for k in bands.keys() if k not in ("TCI", "SCL")]
    if not band_keys:
        raise ValueError("No spectral bands provided in 'bands' (excluding TCI/SCL).")
    if ref_band is None:
        ref_band = "B4" if "B4" in bands else band_keys[0]
    if ref_band not in bands:
        raise ValueError(f"ref_band '{ref_band}' not found in provided bands.")

    print("Step 1: Preparing data...")
    with _as_rio_dataset(bands[ref_band]) as ref:
        profile = ref.profile
        height, width = ref.height, ref.width
        crs = ref.crs
        transform = ref.transform

        # Print CRS information
        print("\nCoordinate Reference System Information:")
        print(f"CRS: {crs}")
        print(f"Transform: {transform}")
        print(f"Image dimensions: {width}x{height} pixels")
        print(f"Resolution: {ref.res[0]}m x {ref.res[1]}m")

    # Build stack from provided bands (exclude TCI and SCL)
    stack = []
    band_names = sorted([k for k in bands.keys() if k not in ("TCI", "SCL")])
    first_shape = None
    for name in band_names:
        source = bands[name]
        with _as_rio_dataset(source) as src:
            if src.count != 1:
                raise ValueError(f"Band '{name}' must be single-band (count=1), got count={src.count}")
            band_data = crop_to_bounds(src, bounds, src.transform)
            # Ensure 2D per-band array
            if band_data.ndim == 3:
                if band_data.shape[0] == 1:
                    band_data = band_data[0]
                else:
                    raise ValueError(f"Band '{name}' returned {band_data.shape} after crop; expected single 2D layer")
            if first_shape is None:
                first_shape = band_data.shape
            elif band_data.shape != first_shape:
                raise ValueError(
                    f"Band '{name}' shape {band_data.shape} does not match first band shape {first_shape}. "
                    f"All bands must have identical size and coverage."
                )
            stack.append(band_data)

    # Stack into (rows, cols, bands) as required by BaySeg (Y,X,F)
    img_stack = np.stack(stack, axis=-1).astype(np.float64)
    print(f"Image stack shape: {img_stack.shape} (rows, cols, features). Using stencil={stencil}.")

    # Apply soil mask if requested and SCL available
    if use_soil_mask:
        scl_path = bands.get("SCL", None)
        if scl_path is None:
            print("[WARN] Soil mask requested but 'SCL' path not provided. Skipping soil mask.")
        else:
            img_stack = apply_soil_mask(img_stack, scl_path, bounds)

    # Save intermediate result
    if save_npy:
        np.save(f"{output_prefix}_input.npy", img_stack)
        print(f"Intermediate data saved to {output_prefix}_input.npy")

    # Plot TCI for reference if available and requested
    if plot_tci and ("TCI" in bands):
        with _as_rio_dataset(bands["TCI"]) as src:
            if src.count >= 3:  # If TCI has 3 bands
                tci_data = crop_to_bounds(src, bounds, src.transform)  # (bands, rows, cols)
                tci_data = tci_data.transpose(1, 2, 0)  # Change to (rows, cols, bands)

                plt.figure(figsize=(10, 10))
                plt.imshow(tci_data)
                plt.title('True Color Image (TCI) of Cropped Region')
                plt.axis('off')
                plt.show()
            else:
                print("[WARN] TCI provided but does not have 3 bands; skipping plot.")
    elif plot_tci:
        print("[INFO] No TCI provided; skipping TCI plot.")

    print("\nStep 2: Running segmentation...")
    # Initialize segmenter
    print(f"Running segmentation with {n_classes} classes...")
    start_time = time.time()
    seg = BaySeg(data=img_stack, n_labels=n_classes, beta_init=beta_init, stencil=stencil)

    # Fit the model and get the labels
    seg.fit(n=n_iterations, beta_jump_length=beta_jump_length)

    # Get the final labels (MAP estimate)
    final_labels = labels_map(seg.labels)
    final_labels = final_labels.reshape(img_stack.shape[0], img_stack.shape[1])

    # Calculate information entropy
    labels_prob = compute_labels_prob(np.array(seg.labels))
    entropy = compute_ie(labels_prob)
    entropy = entropy.reshape(img_stack.shape[0], img_stack.shape[1])

    print(f"Segmentation completed in {time.time() - start_time:.2f} seconds")

    # Save results as numpy arrays
    if save_npy:
        np.save(f"{output_prefix}_labels_n{n_classes}_betajump{beta_jump_length}.npy", final_labels)
        np.save(f"{output_prefix}_entropy_n{n_classes}_betajump{beta_jump_length}.npy", entropy)
        print(f"Results saved to {output_prefix}_labels_n{n_classes}_betajump{beta_jump_length}.npy and {output_prefix}_entropy_n{n_classes}_betajump{beta_jump_length}.npy")

    # Export as georeferenced GeoTIFF
    # Get the transform for the cropped region based on bounds using reference band
    with _as_rio_dataset(bands[ref_band]) as src:
        from rasterio.windows import from_bounds as win_from_bounds, Window, transform as win_transform
        # Recompute the same window used by crop_to_bounds
        window = win_from_bounds(bounds[0], bounds[1], bounds[2], bounds[3], transform=src.transform)
        window = window.round_offsets(op='floor').round_shape(op='ceil')
        raster_window = Window(col_off=0, row_off=0, width=src.width, height=src.height)
        window = window.intersection(raster_window)

        # Enforce even window dimensions to match crop_to_bounds behavior
        if int(round(window.width)) % 2 == 1:
            if window.col_off + window.width < src.width:
                window = Window(col_off=window.col_off, row_off=window.row_off,
                                width=window.width + 1, height=window.height)
            elif window.col_off > 0:
                window = Window(col_off=window.col_off - 1, row_off=window.row_off,
                                width=window.width + 1, height=window.height)
            window = window.intersection(raster_window)

        if int(round(window.height)) % 2 == 1:
            if window.row_off + window.height < src.height:
                window = Window(col_off=window.col_off, row_off=window.row_off,
                                width=window.width, height=window.height + 1)
            elif window.row_off > 0:
                window = Window(col_off=window.col_off, row_off=window.row_off - 1,
                                width=window.width, height=window.height + 1)
            window = window.intersection(raster_window)

        transform = win_transform(window, src.transform)

        # Create the output profile
        profile = src.profile.copy()
        profile.update({
            'driver': 'GTiff',
            'dtype': 'uint8',
            'count': 1,
            'height': final_labels.shape[0],
            'width': final_labels.shape[1],
            'transform': transform,
            'nodata': None
        })

        # Save the segmentation result as a GeoTIFF
        output_path = f"{output_prefix}_result_n{n_classes}_betajump{beta_jump_length}.tif"
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(final_labels.astype(np.uint8), 1)
        print(f"Georeferenced segmentation results saved to {output_path}")

        # Save the entropy as a GeoTIFF
        entropy_profile = profile.copy()
        entropy_profile.update({
            'dtype': 'float32',
            'nodata': -9999
        })
        entropy_path = f"{output_prefix}_entropy_n{n_classes}_betajump{beta_jump_length}.tif"
        with rasterio.open(entropy_path, 'w', **entropy_profile) as dst:
            dst.write(entropy.astype(np.float32), 1)
        print(f"Georeferenced entropy results saved to {entropy_path}")

    # Plot diagnostic information
    seg.diagnostics_plot(transpose=False)

    # Plot acceptance ratios
    seg.plot_acc_ratios()


def parse_band_arg(band_arg: str):
    """Parse a band argument of the form NAME=PATH."""
    if "=" not in band_arg:
        raise argparse.ArgumentTypeError("Band must be in NAME=PATH format, e.g., B4=/path/to/B04.jp2")
    name, path = band_arg.split("=", 1)
    name = name.strip()
    path = path.strip()
    if not name:
        raise argparse.ArgumentTypeError("Band name cannot be empty")
    if not os.path.exists(path):
        print(f"[WARN] Provided path does not exist: {path}")
    return name, path


def main():
    parser = argparse.ArgumentParser(description="Generalized segmentation workflow for Sentinel-2 or similar rasters.")
    parser.add_argument('--band', dest='bands', action='append', type=parse_band_arg, required=False,
                        help='Band specification as NAME=PATH. Repeat for multiple bands (exclude TCI/SCL here).')
    parser.add_argument('--scl', type=str, help='Path to SCL raster for soil masking (optional).')
    parser.add_argument('--tci', type=str, help='Path to TCI raster for plotting (optional).')
    parser.add_argument('--bounds', type=float, nargs=4, metavar=('XMIN','YMIN','XMAX','YMAX'),
                        help='Bounds in raster CRS: xmin ymin xmax ymax')
    parser.add_argument('--n-classes', type=int, default=4)
    parser.add_argument('--beta-init', type=float, default=30.0)
    parser.add_argument('--iterations', type=int, default=400)
    parser.add_argument('--beta-jump', type=float, default=0.1)
    parser.add_argument('--ref-band', type=str, default=None, help='Band name to use as reference for georeferencing.')
    parser.add_argument('--no-soil-mask', action='store_true', help='Disable soil masking even if SCL provided.')
    parser.add_argument('--no-plot-tci', action='store_true', help='Disable TCI plotting.')
    parser.add_argument('--no-save-npy', action='store_true', help='Do not save intermediate/results NPY files.')
    parser.add_argument('--output-prefix', type=str, default='segmentation', help='Prefix for output files.')

    # Backward-compatible defaults: if no args provided, use previous hardcoded setup
    args = parser.parse_args()

    # Build bands dict
    bands = {}
    if args.bands:
        for name, path in args.bands:
            bands[name] = path
    # Inject optional SCL/TCI if provided
    if args.scl:
        bands['SCL'] = args.scl
    if args.tci:
        bands['TCI'] = args.tci

    if not bands:
        # Fallback to previous hardcoded defaults to preserve old behavior
        bands = {
            "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B04_60m.jp2",
            "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B06_60m.jp2",
            "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B07_60m.jp2",
            "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B8A_60m.jp2",
            "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B11_60m.jp2",
            "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B12_60m.jp2",
            "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_TCI_60m.jp2",
            "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_SCL_60m.jp2",
        }

    # Bounds
    if args.bounds:
        bounds = tuple(args.bounds)
    else:
        bounds = (733891.6, 4168988.9, 756038.65, 4186614.52)

    run_workflow(
        bands=bands,
        bounds=bounds,
        n_classes=args.n_classes,
        beta_init=args.beta_init,
        n_iterations=args.iterations,
        beta_jump_length=args.beta_jump,
        use_soil_mask=(not args.no_soil_mask),
        save_npy=(not args.no_save_npy),
        plot_tci=(not args.no_plot_tci),
        ref_band=args.ref_band,
        output_prefix=args.output_prefix,
    )


if __name__ == "__main__":
    main()

