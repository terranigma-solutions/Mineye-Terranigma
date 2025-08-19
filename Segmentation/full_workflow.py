
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

def crop_to_bounds(src, bounds, transform):
    """
    Crop an image to specified geographic bounds.

    Args:
        src: An open rasterio dataset (used to read and to access size/bounds).
        bounds: Tuple (xmin, ymin, xmax, ymax) in the same CRS as the image.
        transform: Affine transform of the raster (usually src.transform).

    Returns:
        data: np.ndarray of shape (bands, rows, cols) for the cropped area.
        out_transform: Affine transform for the cropped data.
        window: rasterio.windows.Window used for the read.

    Raises:
        ValueError: If there is no overlap between the requested bounds and the raster.
    """
    # Local imports to avoid changing module-level imports
    from math import floor, ceil
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
        raise ValueError("Requested bounds do not overlap the raster.")

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

    # Read the data using the window
    data = src.read(window=window)

    # Compute the new transform for the cropped data
    out_transform = win_transform(window, transform)

    # Optional: helpful debug print
    print(f"Crop window -> rows {int(window.row_off)}:{int(window.row_off + window.height)}, "
          f"cols {int(window.col_off)}:{int(window.col_off + window.width)}")

    return data, out_transform, window

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

def apply_soil_mask(img_stack, scl_path, crop_rectangle):
    """
    Apply soil mask to the image stack using SCL layer
    Args:
        img_stack: numpy array of shape (rows, cols, bands)
        scl_path: path to the SCL layer
        crop_rectangle: tuple of (row_start, row_end, col_start, col_end)
    Returns:
        masked image stack with -1 for non-soil pixels
    """
    # Load and process the SCL layer to create a soil mask
    with rasterio.open(scl_path) as src:
        scl_data = src.read(1)
        scl_data = crop_by_rectangle(scl_data, *crop_rectangle)
        # Create soil mask (SCL class 5 represents bare soil)
        soil_mask = (scl_data == 5)

    # Apply soil mask to the image stack
    soil_mask_3d = np.repeat(soil_mask[:, :, np.newaxis], img_stack.shape[2], axis=2)
    masked_stack = np.where(soil_mask_3d, img_stack, -9999)
    
    print(f"Number of soil pixels: {np.sum(soil_mask)}")
    print(f"Percentage of soil pixels: {np.sum(soil_mask) / soil_mask.size * 100:.2f}%")
    
    return masked_stack

def main():
    # Define the crop rectangle
    crop_rectangle = (246, 500, 554, 954)  # (row_start, row_end, col_start, col_end) the width and height needs to be even

    # Segmentation parameters
    n_classes = 4
    beta_init = 30
    n_iterations = 400
    beta_jump_length = 0.1

    # Choose bands best suited for lithology
    bands = {
        "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B04_60m.jp2",   # Red
        "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B06_60m.jp2",   # Red Edge 2
        "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B07_60m.jp2",   # Red Edge 3
        "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B8A_60m.jp2",  # Narrow NIR
        "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B11_60m.jp2",  # SWIR 1
        "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B12_60m.jp2",  # SWIR 2
        # Not used for segmentation:
        "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_TCI_60m.jp2",
        "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_SCL_60m.jp2"  # Scene Classification Layer
    }

    # bounds = (xmin, ymin, xmax, ymax)
    bounds = (733891.6, 4168988.9, 756038.65, 4186614.52)  # in the same CRS as the image

    print("Step 1: Preparing data...")
    with rasterio.open(bands["B4"]) as ref:
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

    stack = []
    for name, path in bands.items():
        if name not in ["TCI", "SCL"]:  # Skip TCI and SCL for the stack
            with rasterio.open(path) as src:
                # Simply read the band data
                band_data = src.read(1)
                band_data = crop_by_rectangle(band_data, *crop_rectangle)
                stack.append(band_data)

    # Stack into (rows, cols, bands)
    img_stack = np.stack(stack, axis=-1).astype(np.float64)
    print(f"Image stack shape: {img_stack.shape}")

    # Apply soil mask
    img_stack = apply_soil_mask(img_stack, bands["SCL"], crop_rectangle)

    # Save intermediate result
    np.save("sentinel2_bayseg_input.npy", img_stack)
    print("Intermediate data saved to sentinel2_bayseg_input.npy")

    # Plot TCI for reference
    with rasterio.open(bands["TCI"]) as src:
        if src.count == 3:  # If TCI has 3 bands
            tci_data = src.read()  # Read all bands
            tci_data = tci_data.transpose(1, 2, 0)  # Change to (rows, cols, bands)
            tci_data = crop_by_rectangle(tci_data, *crop_rectangle)

            plt.figure(figsize=(10, 10))
            plt.imshow(tci_data)
            plt.title('True Color Image (TCI) of Cropped Region')
            plt.axis('off')
            plt.show()

    print("\nStep 2: Running segmentation...")
    # Initialize segmenter
    print(f"Running segmentation with {n_classes} classes...")
    start_time = time.time()
    seg = BaySeg(data=img_stack, n_labels=n_classes, beta_init=beta_init)

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
    np.save(f"bayseg_lithology_labels_n{n_classes}_beta{beta_jump_length}.npy", final_labels)
    np.save(f"bayseg_entropy_n{n_classes}_beta{beta_jump_length}.npy", entropy)
    print(f"Results saved to bayseg_lithology_labels_n{n_classes}_beta{beta_jump_length}.npy and bayseg_entropy_n{n_classes}_beta{beta_jump_length}.npy")

    # Export as georeferenced GeoTIFF
    # Get the transform for the cropped region
    with rasterio.open(bands["B4"]) as src:
        # Calculate the transform for the cropped region using the actual crop rectangle
        transform = src.transform * rasterio.Affine.translation(crop_rectangle[2], crop_rectangle[0])
        
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
        output_path = f"segmentation_result_n{n_classes}_beta{beta_jump_length}.tif"
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(final_labels.astype(np.uint8), 1)
        print(f"Georeferenced segmentation results saved to {output_path}")

        # Save the entropy as a GeoTIFF
        entropy_profile = profile.copy()
        entropy_profile.update({
            'dtype': 'float32',
            'nodata': -9999
        })
        entropy_path = f"segmentation_entropy_n{n_classes}_beta{beta_jump_length}.tif"
        with rasterio.open(entropy_path, 'w', **entropy_profile) as dst:
            dst.write(entropy.astype(np.float32), 1)
        print(f"Georeferenced entropy results saved to {entropy_path}")

    # Plot diagnostic information
    seg.diagnostics_plot(transpose=False)

    # Plot acceptance ratios
    seg.plot_acc_ratios()

if __name__ == "__main__":
    main()