import rasterio
from rasterio.enums import Resampling
import numpy as np
import os
from rasterio.windows import from_bounds
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from bayseg import BaySeg
from bayseg.bayseg import labels_map, compute_labels_prob, compute_ie
import time

def crop_to_bounds(src, bounds, transform):
    """
    Crop an image to specified geographic bounds
    bounds: tuple of (xmin, ymin, xmax, ymax) in the same CRS as the image
    """
    # Convert bounds to pixel coordinates
    row_start, col_start = ~transform * (bounds[0], bounds[3])
    row_stop, col_stop = ~transform * (bounds[2], bounds[1])
    
    # Convert to integers and ensure they're within bounds
    row_start, col_start = int(row_start), int(col_start)
    row_stop, col_stop = int(row_stop), int(col_stop)
    
    # Ensure we don't go out of bounds
    row_start = max(0, row_start)
    col_start = max(0, col_start)
    row_stop = min(src.height, row_stop)
    col_stop = min(src.width, col_stop)
    
    print(f"Crop window: rows {row_start}:{row_stop}, cols {col_start}:{col_stop}")
    
    # Read the data using the window
    data = src.read(1, window=((row_start, row_stop), (col_start, col_stop)))
    print(f"Cropped data shape: {data.shape}")
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

def main():
    # Choose bands best suited for lithology
    bands = {
        "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B04_60m.jp2",   # Red
        "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B06_60m.jp2",   # Red Edge 2
        "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B07_60m.jp2",   # Red Edge 3
        "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B8A_60m.jp2",  # Narrow NIR
        "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B11_60m.jp2",  # SWIR 1
        "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_B12_60m.jp2",  # SWIR 2
        "TCI": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis AOI 1/S2A_MSIL2A_20230829T110621_N0509_R137_T29SQB_20230829T152901.SAFE/GRANULE/L2A_T29SQB_A042748_20230829T111659/IMG_DATA/R60m/T29SQB_20230829T110621_TCI_60m.jp2"
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
        if name != "TCI":  # Skip TCI for the stack
            with rasterio.open(path) as src:
                # Simply read the band data
                band_data = src.read(1)
                band_data = crop_by_rectangle(band_data, 226, 516, 565, 935)
                stack.append(band_data)

    # Stack into (rows, cols, bands)
    img_stack = np.stack(stack, axis=-1).astype(np.float64)
    print(f"Image stack shape: {img_stack.shape}")

    # Save intermediate result
    np.save("sentinel2_bayseg_input.npy", img_stack)
    print("Intermediate data saved to sentinel2_bayseg_input.npy")

    # Plot TCI for reference
    with rasterio.open(bands["TCI"]) as src:
        if src.count == 3:  # If TCI has 3 bands
            tci_data = src.read()  # Read all bands
            tci_data = tci_data.transpose(1, 2, 0)  # Change to (rows, cols, bands)
            tci_data = crop_by_rectangle(tci_data, 226, 516, 565, 935)

            plt.figure(figsize=(10, 10))
            plt.imshow(tci_data)
            plt.title('True Color Image (TCI) of Cropped Region')
            plt.axis('off')
            plt.show()

    print("\nStep 2: Running segmentation...")
    # Initialize segmenter
    n_classes = 6
    print(f"Running segmentation with {n_classes} classes...")
    start_time = time.time()
    seg = BaySeg(data=img_stack, n_labels=n_classes, beta_init=10)

    # Fit the model and get the labels
    seg.fit(n=200, beta_jump_length=0.1)
    
    # Get the final labels (MAP estimate)
    final_labels = labels_map(seg.labels)
    final_labels = final_labels.reshape(img_stack.shape[0], img_stack.shape[1])
    
    # Calculate information entropy
    labels_prob = compute_labels_prob(np.array(seg.labels))
    entropy = compute_ie(labels_prob)
    entropy = entropy.reshape(img_stack.shape[0], img_stack.shape[1])
    
    print(f"Segmentation completed in {time.time() - start_time:.2f} seconds")

    # Save results as numpy arrays
    np.save(f"bayseg_lithology_labels_n{n_classes}.npy", final_labels)
    np.save(f"bayseg_entropy_n{n_classes}.npy", entropy)
    print(f"Results saved to bayseg_lithology_labels_n{n_classes}.npy and bayseg_entropy_n{n_classes}.npy")

    # Export as georeferenced GeoTIFF
    # Get the transform for the cropped region
    with rasterio.open(bands["B4"]) as src:
        # Calculate the transform for the cropped region
        # The crop window was: rows 226:516, cols 565:935
        transform = src.transform * rasterio.Affine.translation(565, 226)
        
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
        output_path = f"segmentation_result_n{n_classes}.tif"
        with rasterio.open(output_path, 'w', **profile) as dst:
            dst.write(final_labels.astype(np.uint8), 1)
        print(f"Georeferenced segmentation results saved to {output_path}")

        # Save the entropy as a GeoTIFF
        entropy_profile = profile.copy()
        entropy_profile.update({
            'dtype': 'float32',
            'nodata': -9999
        })
        entropy_path = f"segmentation_entropy_n{n_classes}.tif"
        with rasterio.open(entropy_path, 'w', **entropy_profile) as dst:
            dst.write(entropy.astype(np.float32), 1)
        print(f"Georeferenced entropy results saved to {entropy_path}")

    # Plot diagnostic information
    seg.diagnostics_plot(transpose=False)

    # Plot acceptance ratios
    seg.plot_acc_ratios()

if __name__ == "__main__":
    main() 