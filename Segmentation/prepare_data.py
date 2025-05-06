import rasterio
from rasterio.enums import Resampling
import numpy as np
import os
from rasterio.windows import from_bounds
import matplotlib.pyplot as plt

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
    with rasterio.open(path) as src:
        # Simply read the band data
        band_data = src.read(1)
        band_data = crop_by_rectangle(band_data,  226, 516,565, 935)
        stack.append(band_data)

# Stack into (rows, cols, bands)
img_stack = np.stack(stack, axis=-1).astype(np.float64)

# Save for BaySeg
print(f"Image stack shape: {img_stack.shape}")
np.save("sentinel2_bayseg_input.npy", img_stack)



with rasterio.open(bands["TCI"]) as src:
    # For TCI, we need to read all bands (RGB)
    if src.count == 3:  # If TCI has 3 bands
        tci_data = src.read()  # Read all bands
        tci_data = tci_data.transpose(1, 2, 0)  # Change to (rows, cols, bands)
        tci_data = crop_by_rectangle(tci_data,  226, 516, 565, 935)

        # Plot TCI properly
        plt.figure(figsize=(10, 10))
        plt.imshow(tci_data)
        plt.title('True Color Image (TCI) of Cropped Region')
        plt.axis('off')
        plt.show()

