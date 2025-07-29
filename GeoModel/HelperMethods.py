import os
import gempy_viewer as gpv
import numpy as np
import rasterio
from skimage.transform import resize
#import liquid_earth_sdk.api.le_api as le_api

def clean_topo_file(input_path: str, output_path: str, invalid_below: float = -100):
    # Load, clean, and save cleaned topography TIFF

    with rasterio.open(input_path) as src:
        profile = src.profile
        data = src.read(1)
        nodata = src.nodata

    data_cleaned = np.where(
        (data == nodata) | (data <= invalid_below) | np.isnan(data),1,data)

    profile.update(nodata=1)
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(data_cleaned, 1)

def reduce_tif_resolution(input_path, output_path, scale_factor):
    # Add scale factor to output file name
    base, ext = os.path.splitext(output_path)
    output_path_sf = f"{base}_sf{scale_factor}{ext}"

    with rasterio.open(input_path) as src:
        data = src.read(1)
        profile = src.profile

        # Calculate new shape
        new_shape = (int(data.shape[0] * scale_factor), int(data.shape[1] * scale_factor))
        # Resample data
        data_resampled = resize(data, new_shape, preserve_range=True, anti_aliasing=True).astype(data.dtype)

        # Update profile
        profile.update(
            height=new_shape[0],
            width=new_shape[1],
            transform=src.transform * src.transform.scale(
                (src.width / new_shape[1]),
                (src.height / new_shape[0])
            )
        )

        with rasterio.open(output_path_sf, "w", **profile) as dst:
            dst.write(data_resampled, 1)


def push_model_to_local_server(gempy_model, le_api):
    unstructured_data = gempy_model.raw_arrays.meshes_to_subsurface()
    le_api.upload_mesh_to_new_space(
        space_name="[TEMP] Test python api",
        data=unstructured_data,
        file_name="ROI1",
        token=os.environ.get("APIKEY")
    )

def create_cross_section(geo_model, cross_section: int):
    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=6,
                        cell_number=i,
                        show_topography=True,
                        legend=True,
                        show_data=False,
                        direction="y")
        # Increment i to avoid infinite loop
        i += 1

    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=6,
                        cell_number=i,
                        show_topography=True,
                        legend=True,
                        show_data=False,
                        direction="x")
        # Increment i to avoid infinite loop
        i += 1