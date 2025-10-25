import os
from raster_to_obj import main

def test_raster_to_obj():
    # Replace these paths with your actual file paths
    raster_path = "/Users/simonvirgo/Downloads/spremberg-srtm-and-map/Maps/spremberg-google-sat.png"
    dem_path = "/Users/simonvirgo/Downloads/spremberg-srtm-and-map/SRTM/spremberg-srtm.tif"
    output_path = "output/spremberg/spremberg-google-sat.obj"
    base_height = 220.0 # height to set null values
    raster_epsg = 31469  # EPSG code for the raster
    dem_epsg = 31469
    scale_factor = 2#

    # Check if input files exist
    if not os.path.exists(raster_path):
        print(f"Error: Raster file not found at {raster_path}")
        return

    if not os.path.exists(dem_path):
        print(f"Error: DEM file not found at {dem_path}")
        return

    try:
        # Run the conversion
        main(raster_path, dem_path, output_path, base_height, raster_epsg, dem_epsg, scale_factor)
        print(f"Successfully created OBJ file at {output_path}")
    except Exception as e:
        print(f"Error during conversion: {str(e)}")

if __name__ == "__main__":
    test_raster_to_obj()