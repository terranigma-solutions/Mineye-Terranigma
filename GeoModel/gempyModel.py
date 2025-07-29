import os
import gempy as gp
import pandas as pd
import gempy_viewer as gpv

from GeoModel.HelperMethods import clean_topo_file
from GeoModel.HelperMethods import create_cross_section
from GeoModel.HelperMethods import reduce_tif_resolution

# -------------------------------
# CONFIGURATION
# -------------------------------
path_to_data = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints.csv"
path_to_orientations = r"C:\Users\maxha\OneDrive\Desktop\orientations.csv"
path_to_topography = r"C:\Users\maxha\OneDrive\Desktop\topo.tif"
path_to_topography_cleaned = r"C:\Users\maxha\OneDrive\Desktop\topo_cleaned.tif"
path_to_topography_reduced = r"C:\Users\maxha\OneDrive\Desktop\topo_reduced.tif"

invalid_below = -100  # Define what we consider invalid topography
model_depth = -500  # Define the depth of the model
topo_reduction_factor = 0.1  # Define the factor by which to reduce the topography resolution
trigger_create_cross_section = False # Set to True if you want to create a cross section

points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')

# min max values from the data
min_x = points_df['X'].min()
max_x = points_df['X'].max()
min_y = points_df['Y'].min()
max_y = points_df['Y'].max()
max_z = points_df['Z'].max()

# -------------------------------
# CLEAN AND SCALE TOPOGRAPHY IF NEEDED
# -------------------------------
clean_topo_file(path_to_topography, path_to_topography_cleaned)
reduce_tif_resolution(
    input_path=path_to_topography_cleaned,
    output_path=path_to_topography_reduced,
    scale_factor=topo_reduction_factor  # Adjust scale factor as needed
)

base, ext = os.path.splitext(path_to_topography_reduced)
path_to_topography_reduced = f"{base}_sf{topo_reduction_factor}{ext}"

# -------------------------------
# MODEL CREATION
# -------------------------------

geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, max_y, min_y, model_depth, max_z],
    refinement=4,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=path_to_orientations,
        path_to_surface_points=path_to_data
    )
)

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Strat_Series4": "Mid Carboniferous Shales",
        "Strat_Series3": ("Upper Carboniferous Volcanics", "Visean Shales"),
        "Strat_Series2": "Tournaisian Plutonites",
        "Strat_Series1": "Upper Devonian Siliciclastics",
    }
)

for element in geo_model.structural_frame.structural_elements:
    if element.name == "Mid Carboniferous Shales":
        element.color = "#b2d9b2" #done - light green
    elif element.name == "Upper Carboniferous Volcanics":
        element.color = "#ff8000" #done - orange
    elif element.name == "Visean Shales":
        element.color = "#b200b2" #done - purple
    elif element.name == "Tournaisian Plutonites":
        element.color = "#e37ecd" #done - pink
    elif element.name == "Upper Devonian Siliciclastics":
        element.color = "#d9b280" #done - light brown

gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_reduced)
gempy_model = gp.compute_model(geo_model)

# -------------------------------
# PLOTTING
# -------------------------------

gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=10, legend=False, show_data=True)
if trigger_create_cross_section:
    create_cross_section(geo_model, cross_section=5)

pass
# -------------------------------
# MISCELLANEOUS
# -------------------------------

"""
gpv.plot_2d(geo_model, show_value=True, show_lith=False, show_scalar=True, legend=False)
gpv.plot_2d(geo_model, show_boundaries=False, show_data=False, direction="z", legend=False)
gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_cleaned)
gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=20, legend=False, show_data=True)
"""