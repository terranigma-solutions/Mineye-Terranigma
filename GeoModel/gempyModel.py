import os
import gempy as gp
import pandas as pd
import gempy_viewer as gpv
import GeoModel.HelperMethods as helper

# -------------------------------
# CONFIGURATION
# -------------------------------
path_to_data = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints_reduced.csv"
path_to_orientations = r"C:\Users\maxha\OneDrive\Desktop\orientations_reduced.csv"
path_to_topography = r"C:\Users\maxha\OneDrive\Desktop\topo.tif"
path_to_topography_cleaned = r"C:\Users\maxha\OneDrive\Desktop\topo_cleaned.tif"
path_to_topography_reduced = r"C:\Users\maxha\OneDrive\Desktop\topo_reduced.tif"

invalid_below = -100  # Define what we consider invalid topography
model_depth = -500  # Define the depth of the model
topo_reduction_factor = 0.1  # Define the factor by which to reduce the topography resolution

trigger_create_cross_section = False # Set to True if you want to create a cross section
trigger_drop_lithologies = False  # Set to True if you want to drop certain lithologies

lithologies_to_drop = [
    "Upper Carboniferous Volcanics",
    "Tournaisian Plutonites"
]

points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')
orientations_df = pd.read_csv(path_to_orientations, encoding='latin1', engine='python')

# min max values from the data
min_x = points_df['X'].min()
max_x = points_df['X'].max()
min_y = points_df['Y'].min()
max_y = points_df['Y'].max()
max_z = points_df['Z'].max()

# -------------------------------
# CLEAN AND PREPROCESS DATA
# -------------------------------
helper.clean_topo_file(path_to_topography, path_to_topography_cleaned)
helper.reduce_tif_resolution(
    input_path=path_to_topography_cleaned,
    output_path=path_to_topography_reduced,
    scale_factor=topo_reduction_factor  # Adjust scale factor as needed
   )

base, ext = os.path.splitext(path_to_topography_reduced)
path_to_topography_reduced = f"{base}_sf{topo_reduction_factor}{ext}"

if trigger_drop_lithologies:
    helper.drop_lithologies(points_df, orientations_df, lithologies_to_drop)

# -------------------------------
# MODEL CREATION
# -------------------------------

geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, max_y, min_y, model_depth, max_z],
    refinement=4,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations="temp_orientations.csv",
        path_to_surface_points="temp_points.csv",
    )
)

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        #"Strat_Series2": ("Visean Shales", "Mid Carboniferous Shales"),
        #"Strat_Series2": ("Upper Carboniferous Volcanics", "Visean Shales"),
        #"Strat_Series2": "Tournaisian Plutonites",
        "Strat_Series1": ("Upper Devonian Siliciclastics", "Visean Shales", "Mid Carboniferous Shales")
    }

)

helper.color_lithology(geo_model.structural_frame.structural_elements)
gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_reduced)
gempy_model = gp.compute_model(geo_model)

# -------------------------------
# PLOTTING
# -------------------------------

gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=10, legend=False, show_data=True, show_topography=True, transformed_data=False)
if trigger_create_cross_section:
    helper.create_cross_section(geo_model, cross_section=5)

pass
# -------------------------------
# MISCELLANEOUS
# -------------------------------

"""
gpv.plot_2d(geo_model, show_value=True, show_lith=False, show_scalar=True, legend=False)
gpv.plot_2d(geo_model, show_boundaries=False, show_data=False, direction="z", legend=False)
gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_cleaned)
gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=20, legend=False, show_data=True)

geo_model.structural_frame.structural_groups = geo_model.structural_frame.structural_groups[1:2]
geo_model.input_transform.scale[2] *= 2
"""