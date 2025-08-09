import os
import threading
import gempy as gp
import pandas as pd
import HelperMethods

# -------------------------------
# CONFIGURATION
# -------------------------------
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(BASE_DIR, '..', 'Data', 'Tharsis AOI 1')

path_to_data = os.path.join(data_dir, 'formationinputpoints_reduced.csv')
path_to_orientations = os.path.join(data_dir, 'orientations_reduced.csv')
path_to_topography = os.path.join(data_dir, 'topo.tif')
path_to_topography_cleaned = os.path.join(data_dir, 'topo_cleaned.tif')
path_to_topography_reduced = os.path.join(data_dir, 'topo_reduced.tif')

invalid_below = -100  # Define what we consider invalid topography
model_depth = -500  # Define the depth of the model
topo_reduction_factor = 0.1  # Define the factor by which to reduce the topography resolution

trigger_create_cross_section = True # Set to True if you want to create a cross section
trigger_drop_lithologies = False  # Set to True if you want to drop certain lithologies

lithologies_to_drop = [
    "Upper Carboniferous Volcanics",
    "Tournaisian Plutonites"
]

points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')
orientations_df = pd.read_csv(path_to_orientations, encoding='latin1', engine='python')

# min max values from the data
top_right = (-678192, 4549947)
bottom_right = (-677518, 4529455)
bottom_left = (-706792, 4528453)
top_left = (-707522, 4548927)

# Define the bounding box of the area of interest
min_x = min(top_left[0], bottom_left[0])
max_x = max(top_right[0], bottom_right[0])
min_y = min(bottom_left[1], bottom_right[1])
max_y = max(top_left[1], top_right[1])
max_z = points_df['Z'].max()


# -------------------------------
# CLEAN AND PREPROCESS DATA
# -------------------------------
HelperMethods.clean_topo_file(path_to_topography, path_to_topography_cleaned)
HelperMethods.reduce_tif_resolution(
    input_path=path_to_topography_cleaned,
    output_path=path_to_topography_reduced,
    scale_factor=topo_reduction_factor  # Adjust scale factor as needed
   )

base, ext = os.path.splitext(path_to_topography_reduced)
path_to_topography_reduced = f"{base}_sf{topo_reduction_factor}{ext}"

if trigger_drop_lithologies:
    HelperMethods.drop_lithologies(points_df, orientations_df, lithologies_to_drop)

# -------------------------------
# MODEL CREATION
# -------------------------------

geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, max_y, min_y, model_depth, max_z],
    refinement=6,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=os.path.join(BASE_DIR, "temp_orientations.csv"),
        path_to_surface_points=os.path.join(BASE_DIR, "temp_points.csv"),
    )
)

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        #"Strat_Series2": ("Visean Shales", "Mid Carboniferous Shales"),
        #"Strat_Series2": ("Upper Carboniferous Volcanics", "Visean Shales"),
        #"Strat_Series2": "Tournaisian Plutonites",
        "Strat_Series1": ( "Mid Carboniferous Shales", "Visean Shales", "Upper Devonian Siliciclastics")
    }
)

HelperMethods.color_lithology(geo_model.structural_frame.structural_elements)
gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_reduced)
gempy_model = gp.compute_model(geo_model)

# -------------------------------
# PLOTTING
# -------------------------------

threading.Thread(target=HelperMethods.plot_3d_async, args=(geo_model,)).start()
if trigger_create_cross_section:
    HelperMethods.create_cross_section(geo_model, cross_section=5)

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