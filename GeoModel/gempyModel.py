import os
import gempy as gp
import pandas as pd
import gempy_viewer as gpv

from GeoModel.HelperMethods import clean_topo_file
from GeoModel.HelperMethods import create_cross_section

# -------------------------------
# CONFIGURATION
# -------------------------------
path_to_data = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints.csv"
path_to_orientations = r"C:\Users\maxha\OneDrive\Desktop\orientations.csv"
path_to_topography = r"C:\Users\maxha\OneDrive\Desktop\topo.tif"
path_to_topography_cleaned = r"C:\Users\maxha\OneDrive\Desktop\topo_cleaned.tif"

invalid_below = -100  # Define what we consider invalid topography
points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')

# min max values from the data
min_x = points_df['X'].min()
max_x = points_df['X'].max()
min_y = points_df['Y'].min()
max_y = points_df['Y'].max()
max_z = points_df['Z'].max()

# -------------------------------
# CLEAN TOPOGRAPHY IF NEEDED
# -------------------------------
if not os.path.exists(path_to_topography_cleaned): clean_topo_file(path_to_topography, path_to_topography_cleaned)
else: print("Cleaned topography file already exists — skipping cleaning.")

# -------------------------------
# MODEL CREATION
# -------------------------------

geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, min_y, max_y, -1000, max_z],
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

structural_frame:gp.data.StructuralFrame = geo_model.structural_frame
structural_frame.structural_groups = structural_frame.structural_groups[0:2]
gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=20, legend=False, show_data=True)

gempy_model = gp.compute_model(geo_model)

# -------------------------------
# PLOTTING
# -------------------------------

#gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_cleaned)
gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=20, legend=False, show_data=True)
#create_cross_section(geo_model, cross_section=10)

pass
# -------------------------------
# MISCELLANEOUS
# -------------------------------

"""
gpv.plot_2d(geo_model, show_value=True, show_lith=False, show_scalar=True, legend=False)
gpv.plot_2d(geo_model, show_boundaries=False, show_data=False, direction="z", legend=False)
gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography_cleaned)
gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=20, legend=False, show_data=True)


After creating the model, modify colors of structural elements
for element in geo_model.structural_frame.structural_elements:
    if element.name == "Visean Conglomerates":
        element.color = "#ff4c00"
    elif element.name == "Visean Mudstones":
        element.color = "#b200b2"
    elif element.name == "Famennian Siliciclastics":
        element.color = "#d9b280"
    elif element.name == "Frasnian Siliciclastics":
        element.color = "#ffedb2"
"""