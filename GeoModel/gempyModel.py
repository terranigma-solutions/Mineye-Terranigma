import gempy as gp
import pandas as pd
import gempy_viewer as gpv

# Load the data
path_to_data = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints.csv"
path_to_orientations = r"C:\Users\maxha\OneDrive\Desktop\orientations.csv"
path_to_topography = r"C:\Users\maxha\OneDrive\Desktop\topo.tif"

# First, let's look at what formations exist in your data
points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')

# min max values from the data
min_x = points_df['X'].min()
max_x = points_df['X'].max()
min_y = points_df['Y'].min()
max_y = points_df['Y'].max()
max_z = points_df['Z'].max()

# First, print the unique formation names in both files
points_df = pd.read_csv(path_to_data)


geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, min_y, max_y, -1000, max_z],
    refinement=4,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=path_to_orientations,
        path_to_surface_points=path_to_data
    )
)

geologicalStack = [29, 40, 1, 51, 67, 76]
idToNameDict = {
    29: "Mid Devonian Siliciclastics",
    40: "Upper Devonian Siliciclastics",
    1: "Tournaisian Plutonites",
    51: "Upper Carboniferous Volcanics",
    67: "Visean Shales",
    76: "Mid Carboniferous Shales"
}

"""
# After creating the model, modify colors of structural elements
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

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Strat_Series4": "Mid Carboniferous Shales",
        "Strat_Series3": ("Upper Carboniferous Volcanics", "Visean Shales"),
        "Strat_Series2": "Tournaisian Plutonites",
        "Strat_Series1": "Upper Devonian Siliciclastics",
    }
)

gpv.plot_2d(geo_model)

#print(geo_model.structural_frame)
structural_frame:gp.data.StructuralFrame = geo_model.structural_frame

gp.compute_model(geo_model)
#gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=None, legend=False, show_data=False)

gp.set_topography_from_file(grid=geo_model.grid, filepath=path_to_topography)

i = 0
while i < 15:
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
while i < 15:
    if i % 2 == 0:
        gpv.plot_2d(geo_model, ve=6,
                   cell_number=i,
                   show_topography=True,
                   legend=True,
                   show_data=False,
                   direction="x")
    # Increment i to avoid infinite loop
    i += 1

#gpv.plot_2d(geo_model, show_value=True, show_lith=False, show_scalar=True, cell_number=15, legend=False)
#gpv.plot_2d(geo_model, show_boundaries=False, show_data=False, direction="z", legend=False)
pass