import gempy as gp
import pandas as pd
import gempy_viewer as gpv

# Load the data
path_to_data = r"C:\Users\maxha\OneDrive\Desktop\formationinputpoints.csv"
path_to_orientations = r"C:\Users\maxha\OneDrive\Desktop\orientations.csv"

# First, let's look at what formations exist in your data
points_df = pd.read_csv(path_to_data, encoding='latin1', engine='python')
idMapping = {"40": "Famennian Siliciclastics", "29": "Frasnian Siliciclastics", "67": "Visean Mudstones", "65": "Visean Conglomerates"}



# min max values from the data
min_x = points_df['X'].min()
max_x = points_df['X'].max()
min_y = points_df['Y'].min()
max_y = points_df['Y'].max()
max_z = points_df['Z'].max()
min_z = points_df['Z'].min()

geo_model = gp.create_geomodel(
    project_name='AOI',
    extent=[min_x, max_x, min_y, max_y, -2000, max_z],
    refinement=4,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=path_to_orientations,
        path_to_surface_points=path_to_data
    )
)

print(geo_model.structural_frame)

gpv.plot_2d(geo_model)

print("\nStructural elements in the model:")
for element in geo_model.structural_frame.structural_elements:
    print(f"Name={element.name}, ID={element.id}")

gp.map_stack_to_surfaces(
    gempy_model=geo_model,
    mapping_object={
        "Strat_Series3": ("Visean Mudstones", "Visean Conglomerates"),
        "Strat_Series2": "Famennian Siliciclastics",
        "Strat_Series1": "Frasnian Siliciclastics"
    }
)

structural_frame:gp.data.StructuralFrame = geo_model.structural_frame
pass