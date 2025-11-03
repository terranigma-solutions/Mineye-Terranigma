import sys
import os
# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import pandas as pd
import gempy as gp
import gempy_viewer as gpv

from mineye.GeoModel import HelperMethods, Plotter
from mineye.config import paths


# ========== CONFIG ==========
BASE_DIR = paths.get_base_dir()
geomodel_dir = paths.get_geomodel_dir(BASE_DIR)
tmp_dir = paths.get_tmp_dir(BASE_DIR)

gis_map_info = paths.get_gis_map_info_path(BASE_DIR)
topography_path = paths.get_topography_path(BASE_DIR)
pickle_model_path = paths.get_pickle_model_path(BASE_DIR)
mod_or_path = paths.get_orientations_path(BASE_DIR)
mod_pts_path = paths.get_points_path(BASE_DIR)

# Load original full points once (for optional comparison)
original_points = pd.read_csv(paths.get_contact_points_path(BASE_DIR))
points_df = original_points.copy()


min_x = -709521
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
simplification_level = 0.9  # 0=no simplification, 1=maximum simplification
dip_values = 10
azimuth_values = 0

# Set the extent of the model
extent=[min_x, max_x, min_y, max_y, model_depth, max_z]  # Fixed Y coordinate order

# Set model resolution
nx = ny = nz = 64
resolution = [nx, ny, nz]

# ========== DATA PROCESSING ==========
orientations_df, points_df = HelperMethods.process_geological_data(
    points_df=points_df,
    default_dip=dip_values,
    default_azimuth=azimuth_values,
    use_default_azimuth=False,
    boundary_tolerance=800,
    formation_id=34,
    simplification_level=simplification_level,
    manual_orientations_at_points=[11, 19],
    azimuth_flip_by_id={
        0: True, 1: True, 2: True, 3: True, 4: True, 5: False, 6: False
    },
    manual_dip_by_id={}  # Add manual dip values as needed: {3: 60, 7: 30}
)

# ========== SAVE PROCESSED DATA ==========
orientations_df.to_csv(mod_or_path, index=False)
points_df.to_csv(mod_pts_path, index=False)


# ========== BUILD GEMPY MODEL ==========
simple_geo_model = gp.create_geomodel(
    project_name='simple_model',
    extent=extent,
    refinement=4,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }
)

gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topography_path)
geo_model = gp.compute_model(simple_geo_model)

p = gpv.plot_2d(
    simple_geo_model,
    section_names=['topography'],   # this triggers the top-down geological map
    show_topography=True,
    show_lith=True,
    show_boundaries=True,
    show_data=True,
    legend=False
)

Plotter.create_cross_section(simple_geo_model, cross_section=5)
gpv.plot_3d(simple_geo_model, section_names=['topography'], show_topography=True, ve=12)