import sys
import os
# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import gempy as gp
import gempy_viewer as gpv

from mineye.GeoModel import Plotter
from mineye.config import paths


# ========== CONFIG ==========
BASE_DIR = paths.get_base_dir()
topography_path = paths.get_topography_path(BASE_DIR)

# Load pre-processed orientations and points from temp_inputs
mod_or_path = paths.get_orientations_path(BASE_DIR)
mod_pts_path = paths.get_points_path(BASE_DIR)

# Model extent
min_x = -709521
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
extent = [min_x, max_x, min_y, max_y, model_depth, max_z]

# Model resolution
nx = ny = nz = 64
resolution = [nx, ny, nz]



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
gp.compute_model(simple_geo_model)

# ========== VISUALIZATION ==========
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

# Note: 3D plotting currently has compatibility issues with gempy_viewer
# gpv.plot_3d(simple_geo_model, section_names=['topography'], show_topography=True, ve=12)
