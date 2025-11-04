import sys
import os
import gempy as gp
import gempy_viewer as gpv

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from mineye.GeoModel import Plotter
from mineye.config.paths import (
    get_base_dir,
    get_topography_path,
    get_orientations_path,
    get_points_path,
)
from mineye.config.example_parameters import TharsisModelConfig

# ========== CONFIG ==========
BASE_DIR = get_base_dir()
topography_path = get_topography_path(BASE_DIR)

# Load pre-processed orientations and points from temp_inputs
mod_or_path = get_orientations_path(BASE_DIR)
mod_pts_path = get_points_path(BASE_DIR)

# ========== BUILD GEMPY MODEL ==========
simple_geo_model = gp.create_geomodel(
    project_name=TharsisModelConfig.PROJECT_NAME,
    extent=TharsisModelConfig.EXTENT,
    refinement=TharsisModelConfig.REFINEMENT,
    resolution=TharsisModelConfig.RESOLUTION,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object=TharsisModelConfig.SURFACE_MAPPING
)

gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topography_path)
gp.compute_model(simple_geo_model)

# ========== VISUALIZATION ==========
p = gpv.plot_2d(
    simple_geo_model,
    section_names=['topography'],   # this triggers the top-down geological map
    show_topography=TharsisModelConfig.SHOW_TOPOGRAPHY,
    show_lith=TharsisModelConfig.SHOW_LITH,
    show_boundaries=TharsisModelConfig.SHOW_BOUNDARIES,
    show_data=TharsisModelConfig.SHOW_DATA,
    legend=TharsisModelConfig.SHOW_LEGEND
)

Plotter.create_cross_section(
    simple_geo_model,
    cross_section=TharsisModelConfig.VIZ_CROSS_SECTION_COUNT,
    vertical_exaggeration=TharsisModelConfig.VIZ_VERTICAL_EXAGGERATION
)

gpv.plot_3d(
    simple_geo_model,
    section_names=['topography'],
    show_topography=TharsisModelConfig.SHOW_TOPOGRAPHY,
    ve=TharsisModelConfig.VIZ_3D_VERTICAL_EXAGGERATION
)
