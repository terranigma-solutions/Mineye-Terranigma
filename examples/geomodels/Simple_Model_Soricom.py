import sys
import os
# Add project root to path for imports
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
import gempy as gp
import gempy_viewer as gpv

from mineye.GeoModel import Plotter
from mineye.config.paths import (
    get_base_dir,
    get_soricom_orientations_with_fault_path,
    get_soricom_formation_points_with_fault_path,
)
from mineye.config.example_parameters import SoricomModelConfig

# ========== CONFIG ==========
BASE_DIR = get_base_dir()

# Load Soricom data paths (with fault data already included)
orientation_path = get_soricom_orientations_with_fault_path(BASE_DIR)
formation_points_path = get_soricom_formation_points_with_fault_path(BASE_DIR)

# ========== BUILD GEMPY MODEL ==========
soricom_geo_model = gp.create_geomodel(
    project_name=SoricomModelConfig.PROJECT_NAME,
    extent=SoricomModelConfig.EXTENT,
    refinement=SoricomModelConfig.REFINEMENT,
    resolution=SoricomModelConfig.RESOLUTION,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=orientation_path,
        path_to_surface_points=formation_points_path,
    )
)

# Map geological formations including the fault
# The fault should be defined first in the structural frame
gp.map_stack_to_surfaces(
    gempy_model=soricom_geo_model,
    mapping_object=SoricomModelConfig.SURFACE_MAPPING
)

# Mark the Fault_Series as a fault
soricom_geo_model.structural_frame.structural_groups[SoricomModelConfig.FAULT_GROUP_INDEX].structural_relation = gp.data.StackRelationType.FAULT

# Set fault relations - the fault affects the lithological units
soricom_geo_model.structural_frame.fault_relations = SoricomModelConfig.FAULT_RELATIONS_MATRIX

# Compute the model
gp.compute_model(soricom_geo_model)

# ========== VISUALIZATION ==========
# 2D plot - multiple sections
p = gpv.plot_2d(
    soricom_geo_model,
    cell_number=SoricomModelConfig.VIZ_CELL_NUMBER,
    direction=SoricomModelConfig.VIZ_DIRECTION,
    show_lith=SoricomModelConfig.SHOW_LITH,
    show_boundaries=SoricomModelConfig.SHOW_BOUNDARIES,
    show_data=SoricomModelConfig.SHOW_DATA,
    legend=SoricomModelConfig.SHOW_LEGEND
)

Plotter.create_cross_section(
    soricom_geo_model,
    cross_section=SoricomModelConfig.VIZ_CROSS_SECTION_COUNT,
    vertical_exaggeration=SoricomModelConfig.VIZ_VERTICAL_EXAGGERATION
)
# 3D plot (commented out due to compatibility issue)
# gpv.plot_3d(soricom_geo_model, show_lith=True, ve=1)
