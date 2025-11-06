import sys
import os
import gempy as gp
import gempy_viewer as gpv
from mineye.GeoModel import Plotter
from mineye.config import paths
from mineye.config.example_parameters import TharsisModelConfig

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))


def create_simple_model():
    """
    Create and compute a simple Tharsis geological model without topography or visualization.

    Returns:
        geo_model: Computed GemPy geological model
    """
    geo_model = gp.create_geomodel(
        project_name=TharsisModelConfig.PROJECT_NAME,
        extent=TharsisModelConfig.EXTENT,
        refinement=TharsisModelConfig.REFINEMENT,
        resolution=TharsisModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_orientations_path(),
            path_to_surface_points=paths.get_points_path(),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=TharsisModelConfig.SURFACE_MAPPING
    )

    gp.set_topography_from_file(grid=geo_model.grid, filepath=paths.get_topography_path())
    gp.compute_model(geo_model)

    return geo_model


def test_simple_model_with_plots():
    geo_model = gp.create_geomodel(
        project_name=TharsisModelConfig.PROJECT_NAME,
        extent=TharsisModelConfig.EXTENT,
        refinement=TharsisModelConfig.REFINEMENT,
        resolution=TharsisModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_orientations_path(),
            path_to_surface_points=paths.get_points_path(),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=TharsisModelConfig.SURFACE_MAPPING
    )

    gp.set_topography_from_file(grid=geo_model.grid, filepath=paths.get_topography_path())
    gp.compute_model(geo_model)

    # ========== VISUALIZATION ==========
    gpv.plot_2d(
        geo_model,
        section_names=['topography'],
        show_topography=TharsisModelConfig.SHOW_TOPOGRAPHY,
        show_lith=TharsisModelConfig.SHOW_LITH,
        show_boundaries=TharsisModelConfig.SHOW_BOUNDARIES,
        show_data=TharsisModelConfig.SHOW_DATA,
        legend=TharsisModelConfig.SHOW_LEGEND
    )

    Plotter.create_cross_section(
        geo_model,
        cross_section=TharsisModelConfig.VIZ_CROSS_SECTION_COUNT,
        vertical_exaggeration=TharsisModelConfig.VIZ_VERTICAL_EXAGGERATION
    )

    gpv.plot_3d(
        geo_model,
        section_names=['topography'],
        show_topography=TharsisModelConfig.SHOW_TOPOGRAPHY,
        ve=TharsisModelConfig.VIZ_3D_VERTICAL_EXAGGERATION
    )
