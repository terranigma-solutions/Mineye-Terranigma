import sys
import os
import gempy as gp
import gempy_viewer as gpv
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig, SoricomSimpleModelConfig

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

def test_simple_model_sor():
    geo_model = gp.create_geomodel(
        project_name=SoricomSimpleModelConfig.PROJECT_NAME,
        extent=SoricomSimpleModelConfig.EXTENT,
        refinement=SoricomSimpleModelConfig.REFINEMENT,
        resolution=SoricomSimpleModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_soricom_orientations(),
            path_to_surface_points=paths.get_soricom_formation_points(),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=SoricomSimpleModelConfig.SURFACE_MAPPING
    )

    assert geo_model is not None

    gp.compute_model(geo_model)
    gpv.plot_2d(geo_model)
    gpv.plot_3d(geo_model, image=False)


def test_model_with_fault_sor():
    geo_model = gp.create_geomodel(
        project_name=SoricomModelConfig.PROJECT_NAME,
        extent=SoricomModelConfig.EXTENT,
        refinement=SoricomModelConfig.REFINEMENT,
        resolution=SoricomModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_soricom_orientations(),
            path_to_surface_points=paths.get_soricom_formation_points(),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=SoricomModelConfig.SURFACE_MAPPING
    )

    geo_model.structural_frame.structural_groups[SoricomModelConfig.FAULT_GROUP_INDEX].structural_relation = gp.data.StackRelationType.FAULT
    geo_model.structural_frame.fault_relations = SoricomModelConfig.FAULT_RELATIONS_MATRIX

    assert geo_model is not None

    gp.compute_model(geo_model)
    gpv.plot_2d(geo_model)
    gpv.plot_3d(geo_model, image=False)