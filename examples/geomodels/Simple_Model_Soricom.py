import sys
import os
import gempy as gp
import gempy_viewer as gpv
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig, SoricomSimpleModelConfig, SoricomErosiveModelConfig

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))

def test_simple_model_sor():
    geo_model = gp.create_geomodel(
        project_name=SoricomSimpleModelConfig.PROJECT_NAME,
        extent=SoricomSimpleModelConfig.EXTENT,
        refinement=SoricomSimpleModelConfig.REFINEMENT,
        resolution=SoricomSimpleModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=os.path.join(paths.get_soricom_data_dir(), 'orientation_new.csv'),
            path_to_surface_points=os.path.join(paths.get_soricom_data_dir(), 'formation_points_new.csv'),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=SoricomSimpleModelConfig.SURFACE_MAPPING
    )

    assert geo_model is not None

    # Add topography
    gp.set_topography_from_file(
        grid=geo_model.grid,
        filepath=paths.get_soricom_dem_path()
    )

    gp.compute_model(geo_model)
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

    # Add topography
    gp.set_topography_from_file(
        grid=geo_model.grid,
        filepath=paths.get_soricom_dem_path()
    )

    gp.compute_model(geo_model)
    gpv.plot_3d(geo_model, image=False, show_topography=False)


def test_erosive_chromite_model_sor():
    geo_model = gp.create_geomodel(
        project_name=SoricomErosiveModelConfig.PROJECT_NAME,
        extent=SoricomErosiveModelConfig.EXTENT,
        refinement=SoricomErosiveModelConfig.REFINEMENT,
        resolution=SoricomErosiveModelConfig.RESOLUTION,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_soricom_erosive_orientations(),
            path_to_surface_points=paths.get_soricom_erosive_formation_points(),
        )
    )

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=SoricomErosiveModelConfig.SURFACE_MAPPING
    )

    # Set the fault relation
    geo_model.structural_frame.structural_groups[SoricomErosiveModelConfig.FAULT_GROUP_INDEX].structural_relation = gp.data.StackRelationType.FAULT

    # Set the chromite lense as erosive (intrusion)
    geo_model.structural_frame.structural_groups[SoricomErosiveModelConfig.CHROMITE_GROUP_INDEX].structural_relation = gp.data.StackRelationType.ERODE

    # Set fault relations matrix
    geo_model.structural_frame.fault_relations = SoricomErosiveModelConfig.FAULT_RELATIONS_MATRIX

    assert geo_model is not None

    # Add topography
    gp.set_topography_from_file(
        grid=geo_model.grid,
        filepath=paths.get_soricom_dem_path()
    )

    gp.compute_model(geo_model)

    gpv.plot_3d(geo_model, image=False)