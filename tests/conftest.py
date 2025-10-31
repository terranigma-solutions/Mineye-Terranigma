import dotenv
import pyro
import torch

dotenv.load_dotenv()

import os
import pytest
import gempy as gp

seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)

@pytest.fixture(scope="session")
def base_dir():
    """Base directory for the project."""
    return os.getcwd()


@pytest.fixture(scope="session")
def geomodel_dir(base_dir):
    """Directory containing geomodel data."""
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'Output_Data'))


@pytest.fixture(scope="session")
def data_dir(base_dir):
    """Directory containing geomodel data."""
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'Input_Data'))


# topography_path = os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')
@pytest.fixture(scope="session")
def topography_dir(data_dir):
    """Directory containing topographic data."""
    return os.path.join(data_dir, 'Topographic_Data')


@pytest.fixture(scope="session")
def tmp_dir(geomodel_dir):
    """Temporary directory for model inputs."""
    return os.path.join(geomodel_dir, 'Simple-Models', 'temp_inputs')


@pytest.fixture(scope="session")
def geophysical_dir(data_dir):
    """Directory containing geophysical data."""
    return os.path.join(data_dir, 'Geophysical_Cleaned_Data')


@pytest.fixture(scope="session")
def model_paths(tmp_dir):
    """Paths to orientation and points CSV files."""
    return {
        'orientations': os.path.join(tmp_dir, 'orientations_mod.csv'),
        'points': os.path.join(tmp_dir, 'points_mod.csv')
    }


@pytest.fixture(scope="session")
def model_extent():
    """Model extent coordinates [min_x, max_x, min_y, max_y, min_z, max_z]."""
    min_x = -709521
    max_x = -675558
    min_y = 4526832
    max_y = 4551949
    max_z = 505
    model_depth = -500
    return [min_x, max_x, min_y, max_y, model_depth, max_z]


@pytest.fixture(scope="session")
def model_resolution():
    """Model resolution [nx, ny, nz]."""
    nx = ny = nz = 64
    return [nx, ny, nz]


@pytest.fixture
def simple_geo_model(model_extent, model_resolution, model_paths):
    """Factory for creating simple geological models with custom parameters."""
    def _create_model(project_name='simple_model', refinement=4, extent=None, resolution=None):
        geo_model = gp.create_geomodel(
            project_name=project_name,
            extent=extent or model_extent,
            refinement=refinement,
            resolution=resolution or model_resolution,
            importer_helper=gp.data.ImporterHelper(
                path_to_orientations=model_paths['orientations'],
                path_to_surface_points=model_paths['points'],
            )
        )

        gp.map_stack_to_surfaces(
            gempy_model=geo_model,
            mapping_object={
                    "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
            }
        )

        return geo_model

    return _create_model()