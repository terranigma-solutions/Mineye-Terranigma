import dotenv
import pandas as pd
orig_read_csv = pd.read_csv
pd.read_csv = lambda *args, **kwargs: orig_read_csv(*args, **kwargs).apply(lambda x: x.astype(object) if isinstance(x.array, pd.arrays.StringArray) else x)

import pyro
import torch

dotenv.load_dotenv()
import os
os.environ["VALIDATE_SERIALIZATION"] = ""

seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)

import pytest

# --- GemPy-specific patches (may not be available in all envs) ---
try:
    import gempy.core.data.encoders.json_geomodel_encoder as jge
    def robust_encode_numpy_array(array):
        size = array.size() if callable(getattr(array, 'size', None)) else getattr(array, 'size', 0)
        if isinstance(size, (list, tuple)):
            import numpy as np
            size = int(np.prod(size))
        if size > 10:
            return []
        if hasattr(array, 'tolist'):
            try:
                return array.tolist()
            except Exception:
                pass
        return []
    jge.encode_numpy_array = robust_encode_numpy_array
    import gempy.core.data.geo_model as gm
    gm.encode_numpy_array = robust_encode_numpy_array
    if hasattr(gm.GeoModel, 'model_config') and gm.GeoModel.model_config:
        if 'json_encoders' in gm.GeoModel.model_config:
            import numpy as np
            gm.GeoModel.model_config['json_encoders'][np.ndarray] = robust_encode_numpy_array
    import gempy as gp
    _HAS_GEMPY = True
except ImportError:
    gp = None
    _HAS_GEMPY = False

@pytest.fixture(scope="session")
def base_dir():
    """Base directory for the project."""
    conftest_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(conftest_dir, '..', '..'))
    return project_root


@pytest.fixture(scope="session")
def geomodel_dir(base_dir):
    """Directory containing geomodel data."""
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'Model_Input_Data'))


@pytest.fixture(scope="session")
def data_dir(base_dir):
    """Directory containing geomodel data."""
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'General_Input_Data'))


# topography_path = os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')
@pytest.fixture(scope="session")
def topography_dir(data_dir):
    """Directory containing topographic data."""
    return os.path.join(data_dir, 'Topographic_Data')


@pytest.fixture(scope="session")
def tmp_dir(geomodel_dir):
    """Temporary directory for model inputs."""
    return os.path.join(geomodel_dir, 'Simple-Models')


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
    # min_x = -707521
    min_x = -707_521  # * Cropping the corrupted area of the geotiff 
    max_x = -675558
    min_y = 4_526_832
    max_y = 4_551_949
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
    model_resolution =None # ! Let's use octrees 
    
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