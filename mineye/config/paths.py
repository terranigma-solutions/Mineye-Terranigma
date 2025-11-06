"""Shared path utilities for the project."""
import os

# Define input paths similar to GemPy's conftest pattern
# Primary path: relative to this file (for installed/development use)
input_path = os.path.dirname(__file__) + '/../../examples/Data/Input_Data'
# Alternative path: for test/example contexts
input_path2 = os.path.dirname(__file__) + '/../../examples/Data/Input_Data/'

def get_base_dir():
    """Get the base directory for the project (project root)."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(this_dir, '..', '..'))


def get_data_dir():
    """Get the directory containing input data.

    Returns the first existing path from the configured input paths.
    """
    # Try primary path first
    primary = os.path.abspath(input_path)
    if os.path.exists(primary):
        return primary

    # Fall back to base_dir construction
    fallback = os.path.abspath(os.path.join(get_base_dir(), 'examples', 'Data', 'Input_Data'))
    if os.path.exists(fallback):
        return fallback

    # Return primary even if it doesn't exist (will be created or raise error appropriately)
    return primary


def get_geomodel_dir():
    return os.path.abspath(os.path.join(get_base_dir(), 'examples', 'Data', 'Output_Data'))


def get_topography_dir():
    """Get the directory containing topographic data."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Topographic_Data')


def get_tmp_dir():
    """Get the temporary directory for model inputs."""
    geomodel_dir = get_geomodel_dir()
    tmp_dir = os.path.join(geomodel_dir, 'Simple-Models', 'temp_inputs')
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir

def get_geophysical_dir():
    """Get the directory containing geophysical data."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Geophysical_Cleaned_Data')

def get_model_paths():
    """Get paths to orientation and points CSV files."""
    tmp_dir = get_tmp_dir()
    return {
        'orientations': os.path.join(tmp_dir, 'orientations_mod.csv'),
        'points': os.path.join(tmp_dir, 'points_mod.csv')
    }

def get_gis_map_info_path():
    """Get path to the plutonite outline GeoPackage file."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Stratigraphic_Data', 'QGIS', 'plutonite_outline.gpkg')

def get_contact_points_path():
    """Get path to the contact points CSV file."""
    geomodel_dir = get_geomodel_dir()
    return os.path.join(geomodel_dir, 'Simple-Models', 'contact_points.csv')

def get_topography_path(base_dir=None):
    """Get path to the reduced topography TIFF file.

    Args:
        base_dir: Optional base directory (for compatibility). If None, uses get_data_dir().
    """
    if base_dir is not None:
        # Legacy compatibility: base_dir points to project root
        data_dir = os.path.join(base_dir, 'examples', 'Data', 'Input_Data')
    else:
        data_dir = get_data_dir()
    return os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')

def get_pickle_model_path():
    """Get path to the pickled geological model file."""
    tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'simple_geo_model.pkl')

def get_orientations_path(base_dir=None):
    """Get path to the modified orientations CSV file.

    Args:
        base_dir: Optional base directory (for compatibility). If None, uses get_tmp_dir().
    """
    if base_dir is not None:
        # Legacy compatibility: base_dir points to project root
        tmp_dir = os.path.join(base_dir, 'examples', 'Data', 'Output_Data', 'Simple-Models', 'temp_inputs')
    else:
        tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'orientations_mod.csv')


def get_points_path(base_dir=None):
    """Get path to the modified points CSV file.

    Args:
        base_dir: Optional base directory (for compatibility). If None, uses get_tmp_dir().
    """
    if base_dir is not None:
        # Legacy compatibility: base_dir points to project root
        tmp_dir = os.path.join(base_dir, 'examples', 'Data', 'Output_Data', 'Simple-Models', 'temp_inputs')
    else:
        tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'points_mod.csv')


def get_gravity_data_path():
    """Get path to the cleaned gravity data GeoJSON file."""
    geophysical_dir = get_geophysical_dir()
    return os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson')


def get_plutonic_data_dir():
    """Get the directory containing plutonic data."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Plutonic_Data')


def get_plutonic_contact_points_path():
    """Get path to the original plutonic contact points CSV file."""
    plutonic_dir = get_plutonic_data_dir()
    return os.path.join(plutonic_dir, 'plutonic_contact_points.csv')


def get_soricom_data_dir():
    """Get the directory containing Soricom data."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Soricom_Data')

def get_soricom_orientations():
    """Get path to the Soricom orientation CSV file with fault data included."""
    soricom_dir = get_soricom_data_dir()
    return os.path.join(soricom_dir, 'orientations_with_fault.csv')

def get_soricom_formation_points():
    """Get path to the Soricom formation points CSV file with fault data included."""
    soricom_dir = get_soricom_data_dir()
    return os.path.join(soricom_dir, 'formation_points_with_fault.csv')
