"""Shared path utilities for the project."""
import os

def get_base_dir():
    """Get the base directory for the project (project root)."""
    this_dir = os.path.dirname(os.path.abspath(__file__))
    return os.path.abspath(os.path.join(this_dir, '..', '..'))


def get_data_dir():
    """Get the directory containing input data."""
    return os.path.abspath(os.path.join(get_base_dir(), 'examples', 'Data', 'Input_Data'))


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

def get_topography_path():
    """Get path to the reduced topography TIFF file."""
    data_dir = get_data_dir()
    return os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')

def get_pickle_model_path():
    """Get path to the pickled geological model file."""
    tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'simple_geo_model.pkl')

def get_orientations_path():
    """Get path to the modified orientations CSV file."""
    tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'orientations_mod.csv')


def get_points_path():
    """Get path to the modified points CSV file."""
    tmp_dir = get_tmp_dir()
    return os.path.join(tmp_dir, 'points_mod.csv')


def get_gravity_data_path():
    """Get path to the cleaned gravity data GeoJSON file."""
    geophysical_dir = get_geophysical_dir()
    return os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson')


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
