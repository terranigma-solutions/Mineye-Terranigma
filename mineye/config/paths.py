"""Shared path utilities for the project."""
import os


def get_base_dir():
    """Get the base directory for the project."""
    return os.getcwd()


def get_data_dir(base_dir=None):
    """Get the directory containing input data."""
    if base_dir is None:
        base_dir = get_base_dir()
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'Input_Data'))


def get_geomodel_dir(base_dir=None):
    """Get the directory containing geomodel output data."""
    if base_dir is None:
        base_dir = get_base_dir()
    return os.path.abspath(os.path.join(base_dir, 'examples', 'Data', 'Output_Data'))


def get_topography_dir(base_dir=None):
    """Get the directory containing topographic data."""
    data_dir = get_data_dir(base_dir)
    return os.path.join(data_dir, 'Topographic_Data')


def get_tmp_dir(base_dir=None):
    """Get the temporary directory for model inputs."""
    geomodel_dir = get_geomodel_dir(base_dir)
    tmp_dir = os.path.join(geomodel_dir, 'Simple-Models', 'temp_inputs')
    os.makedirs(tmp_dir, exist_ok=True)
    return tmp_dir


def get_geophysical_dir(base_dir=None):
    """Get the directory containing geophysical data."""
    data_dir = get_data_dir(base_dir)
    return os.path.join(data_dir, 'Geophysical_Cleaned_Data')


def get_model_paths(base_dir=None):
    """Get paths to orientation and points CSV files."""
    tmp_dir = get_tmp_dir(base_dir)
    return {
        'orientations': os.path.join(tmp_dir, 'orientations_mod.csv'),
        'points': os.path.join(tmp_dir, 'points_mod.csv')
    }


def get_model_extent():
    """Get model extent coordinates [min_x, max_x, min_y, max_y, min_z, max_z]."""
    min_x = -709521
    max_x = -675558
    min_y = 4526832
    max_y = 4551949
    max_z = 505
    model_depth = -500
    return [min_x, max_x, min_y, max_y, model_depth, max_z]


def get_model_resolution():
    """Get model resolution [nx, ny, nz]."""
    nx = ny = nz = 64
    return [nx, ny, nz]


def get_gis_map_info_path(base_dir=None):
    """Get path to the plutonite outline GeoPackage file."""
    data_dir = get_data_dir(base_dir)
    return os.path.join(data_dir, 'Stratigraphic_Data', 'QGIS', 'plutonite_outline.gpkg')


def get_contact_points_path(base_dir=None):
    """Get path to the contact points CSV file."""
    geomodel_dir = get_geomodel_dir(base_dir)
    return os.path.join(geomodel_dir, 'Simple-Models', 'contact_points.csv')


def get_topography_path(base_dir=None):
    """Get path to the reduced topography TIFF file."""
    data_dir = get_data_dir(base_dir)
    return os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')


def get_pickle_model_path(base_dir=None):
    """Get path to the pickled geological model file."""
    tmp_dir = get_tmp_dir(base_dir)
    return os.path.join(tmp_dir, 'simple_geo_model.pkl')


def get_orientations_path(base_dir=None):
    """Get path to the modified orientations CSV file."""
    tmp_dir = get_tmp_dir(base_dir)
    return os.path.join(tmp_dir, 'orientations_mod.csv')


def get_points_path(base_dir=None):
    """Get path to the modified points CSV file."""
    tmp_dir = get_tmp_dir(base_dir)
    return os.path.join(tmp_dir, 'points_mod.csv')


def get_gravity_data_path(base_dir=None):
    """Get path to the cleaned gravity data GeoJSON file."""
    geophysical_dir = get_geophysical_dir(base_dir)
    return os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson')


