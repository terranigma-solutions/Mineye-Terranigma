import pickle
import pandas as pd
import numpy as np
import os
import gempy as gp
import gempy_viewer as gpv

from mineye.GeoModel import HelperMethods, Plotter

# ========== CONFIG ==========
BASE_DIR = os.getcwd()
data_dir = os.path.abspath(os.path.join(BASE_DIR, '..', 'Data', 'Input_Data'))
geomodel_dir = os.path.abspath(os.path.join(BASE_DIR, '..',  'Data', 'Output_Data'))
forward_model_folder = os.path.abspath(os.path.join(BASE_DIR, '..',  'GeoModel', 'Geological_Forward_Modelling'))

gis_map_info = os.path.join(data_dir, 'Stratigraphic_Data', 'QGIS', 'plutonite_outline.gpkg')
points_df = pd.read_csv(os.path.join(geomodel_dir, 'Simple-Models', 'contact_points.csv'))
topography_path = os.path.join(data_dir, 'Topographic_Data', 'topo_reduced_sf0.1.tif')

tmp_dir = os.path.join(geomodel_dir, 'Simple-Models', 'temp_inputs')
os.makedirs(tmp_dir, exist_ok=True)
pickle_model_path = os.path.join(tmp_dir, 'simple_geo_model.pkl')
mod_or_path = os.path.join(tmp_dir, 'orientations_mod.csv')
mod_pts_path = os.path.join(tmp_dir, 'points_mod.csv')

# Load original full points once (for optional comparison)
original_points = pd.read_csv(os.path.join(geomodel_dir, 'Simple-Models', 'contact_points.csv'))
points_df = original_points.copy()

# Toggle showing original full dataset in the map plot
show_original_full_points = False
show_dip_labels = True  # Toggle to show/hide dip value annotations
show_scalebar = True    # Toggle to show/hide scalebar
scalebar_length = 5000  # Length of scalebar in coordinate units (assumed meters)

trigger_map_with_data_plot = True
trigger_geological_map_plot = True
trigger_create_cross_section = False
trigger_2d_plot = True
trigger_3d_plot = False

save_pickled_model = False
trigger_recreate_data = True

min_x = -709521
max_x = -675558
min_y = 4526832
max_y = 4551949
max_z = 505
model_depth = -500
simplification_level = 0.9  # 0=no simplification, 1=maximum simplification
dip_values = 10
azimuth_values = 0

# Set the extent of the model
extent=[min_x, max_x, min_y, max_y, model_depth, max_z]  # Fixed Y coordinate order

# Set model resolution
nx = ny = nz = 64
resolution = [nx, ny, nz]

# ========== DATA CLEANING ==========
if trigger_recreate_data:
    # TODO: This should be its own function with its own test
    # Generate orientations from contact points
    orientations_df = HelperMethods.generate_orientations_from_points(points_df, default_dip=dip_values, default_azimuth=azimuth_values, use_default_azimuth=False)
    orientations_df, points_df = HelperMethods.remove_boundary_artifacts(
        points_df, orientations_df, boundary_tolerance=800)

    # Apply simplification after cleaning - no formation filtering needed since there's only one formation
    _, points_df = HelperMethods.simplify_formation_data(
        orientations_df, points_df, 34, simplification_level)

# ========== MANIPULATE DATA MANUALLY ==========
    # Ensure an 'id' column exists for per-id edits
    if 'id' not in orientations_df.columns:
        orientations_df['id'] = np.arange(len(orientations_df))

    # Ensure contact points also have an 'id' column for plotting
    if 'id' not in points_df.columns:
        points_df['id'] = np.arange(len(points_df))

    # Optional: add manual orientations at specific contact point IDs
    manual_orientations_at_points = [
        11, 19
        # Example: 5, 12, 18  # Add orientations at these contact point IDs
    ]

    # Add manual orientations if specified
    if manual_orientations_at_points:
        orientations_df = HelperMethods.add_manual_orientations_at_points(
            orientations_df, points_df, manual_orientations_at_points, dip_values)

    # Optional: specify adjustments by orientation id (leave dicts empty to skip)
    azimuth_flip_by_id = {
        0: True, 1: True, 2: True, 3: True, 4: True, 5: False, 6: False
        # Example: 3: True, 7: True  # True = flip by 180°
    }
    manual_dip_by_id = {
        # Example: 3: 60, 7: 30  # degrees, will be clipped to [0, 90]
    }

    # Apply azimuth flips by id (180°)
    for oid, flip in (azimuth_flip_by_id or {}).items():
        if not flip:
            continue
        mask = orientations_df['id'] == oid
        orientations_df.loc[mask, 'azimuth'] = (orientations_df.loc[mask, 'azimuth'] + 180) % 360

    # Apply dip overrides by id
    for oid, dip in (manual_dip_by_id or {}).items():
        mask = orientations_df['id'] == oid
        try:
            dip_val = float(dip)
        except (TypeError, ValueError):
            continue
        orientations_df.loc[mask, 'dip'] = np.clip(dip_val, 0, 90)

# ========== SAVE TEMP FILES ==========
    orientations_df.to_csv(mod_or_path, index=False)
    points_df.to_csv(mod_pts_path, index=False)

else:
    # Load modified (reduced) data instead of using original full points
    orientations_df = pd.read_csv(mod_or_path)
    if os.path.exists(mod_pts_path):
        points_df = pd.read_csv(mod_pts_path)
    # Ensure IDs exist
    if 'id' not in orientations_df.columns:
        orientations_df['id'] = np.arange(len(orientations_df))
    if 'id' not in points_df.columns:
        points_df['id'] = np.arange(len(points_df))

# ========== PLOT DATA ==========
if trigger_map_with_data_plot and False: # ! @Max this is broken
    Plotter.plot_initial_data_on_map()


# ========== BUILD GEMPY MODEL ==========
simple_geo_model = gp.create_geomodel(
    project_name='simple_model',
    extent=extent,
    refinement=4,
    resolution=resolution,
    importer_helper=gp.data.ImporterHelper(
        path_to_orientations=mod_or_path,
        path_to_surface_points=mod_pts_path,
    )
)

gp.map_stack_to_surfaces(
    gempy_model=simple_geo_model,
    mapping_object={
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }
)

gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topography_path)
geo_model = gp.compute_model(simple_geo_model)

if save_pickled_model:
    with open(pickle_model_path, 'wb') as f:
        pickle.dump(simple_geo_model, f)

if trigger_2d_plot:
    p = gpv.plot_2d(
        simple_geo_model,
        section_names=['topography'],   # this triggers the top-down geological map
        show_topography=True,
        show_lith=True,
        show_boundaries=True,
        show_data=True,
        legend=False
    )

if trigger_create_cross_section:
    Plotter.create_cross_section(simple_geo_model, cross_section=5)

if trigger_3d_plot:
    gpv.plot_3d(simple_geo_model, section_names=['topography'], show_topography=True, ve=12)