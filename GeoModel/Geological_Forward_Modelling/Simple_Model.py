from io import BufferedWriter

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import pickle
import HelperMethods
import gempy as gp
import gempy_viewer as gpv
import geopandas as gpd

# ========== CONFIG ==========
BASE_DIR = os.getcwd()
data_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Input_Data'))
geomodel_dir = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'Data', 'Output_Data'))
forward_model_folder = os.path.abspath(os.path.join(BASE_DIR, '..', '..', 'GeoModel', 'Geological_Forward_Modelling'))

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
trigger_recreate_data = False

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
    manual_orientation_flip_azimuth = True  # Set to True to flip azimuth by 180° for manually added orientations

    # Add manual orientations if specified
    if manual_orientations_at_points:
        orientations_df = HelperMethods.add_manual_orientations_at_points(
            orientations_df, points_df, manual_orientations_at_points, dip_values, manual_orientation_flip_azimuth
        )

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
if trigger_map_with_data_plot:
    fig, ax = plt.subplots(1, 1, figsize=(12, 8))

    # --- Add plutonite geological map as background (combined figure) ---
    if os.path.exists(gis_map_info):
        try:
            plutonite_gdf = gpd.read_file(gis_map_info)
            plutonite_color = (125/255, 159/255, 183/255)  # RGB(125,159,183)
            plutonite_gdf.plot(
                ax=ax,
                facecolor=plutonite_color,
                edgecolor='white',  # subtle outline for contrast
                alpha=0.95,
                linewidth=0.6,
                zorder=0
            )
        except Exception as e:
            print(f"Warning: Could not plot plutonite map: {e}")

    if show_original_full_points:
        ax.scatter(original_points['X'], original_points['Y'], c='gray', s=16, alpha=0.35, marker='x', label='Original Points')
        original_orientations = HelperMethods.generate_orientations_from_points(original_points, default_dip=45, default_azimuth=90)
        ax.scatter(original_orientations['X'], original_orientations['Y'], c='gray', s=20, alpha=0.25, label='Original Orientations')

    # Plot ONLY the reduced/cleaned points (current working dataset)
    ax.scatter(points_df['X'], points_df['Y'], c='darkred', s=42, alpha=0.9, marker='x', label='Reduced Contact Points', zorder=5)

    # Orientation arrows from current orientations_df
    arrow_length = 1000
    azimuth_rad = np.deg2rad(orientations_df['azimuth'])
    dx = np.sin(azimuth_rad) * arrow_length
    dy = np.cos(azimuth_rad) * arrow_length
    ax.quiver(orientations_df['X'], orientations_df['Y'], dx, dy,
              angles='xy', scale_units='xy', scale=1, color='navy', width=0.0028, headwidth=3, headlength=4, zorder=6, label='Orientations')

    # Orientation IDs
    for _, row in orientations_df.iterrows():
        ax.text(row['X'] - 300, row['Y'] + 300, f"ID:{int(row['id'])}",
                color='red', fontsize=8, fontweight='bold',
                bbox=dict(boxstyle="round,pad=0.22", facecolor='white', alpha=0.85, edgecolor='red'),
                zorder=7)

    # Dip labels (optional)
    if show_dip_labels:
        for _, row in orientations_df.iterrows():
            azimuth_rad = np.deg2rad(row['azimuth'])
            perpendicular_offset = 600
            text_x = row['X'] + 0.5 * arrow_length * np.sin(azimuth_rad) - perpendicular_offset * np.cos(azimuth_rad)
            text_y = row['Y'] + 0.5 * arrow_length * np.cos(azimuth_rad) + perpendicular_offset * np.sin(azimuth_rad)
            ax.text(text_x, text_y, f"{row['dip']:.0f}°", color='navy', fontsize=9, ha='center', va='center', fontweight='bold', zorder=7,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', alpha=0.65, edgecolor='navy'))

    # Contact point IDs
    for _, row in points_df.iterrows():
        if 'id' in points_df.columns:
            ax.text(row['X'] - 300, row['Y'] + 300, f"ID:{int(row['id'])}",
                    color='darkred', fontsize=8, fontweight='bold',
                    bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8, edgecolor='darkred'),
                    zorder=7)

    # === Beautify ===
    # Build dynamic title including reduction level and point counts
    try:
        reduction_pct = simplification_level * 100
        title_suffix = f"Reduction Level: {reduction_pct:.0f}% | Points: {len(points_df)}/{len(original_points)}"
    except Exception:
        title_suffix = f"Reduction Level: {simplification_level}"
    ax.set_title(f'Reduced Geological Input Data (Contacts & Orientations) – {title_suffix}', fontsize=14, fontweight='bold', pad=12)
    # Remove axis label names but keep ticks/gridlines
    ax.set_xlabel('')
    ax.set_ylabel('')

    ax.set_aspect('equal')
    ax.grid(False)

    # Slim legend
    leg = ax.legend(frameon=True, fontsize=9, loc='upper right')
    if leg:
        leg.get_frame().set_alpha(0.85)
        leg.get_frame().set_facecolor('white')
        leg.get_frame().set_edgecolor('#444444')

    # Add scalebar (lower right)
    if show_scalebar:
        try:
            xmin, xmax = ax.get_xlim(); ymin, ymax = ax.get_ylim()
            sb_len = scalebar_length
            # Place with 2% padding
            x_pad = 0.02 * (xmax - xmin)
            y_pad = 0.02 * (ymax - ymin)
            x_start = xmax - sb_len - x_pad
            x_end = xmax - x_pad
            y_bar = ymin + y_pad + 0.01 * (ymax - ymin)
            ax.plot([x_start, x_end], [y_bar, y_bar], color='k', linewidth=3, solid_capstyle='butt', zorder=10)
            ax.plot([x_start, x_start], [y_bar - 0.004*(ymax-ymin), y_bar + 0.004*(ymax-ymin)], color='k', linewidth=2, zorder=10)
            ax.plot([x_end, x_end], [y_bar - 0.004*(ymax-ymin), y_bar + 0.004*(ymax-ymin)], color='k', linewidth=2, zorder=10)
            label_txt = f"{int(sb_len/1000)} km" if sb_len >= 1000 else f"{int(sb_len)} m"
            ax.text((x_start + x_end)/2, y_bar + 0.012*(ymax-ymin), label_txt,
                    ha='center', va='bottom', fontsize=9, fontweight='bold', color='k', zorder=11,
                    bbox=dict(boxstyle='round,pad=0.15', facecolor='white', edgecolor='k', alpha=0.8))
        except Exception as e:
            print(f"Scalebar failed: {e}")

    plt.tight_layout(); plt.show()

# ========== BUILD GEMPY MODEL (==========
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
    HelperMethods.create_cross_section(simple_geo_model, cross_section=5)

if trigger_3d_plot:
    gpv.plot_3d(simple_geo_model, section_names=['topography'], show_topography=True, ve=12)