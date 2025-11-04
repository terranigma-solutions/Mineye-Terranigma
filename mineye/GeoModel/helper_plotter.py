import gempy_viewer as gpv
import matplotlib.pyplot as plt
import geopandas as gpd
import os
import numpy as np
import pandas as pd
from . import helper_methods


def create_cross_section(geo_model, cross_section: int, vertical_exaggeration: int):
    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=vertical_exaggeration,
                        cell_number=i,
                        show_topography=False,
                        legend=True,
                        show_data=False,
                        direction="y")
        # Increment i to avoid infinite loop
        i += 1

    i = 0
    while i < cross_section:
        if i % 2 == 0:
            gpv.plot_2d(geo_model, ve=vertical_exaggeration,
                        cell_number=i,
                        show_topography=False,
                        legend=True,
                        show_data=False,
                        direction="x")
        # Increment i to avoid infinite loop
        i += 1

def plot_3d_async(geo_model):
    gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=10, legend=False, show_data=True,
                show_topography=True, transformed_data=False)

def plot_initial_data_on_map(gis_map_info: str,
                             original_points: str,
                             points_df: pd.DataFrame,
                             orientations_df: pd.DataFrame,
                             show_original_full_points: bool = False,
                             show_dip_labels: bool = True,
                             show_scalebar: bool = True,
                             simplification_level: int = 0.9,
                             scalebar_length: float = 5000,
                             ):

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
        original_orientations = helper_methods.generate_orientations_from_points(original_points, default_dip=45, default_azimuth=90)
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
