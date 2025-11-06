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


def plot_forward_gravity_model(xy_ravel, grav, gravity_data, gravity_resolution=20, use_actual_locations=True):
    """
    Plot forward gravity model results.

    Args:
        xy_ravel: Array of measurement point coordinates
        grav: Computed gravity values
        gravity_data: GeoDataFrame with actual measurement locations
        gravity_resolution: Grid resolution (used if not using actual locations)
        use_actual_locations: Whether to use actual measurement locations or regular grid
    """
    fig, ax = plt.subplots(figsize=(12, 10))

    if use_actual_locations:
        scatter = ax.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=grav, s=30,
                            cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
        cbar = plt.colorbar(scatter, ax=ax)
        cbar.set_label(r'Forward Model Gravity ($\mu$gal)')
        print(f"Plotting {len(xy_ravel)} actual measurement locations")
    else:
        grav_res = gravity_resolution
        ax.scatter(xy_ravel[:, 0], xy_ravel[:, 1], s=1, c='white', alpha=0.5)
        im = ax.imshow(grav.reshape(grav_res, grav_res),
                   extent=(xy_ravel[:, 0].min() + (xy_ravel[0, 0] - xy_ravel[1, 0]) / 2,
                           xy_ravel[:, 0].max() - (xy_ravel[0, 0] - xy_ravel[1, 0]) / 2,
                           xy_ravel[:, 1].min() + (xy_ravel[0, 1] - xy_ravel[grav_res, 1]) / 2,
                           xy_ravel[:, 1].max() - (xy_ravel[0, 1] - xy_ravel[grav_res, 1]) / 2),
                   cmap='viridis_r', origin='lower', alpha=.8)
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label(r'Forward Model Gravity ($\mu$gal)')
        print(f"Plotting regular grid: {grav_res}x{grav_res} points")

    # Show actual measurement point locations
    actual_measurement_points = np.column_stack([
        np.array(gravity_data.geometry.x.values),
        np.array(gravity_data.geometry.y.values)
    ])
    ax.scatter(actual_measurement_points[:, 0], actual_measurement_points[:, 1],
               s=15, c='red', marker='x', alpha=0.8, linewidth=1.5,
               label=f'Actual Measurement Points (n={len(actual_measurement_points)})')

    ax.legend(loc='upper right', framealpha=0.9)
    ax.set_title('Forward Gravity Model Results', fontsize=14, fontweight='bold')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.show()


def plot_gravity_comparison(xy_ravel, observed_norm, forward_norm, residuals_norm,
                            unit_label, normalize_data=True, normalization_method='minmax'):
    """
    Plot comparison between observed and forward model gravity.

    Args:
        xy_ravel: Array of measurement point coordinates
        observed_norm: Normalized observed gravity values
        forward_norm: Normalized forward model gravity values
        residuals_norm: Normalized residuals (observed - forward)
        unit_label: Label for the units being used
        normalize_data: Whether data was normalized
        normalization_method: Method used for normalization

    Returns:
        correlation: Correlation coefficient between observed and forward model
    """
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

    # Plot 1: Observed gravity
    scatter1 = ax1.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=observed_norm,
                          s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    title_suffix = f" ({normalization_method})" if normalize_data else ""
    ax1.set_title(f'Observed Gravity{title_suffix}', fontsize=12, fontweight='bold')
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_aspect('equal')
    cbar1 = plt.colorbar(scatter1, ax=ax1)
    cbar1.set_label(f'Observed ({unit_label})')

    # Plot 2: Forward model gravity
    scatter2 = ax2.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=forward_norm,
                          s=30, cmap='viridis_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax2.set_title(f'Forward Model Gravity{title_suffix}', fontsize=12, fontweight='bold')
    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.set_aspect('equal')
    cbar2 = plt.colorbar(scatter2, ax=ax2)
    cbar2.set_label(f'Forward Model ({unit_label})')

    # Plot 3: Residuals
    scatter3 = ax3.scatter(xy_ravel[:, 0], xy_ravel[:, 1], c=residuals_norm,
                          s=30, cmap='RdBu_r', alpha=0.8, edgecolors='black', linewidth=0.5)
    ax3.set_title(f'Residuals (Observed - Forward Model){title_suffix}', fontsize=12, fontweight='bold')
    ax3.set_xlabel('X (m)')
    ax3.set_ylabel('Y (m)')
    ax3.set_aspect('equal')
    cbar3 = plt.colorbar(scatter3, ax=ax3)
    cbar3.set_label(f'Residuals ({unit_label})')

    # Plot 4: Correlation plot
    ax4.scatter(observed_norm, forward_norm, alpha=0.7, s=40, edgecolors='black', linewidth=0.5)
    ax4.set_xlabel(f'Observed ({unit_label})')
    ax4.set_ylabel(f'Forward Model ({unit_label})')
    ax4.set_title('Observed vs Forward Model Correlation', fontsize=12, fontweight='bold')

    # Add 1:1 line
    lims = [min(ax4.get_xlim()[0], ax4.get_ylim()[0]), max(ax4.get_xlim()[1], ax4.get_ylim()[1])]
    ax4.plot(lims, lims, 'r--', alpha=0.75, linewidth=2, label='1:1 line')
    ax4.legend()
    ax4.grid(True, alpha=0.3)

    # Calculate and display correlation coefficient
    correlation = np.corrcoef(observed_norm, forward_norm)[0, 1]
    ax4.text(0.05, 0.95, f'R = {correlation:.3f}', transform=ax4.transAxes,
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.show()

    return correlation

