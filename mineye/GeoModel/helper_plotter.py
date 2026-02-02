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
                        show_topography=True,
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
                        show_topography=True,
                        legend=True,
                        show_data=False,
                        direction="x")
        # Increment i to avoid infinite loop
        i += 1

def plot_3d_async(geo_model):
    gpv.plot_3d(geo_model, show_lith=True, show_boundaries=True, ve=10, legend=False, show_data=True,
                show_topography=True, transformed_data=False)


def plot_combined_model(
        lith_block: np.ndarray,
        voxel_coords: np.ndarray,
        formation_id_map: dict,
        formation_colors: dict,
        topography_points: np.ndarray = None,
        title: str = 'Combined Geological Model',
        sample_step: int = 10,
        figsize: tuple = (12, 10),
        elev: float = 90,
        azim: float = 270
):
    """
    Plot a combined 3D geological model with proper legend.

    Parameters
    ----------
    lith_block : np.ndarray
        Flattened lithology block array with formation IDs.
    voxel_coords : np.ndarray
        Array of voxel coordinates (N, 3) for X, Y, Z.
    formation_id_map : dict
        Mapping of lithology IDs to formation names, e.g. {1: 'Formation A', 2: 'Formation B'}.
    formation_colors : dict
        Mapping of formation names to colors, e.g. {'Formation A': '#e74c3c'}.
    topography_points : np.ndarray, optional
        Topography as array of (X, Y, Z) points. If provided, voxels above
        topography will be masked out.
    title : str, optional
        Plot title. Default is 'Combined Geological Model'.
    sample_step : int, optional
        Subsampling step for plotting (every nth point). Default is 10.
    figsize : tuple, optional
        Figure size. Default is (12, 10).
    elev : float, optional
        Elevation angle for 3D view. Default is 90.
    azim : float, optional
        Azimuth angle for 3D view. Default is 270.
    """
    from scipy.interpolate import LinearNDInterpolator

    # Create a copy of coords and lith_block for masking
    coords_to_plot = voxel_coords.copy()
    lith_to_plot = lith_block.copy()

    # Apply topography mask if provided
    if topography_points is not None:
        # Create interpolator from topography X, Y -> Z
        topo_xy = topography_points[:, :2]  # X, Y
        topo_z = topography_points[:, 2]     # Z elevations

        interp = LinearNDInterpolator(topo_xy, topo_z)

        # Get topography elevation at each voxel's X,Y position
        topo_at_voxels = interp(voxel_coords[:, :2])

        # Create mask: keep voxels that are below or at topography
        # Handle NaN values (outside interpolation range) by keeping those voxels
        below_topo_mask = (voxel_coords[:, 2] <= topo_at_voxels) | np.isnan(topo_at_voxels)

        # Apply mask
        coords_to_plot = voxel_coords[below_topo_mask]
        lith_to_plot = lith_block[below_topo_mask]

    # Flip Y coordinate for correct orientation
    y_min, y_max = np.min(coords_to_plot[:, 1]), np.max(coords_to_plot[:, 1])
    voxel_coords_flipped = coords_to_plot.copy()
    voxel_coords_flipped[:, 1] = y_max - (coords_to_plot[:, 1] - y_min)

    # Create 3D plot
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')

    # Plot each formation separately to create proper legend
    for lith_id, formation_name in formation_id_map.items():
        mask = lith_to_plot[::sample_step] == lith_id
        if np.any(mask):
            ax.scatter(
                voxel_coords_flipped[:, 0][::sample_step][mask],
                voxel_coords_flipped[:, 1][::sample_step][mask],
                voxel_coords_flipped[:, 2][::sample_step][mask],
                c=formation_colors.get(formation_name, '#888888'),
                s=15,
                alpha=0.7,
                label=formation_name
            )

    # Add legend with formation names - positioned below the figure
    legend = ax.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.05),
        fontsize=14,
        framealpha=0.95,
        edgecolor='black',
        fancybox=True,
        shadow=True,
        title='Formations',
        title_fontsize=16,
        markerscale=2.5,
        ncol=2  # Arrange in 2 columns for horizontal layout
    )
    legend.get_title().set_fontweight('bold')

    # Set labels and title with padding
    ax.set_xlabel('X Coordinate', fontsize=12, labelpad=10)
    ax.set_ylabel('Y Coordinate', fontsize=12, labelpad=10)
    ax.set_zlabel('Z Elevation', fontsize=12, labelpad=10)
    ax.set_title(title, fontsize=16, fontweight='bold', pad=15)

    # Reduce number of ticks to avoid clutter
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.yaxis.set_major_locator(plt.MaxNLocator(5))
    ax.zaxis.set_major_locator(plt.MaxNLocator(5))

    # Reduce tick label size
    ax.tick_params(axis='x', labelsize=9, pad=2)
    ax.tick_params(axis='y', labelsize=9, pad=2)
    ax.zaxis.set_tick_params(labelsize=9, pad=2)

    # Set view angle
    ax.view_init(elev=elev, azim=azim)
    plt.subplots_adjust(bottom=0.15)  # Make room for legend below
    plt.show()


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


def plot_model_and_gravity_sensors(simple_geo_model, xy_ravel, grav):
    import pyvista as pv
    p3d = gpv.plot_3d(
        simple_geo_model,
        ve=5,
        image=False,
        show=False,
        show_boundaries=True,
        show_lith=False,
        show_octree=False,
        kwargs_pyvista_bounds={
                'show_xlabels': False,
                'show_ylabels': False,
                'show_zlabels': False
        }
    )

    # Create a PolyData object from the detector coordinates
    detectors_poly = pv.PolyData(xy_ravel)
    detectors_poly.point_data['gravity'] = grav

    # Define the inverted pyramid geometry (Cone with 4 faces)
    # direction=(0, 0, -1) points the pyramid downwards (inverted)
    # resolution=4 creates a square base (pyramid)
    # center=(0, 0, 600) shifts the pyramid up so its tip is at the data point (offset = height/2)
    # Adjust radius and height based on your scene scale (e.g., 600m radius, 1200m height)
    pyramid_source = pv.Cone(
        center=(0, 0, 0),
        radius=600,
        height=100,
        direction=(0, 0, -1),
        resolution=4
    )

    # Glyph the points with the pyramid source
    # scale=False ensures all markers are the same size
    glyphs = detectors_poly.glyph(geom=pyramid_source, scale=False, orient=False)

    p3d.p.add_mesh(
        glyphs,
        scalars='gravity',
        cmap='viridis_r',
        name='gravity_observations'
    )
    p3d.p.show()
