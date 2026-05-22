"""Test reading and visualizing raw magnetic raster data for the Ternove model."""

import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio


def test_read_magnetic_raster(data_dir, geophysical_dir):
    """Read the merged/masked B1B2 magnetic raster and plot it alongside the cleaned point data."""

    # Path to the raw .ers raster (TMI, upward-continued 3m)
    raster_path = os.path.join(
        data_dir,
        'Geophysical_Raw_Data',
        'merged_masked_B1B2_mean_median_201_up3m.ers',
    )

    assert os.path.exists(raster_path), f"Raster not found: {raster_path}"

    # --- Read DEM ---
    dem_path = os.path.join(data_dir, 'Topographic_Data', 'soricom_DEM10m.tif')
    assert os.path.exists(dem_path), f"DEM not found: {dem_path}"

    with rasterio.open(dem_path) as dem_src:
        dem_crs = dem_src.crs
        dem_bounds = dem_src.bounds
        dem_data = dem_src.read(1)

    # DEM nodata mask (soricom DEM uses 0 for nodata)
    dem_valid = dem_data > 0
    dem_masked = np.where(dem_valid, dem_data, np.nan)

    print(f"\nDEM shape: {dem_data.shape}")
    print(f"DEM CRS: {dem_crs}")
    print(f"DEM bounds: {dem_bounds}")
    print(f"DEM elevation range: {np.nanmin(dem_masked):.1f} – {np.nanmax(dem_masked):.1f} m")

    # --- Read magnetic raster ---
    with rasterio.open(raster_path) as src:
        crs = src.crs
        bounds = src.bounds
        transform = src.transform
        data = src.read(1)

    # Mask nodata
    valid_mask = data != -999999
    masked = np.where(valid_mask, data, np.nan)

    # Extract all pixel coordinates and values
    rows, cols = np.where(valid_mask)
    lon, lat = rasterio.transform.xy(transform, rows, cols)
    lon = np.array(lon)
    lat = np.array(lat)
    valid_vals = data[valid_mask]

    print(f"\nRaster shape: {data.shape}")
    print(f"Raster CRS (declared): {crs}")
    print(f"Raster bounds: {bounds}")
    print(f"Valid pixels: {valid_mask.sum():,} / {data.size:,}")
    print(f"Magnetic range: {np.nanmin(masked):.2f} – {np.nanmax(masked):.2f} nT")
    print(f"Mean: {np.nanmean(masked):.2f} nT, Std: {np.nanstd(masked):.2f} nT")

    # Note: the .ers file declares EPSG:4326 but stores UTM eastings/northings
    #   (441878–442444 mE, 4586162–4586735 mN) — treat as EPSG:32634 for overlay.
    raster_crs_actual = 'EPSG:32634'

    # --- Read magnetic susceptibility lab data ---
    sus_path = os.path.join(data_dir, 'Soricom_Data', 'magnetic_susceptibility.csv')
    if os.path.exists(sus_path):
        sus_df = pd.read_csv(sus_path)
        print(f"\nMagnetic susceptibility data: {len(sus_df)} samples")
        print(f"  sus_SI range: {sus_df['sus_SI'].min():.6f} – {sus_df['sus_SI'].max():.6f}")
        print(f"  Mean: {sus_df['sus_SI'].mean():.6f}, Median: {sus_df['sus_SI'].median():.6f}")
        print(f"  By lithology:")
        for litho, grp in sus_df.groupby('lithology'):
            if litho:
                print(f"    {litho:25s}: n={len(grp):2d}, mean={grp['sus_SI'].mean():.6f}, "
                      f"min={grp['sus_SI'].min():.6f}, max={grp['sus_SI'].max():.6f}")
    else:
        sus_df = None
        print(f"\nSusceptibility data not found at: {sus_path}")

    # --- Read cleaned magnetic point data ---
    cleaned_path = os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson')

    if os.path.exists(cleaned_path):
        points = gpd.read_file(cleaned_path)
        print(f"\nCleaned point data: {len(points)} points")
        print(f"Point data CRS: {points.crs}")
        print(f"Point MAG range: {points['MAG'].min():.2f} – {points['MAG'].max():.2f} nT")

        # Reproject points to actual raster CRS (EPSG:32634)
        if points.crs is not None and points.crs.to_string() != raster_crs_actual:
            points = points.to_crs(raster_crs_actual)
    else:
        points = None
        print(f"\nCleaned point data not found at: {cleaned_path}")

    # --- Subsample raster points for scatter plot ---
    n_sample = min(5000, len(lon))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(lon), size=n_sample, replace=False)

    # --- Read borehole locations from formation points ---
    borehole_path = os.path.join(data_dir, 'Soricom_Data', 'formation_points_with_fault.csv')
    if os.path.exists(borehole_path):
        bh_df = pd.read_csv(borehole_path)
        # Borehole collars = host_rock entries (each is a unique collar location)
        collars = bh_df[bh_df['formation'] == 'host_rock'][['X', 'Y', 'Z']].drop_duplicates().values
        # Chromite lense intersection midpoints
        lenses = bh_df[bh_df['formation'] == 'chromite lense'][['X', 'Y', 'Z']].values
        # Fault point
        fault_pt = bh_df[bh_df['formation'] == 'Main_Fault'][['X', 'Y', 'Z']].values
        print(f"\nBorehole collars (host_rock): {len(collars)} points")
        print(f"  X range: {collars[:,0].min():.1f} – {collars[:,0].max():.1f}")
        print(f"  Y range: {collars[:,1].min():.1f} – {collars[:,1].max():.1f}")
        print(f"  Z range: {collars[:,2].min():.1f} – {collars[:,2].max():.1f}")
        print(f"Chromite lenses: {len(lenses)} intersections")
        print(f"Fault point: {fault_pt[0] if len(fault_pt) > 0 else 'none'}")
    else:
        collars, lenses, fault_pt = None, None, None
        print(f"\nBorehole data not found at: {borehole_path}")

    # --- Plot ---
    fig, axes = plt.subplots(2, 3, figsize=(24, 14))

    # 1) Raster image
    ax0 = axes[0, 0]
    im = ax0.imshow(
        masked,
        extent=(bounds.left, bounds.right, bounds.bottom, bounds.top),
        cmap='RdYlBu_r',
        aspect='equal',
    )
    cbar = fig.colorbar(im, ax=ax0, shrink=0.78)
    cbar.set_label('Magnetic TMI (nT)')
    ax0.set_title('Ternove – Merged/Masked B1B2 Magnetics (up 3m)')
    ax0.set_xlabel('Easting (EPSG:32634)')
    ax0.set_ylabel('Northing (EPSG:32634)')

    # 2) Raster sampled points (scatter)
    ax1 = axes[0, 1]
    sc = ax1.scatter(
        lon[idx], lat[idx], c=valid_vals[idx], cmap='RdYlBu_r',
        s=1, alpha=0.8,
    )
    cbar2 = fig.colorbar(sc, ax=ax1, shrink=0.78)
    cbar2.set_label('Magnetic TMI (nT)')
    ax1.set_title(f'Raster Pixel Values (n={n_sample:,} sampled)')
    ax1.set_xlabel('Easting (EPSG:32634)')
    ax1.set_ylabel('Northing (EPSG:32634)')
    ax1.set_aspect('equal')

    # 3) Raster histogram
    ax2 = axes[0, 2]
    ax2.hist(valid_vals, bins=100, color='steelblue', edgecolor='white', alpha=0.85)
    ax2.axvline(np.mean(valid_vals), color='red', linestyle='--',
                label=f'Mean = {np.mean(valid_vals):.1f}')
    ax2.axvline(np.median(valid_vals), color='orange', linestyle='--',
                label=f'Median = {np.median(valid_vals):.1f}')
    ax2.set_xlabel('Magnetics (nT)')
    ax2.set_ylabel('Pixel count')
    ax2.set_title('Raster Value Distribution')
    ax2.legend()

    # 4) Magnetic data on top of DEM with boreholes
    ax3 = axes[1, 0]
    ax3.imshow(
        dem_masked,
        extent=(dem_bounds.left, dem_bounds.right, dem_bounds.bottom, dem_bounds.top),
        cmap='terrain',
        aspect='equal',
        alpha=0.7,
    )
    sc3 = ax3.scatter(
        lon[idx], lat[idx], c=valid_vals[idx], cmap='RdYlBu_r',
        s=2, alpha=0.6,
    )
    # Overlay borehole locations
    if collars is not None:
        ax3.scatter(collars[:, 0], collars[:, 1], marker='v', c='black',
                     s=60, label='BH collar', zorder=5, edgecolors='white', linewidth=0.5)
    if lenses is not None and len(lenses) > 0:
        ax3.scatter(lenses[:, 0], lenses[:, 1], marker='o', c='gold',
                     s=40, label='chromite lense', zorder=5, edgecolors='black', linewidth=0.5)
    if fault_pt is not None and len(fault_pt) > 0:
        ax3.scatter(fault_pt[:, 0], fault_pt[:, 1], marker='*', c='red',
                     s=120, label='Main_Fault', zorder=6, edgecolors='black', linewidth=0.5)
    ax3.legend(fontsize=7, loc='lower left')
    cbar3 = fig.colorbar(sc3, ax=ax3, shrink=0.78)
    cbar3.set_label('Magnetic TMI (nT)')
    ax3.set_title('Magnetics over Soricom DEM + Boreholes')
    ax3.set_xlabel('Easting (EPSG:32634)')
    ax3.set_ylabel('Northing (EPSG:32634)')
    ax3.set_aspect('equal')

    # 5) Susceptibility histogram by lithology
    ax4 = axes[1, 1]
    if sus_df is not None and len(sus_df) > 0:
        litho_colors = {
            'chromite': '#8B0000',
            'chromite+serpentinite': '#CD5C5C',
            'serpentinite+chromite': '#CD5C5C',
            'harzburgite': '#2E8B57',
            'pyroxenite': '#4682B4',
            'serpentinite': '#DAA520',
            'granite?': '#A9A9A9',
        }
        sus_df['color'] = sus_df['lithology'].map(litho_colors).fillna('#999999')
        sus_df['label'] = sus_df['lithology'].replace({
            'chromite+serpentinite': 'chromite+serp',
            'serpentinite+chromite': 'serp+chromite',
            '': 'unknown',
        })

        for litho, grp in sus_df.groupby('label'):
            color = litho_colors.get(
                sus_df[sus_df['label'] == litho]['lithology'].iloc[0], '#999999'
            ) if len(grp) > 0 else '#999999'
            ax4.hist(grp['sus_SI'], bins=15, alpha=0.65, color=color,
                     label=f'{litho} (n={len(grp)})')

        ax4.set_xlabel('Susceptibility (SI)')
        ax4.set_ylabel('Count')
        ax4.set_title('Lab Susceptibility by Lithology')
        ax4.legend(fontsize=8)
    else:
        ax4.text(0.5, 0.5, 'No susceptibility data', transform=ax4.transAxes, ha='center')
        ax4.set_title('Lab Susceptibility – Not Available')

    # 6) Point data scatter (if available)
    ax5 = axes[1, 2]
    if points is not None:
        n_pt_sample = min(5000, len(points))
        pt_idx = rng.choice(len(points), size=n_pt_sample, replace=False)
        pts_sub = points.iloc[pt_idx]

        sc5 = ax5.scatter(
            pts_sub.geometry.x, pts_sub.geometry.y,
            c=pts_sub['MAG'], cmap='RdYlBu_r',
            s=1, alpha=0.8,
        )
        cbar5 = fig.colorbar(sc5, ax=ax5, shrink=0.78)
        cbar5.set_label('Magnetic (nT)')
        ax5.set_title(f'Cleaned Point Data (n={n_pt_sample:,} sampled / {len(points):,} total)')
    else:
        ax5.text(0.5, 0.5, 'No cleaned point data', transform=ax5.transAxes, ha='center')
        ax5.set_title('Cleaned Point Data – Not Available')
    ax5.set_xlabel('Easting (EPSG:32634)')
    ax5.set_ylabel('Northing (EPSG:32634)')
    ax5.set_aspect('equal')

    fig.tight_layout()
    plt.show()

    # Basic assertions
    assert data is not None
    assert valid_mask.sum() > 0
    if points is not None:
        assert len(points) > 0
    if sus_df is not None:
        assert len(sus_df) > 0