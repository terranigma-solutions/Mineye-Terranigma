"""Test reading and visualizing raw magnetic raster data for the Ternove model."""

import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
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

    # --- Read raster ---
    with rasterio.open(raster_path) as src:
        crs = src.crs
        bounds = src.bounds
        transform = src.transform
        data = src.read(1)

    # Mask nodata
    valid_mask = data != -999999
    masked = np.where(valid_mask, data, np.nan)

    # Extract all pixel coordinates (lon/lat) and values
    rows, cols = np.where(valid_mask)
    lon, lat = rasterio.transform.xy(transform, rows, cols)
    lon = np.array(lon)
    lat = np.array(lat)
    valid_vals = data[valid_mask]

    print(f"\nRaster shape: {data.shape}")
    print(f"CRS: {crs}")
    print(f"Bounds: {bounds}")
    print(f"Valid pixels: {valid_mask.sum():,} / {data.size:,}")
    print(f"Magnetic range: {np.nanmin(masked):.2f} – {np.nanmax(masked):.2f} nT")
    print(f"Mean: {np.nanmean(masked):.2f} nT, Std: {np.nanstd(masked):.2f} nT")

    # --- Read cleaned magnetic point data ---
    cleaned_path = os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson')

    if os.path.exists(cleaned_path):
        points = gpd.read_file(cleaned_path)
        print(f"\nCleaned point data: {len(points)} points")
        print(f"Point data CRS: {points.crs}")
        print(f"Point MAG range: {points['MAG'].min():.2f} – {points['MAG'].max():.2f} nT")

        # Reproject points to raster CRS
        if points.crs != crs:
            points = points.to_crs(crs)
    else:
        points = None
        print(f"\nCleaned point data not found at: {cleaned_path}")

    # --- Subsample raster points for scatter plot (too many to plot all) ---
    n_sample = min(5000, len(lon))
    rng = np.random.default_rng(42)
    idx = rng.choice(len(lon), size=n_sample, replace=False)

    # --- Plot ---
    fig, axes = plt.subplots(2, 2, figsize=(18, 14))

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
    ax0.set_xlabel('Longitude')
    ax0.set_ylabel('Latitude')

    # 2) Raster sampled points (scatter)
    ax1 = axes[0, 1]
    sc = ax1.scatter(
        lon[idx], lat[idx], c=valid_vals[idx], cmap='RdYlBu_r',
        s=1, alpha=0.8,
    )
    cbar2 = fig.colorbar(sc, ax=ax1, shrink=0.78)
    cbar2.set_label('Magnetic TMI (nT)')
    ax1.set_title(f'Raster Pixel Values (n={n_sample:,} sampled)')
    ax1.set_xlabel('Longitude')
    ax1.set_ylabel('Latitude')
    ax1.set_aspect('equal')

    # 3) Raster histogram
    ax2 = axes[1, 0]
    ax2.hist(valid_vals, bins=100, color='steelblue', edgecolor='white', alpha=0.85)
    ax2.axvline(np.mean(valid_vals), color='red', linestyle='--', label=f'Mean = {np.mean(valid_vals):.1f}')
    ax2.axvline(np.median(valid_vals), color='orange', linestyle='--', label=f'Median = {np.median(valid_vals):.1f}')
    ax2.set_xlabel('Magnetics (nT)')
    ax2.set_ylabel('Pixel count')
    ax2.set_title('Raster Value Distribution')
    ax2.legend()

    # 4) Point data scatter (if available)
    ax3 = axes[1, 1]
    if points is not None:
        # Subsample for display
        n_pt_sample = min(5000, len(points))
        pt_idx = rng.choice(len(points), size=n_pt_sample, replace=False)
        pts_sub = points.iloc[pt_idx]

        sc3 = ax3.scatter(
            pts_sub.geometry.x, pts_sub.geometry.y,
            c=pts_sub['MAG'], cmap='RdYlBu_r',
            s=1, alpha=0.8,
        )
        cbar3 = fig.colorbar(sc3, ax=ax3, shrink=0.78)
        cbar3.set_label('Magnetic (nT)')
        ax3.set_title(f'Cleaned Point Data (n={n_pt_sample:,} sampled / {len(points):,} total)')
    else:
        ax3.text(0.5, 0.5, 'No cleaned point data', transform=ax3.transAxes, ha='center')
        ax3.set_title('Cleaned Point Data – Not Available')
    ax3.set_xlabel('Longitude')
    ax3.set_ylabel('Latitude')
    ax3.set_aspect('equal')

    fig.tight_layout()
    plt.show()

    # Basic assertions
    assert data is not None
    assert valid_mask.sum() > 0
    if points is not None:
        assert len(points) > 0