"""Extract magnetic point observations from the Soricom raster for inversion.

Sub-samples the merged B1B2 magnetic raster to create a set of point
observations within the Soricom model extent, with Z from the DEM.

Strategy: sample on a coarse regular grid (step=5 pixels = 5m) within
the model extent, then randomly subset to n_points for the inversion.
"""

import os

import numpy as np
import rasterio
from matplotlib import pyplot as plt
from rasterio.windows import from_bounds

from mineye.config.example_parameters import SoricomModelConfig


def _extract_magnetic_points_from_raster(
    raster_path,
    extent,
    dem_path=None,
    step=5,
    nodata=-999999,
):
    """Extract magnetic point observations from a raster within a given extent.

    Parameters
    ----------
    raster_path : str
        Path to the magnetic raster (.ers or .tif).
    extent : list
        [min_x, max_x, min_y, max_y, min_z, max_z].
    dem_path : str, optional
        Path to DEM raster for Z elevation sampling.
    step : int
        Grid sampling step in pixels (1m = 1 pixel for this raster).
    nodata : float
        Nodata value in the magnetic raster.

    Returns
    -------
    xyz : np.ndarray
        Columns: X, Y, Z (m, EPSG:32634).
    mag_values : np.ndarray
        Magnetic TMI values (nT) at each point.
    """
    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        data = src.read(1, window=window)
        transform = src.window_transform(window)

    rows, cols = data.shape
    row_indices = np.arange(0, rows, step)
    col_indices = np.arange(0, cols, step)

    ii, jj = np.meshgrid(row_indices, col_indices, indexing='ij')
    ii_flat = ii.flatten()
    jj_flat = jj.flatten()

    mag_values = data[ii_flat, jj_flat]

    # Filter out nodata
    valid_mask = (mag_values != nodata) & ~np.isnan(mag_values)
    ii_valid = ii_flat[valid_mask]
    jj_valid = jj_flat[valid_mask]
    mag_values = mag_values[valid_mask]

    # Get XY coordinates (raster stores UTM coordinates despite declaring EPSG:4326)
    xs, ys = rasterio.transform.xy(transform, ii_valid.tolist(), jj_valid.tolist())
    xs = np.array(xs)
    ys = np.array(ys)

    # Get Z from DEM if available, otherwise use model top
    if dem_path and os.path.exists(dem_path):
        with rasterio.open(dem_path) as topo_src:
            coords = list(zip(xs, ys))
            zs = np.array([val[0] for val in topo_src.sample(coords)])
        # Mask nodata Z values (DEM uses 0 for nodata)
        z_valid = zs > 0
        if not np.all(z_valid):
            xs = xs[z_valid]
            ys = ys[z_valid]
            zs = zs[z_valid]
            mag_values = mag_values[z_valid]
    else:
        zs = np.full_like(xs, extent[5])

    xyz = np.column_stack((xs, ys, zs))
    return xyz, mag_values


def test_extract_magnetic_points(data_dir):
    """Extract magnetic point observations from the Soricom raster.

    Extracts regularly-spaced points from the merged B1B2 magnetic raster
    within the Soricom model extent, with Z elevation from the DEM.
    Plots the extracted points for visual verification.
    """
    raster_path = os.path.join(
        data_dir, 'Geophysical_Raw_Data',
        'merged_masked_B1B2_mean_median_201_up3m.ers',
    )
    dem_path = os.path.join(data_dir, 'Topographic_Data', 'soricom_DEM10m.tif')

    assert os.path.exists(raster_path), f"Raster not found: {raster_path}"
    assert os.path.exists(dem_path), f"DEM not found: {dem_path}"

    extent = SoricomModelConfig.EXTENT

    # Extract points at coarse resolution (step=5m)
    xyz, mag_values = _extract_magnetic_points_from_raster(
        raster_path=raster_path,
        extent=extent,
        dem_path=dem_path,
        step=5,
    )

    print(f"\nExtracted {len(xyz)} valid magnetic points within model extent")
    print(f"  XY extent: X={xyz[:, 0].min():.0f}–{xyz[:, 0].max():.0f} m")
    print(f"             Y={xyz[:, 1].min():.0f}–{xyz[:, 1].max():.0f} m")
    print(f"  Z  range: {xyz[:, 2].min():.1f}–{xyz[:, 2].max():.1f} m")
    print(f"  MAG range: {mag_values.min():.1f}–{mag_values.max():.1f} nT")
    print(f"  MAG mean:  {mag_values.mean():.1f} nT, std: {mag_values.std():.1f} nT")

    # Plot extracted points over raster
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    # Read full raster for background
    with rasterio.open(raster_path) as src:
        full_data = src.read(1)
        full_valid = full_data != -999999
        full_masked = np.where(full_valid, full_data, np.nan)

    ax0 = axes[0]
    ax0.imshow(
        full_masked,
        extent=(
            src.bounds.left, src.bounds.right,
            src.bounds.bottom, src.bounds.top,
        ),
        cmap='RdYlBu_r',
        aspect='equal',
        alpha=0.6,
    )
    sc = ax0.scatter(
        xyz[:, 0], xyz[:, 1], c=mag_values, cmap='RdYlBu_r',
        s=8, alpha=0.9, edgecolors='black', linewidth=0.3,
    )
    # Model extent rectangle
    ax0.plot(
        [extent[0], extent[1], extent[1], extent[0], extent[0]],
        [extent[2], extent[2], extent[3], extent[3], extent[2]],
        'r--', linewidth=1.5, label='Model extent',
    )
    cbar = fig.colorbar(sc, ax=ax0, shrink=0.75)
    cbar.set_label('Magnetic TMI (nT)')
    ax0.set_title(f'Extracted Points (n={len(xyz)}) over Raster', fontweight='bold')
    ax0.set_xlabel('Easting (EPSG:32634)')
    ax0.set_ylabel('Northing (EPSG:32634)')
    ax0.legend(fontsize=8)

    # Histogram
    ax1 = axes[1]
    ax1.hist(mag_values, bins=40, color='steelblue', edgecolor='white', alpha=0.85)
    ax1.axvline(mag_values.mean(), color='red', linestyle='--',
                label=f'Mean = {mag_values.mean():.1f} nT')
    ax1.set_xlabel('Magnetic TMI (nT)')
    ax1.set_ylabel('Count')
    ax1.set_title('Extracted Point Values Distribution')
    ax1.legend()

    fig.tight_layout()
    plt.show()

    assert len(xyz) > 0, "No magnetic points extracted"
    assert mag_values is not None

    # Save for later use by inversion test
    out_dir = os.path.join(os.path.dirname(__file__))
    np.save(os.path.join(out_dir, 'soricom_magnetic_xyz.npy'), xyz)
    np.save(os.path.join(out_dir, 'soricom_magnetic_values.npy'), mag_values)
    print(f"\nSaved point data to soricom_magnetic_xyz.npy and soricom_magnetic_values.npy")
