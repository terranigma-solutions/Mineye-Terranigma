"""Extract magnetic point observations from the Soricom raster for inversion.

Sub-samples the merged B1B2 magnetic raster to create a set of point
observations within the Soricom model extent, with Z from the DEM.

Strategies:
- Regular grid (baseline): coarse grid at step=5m
- Spatially reduced (gradient-adaptive): dense at magnetic contacts, sparse elsewhere
- Gradient peak: points at local maxima of gradient magnitude (most informative)

All strategies aim to minimise spatial correlation while preserving information.
"""

import os

import numpy as np
import rasterio
from matplotlib import pyplot as plt
from rasterio.windows import from_bounds
from scipy import ndimage
from skimage.feature import peak_local_max

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


def _extract_magnetic_points_spatially_reduced(
    raster_path,
    extent,
    dem_path=None,
    high_grad_step=10,
    low_grad_step=50,
    grad_threshold_percentile=80,
    nodata=-999999,
    margin=0.02,
):
    """Extract magnetic points using gradient-adaptive sampling.

    Computes local gradient magnitude (Sobel) as a proxy for magnetic
    contacts / structural boundaries. Samples densely in high-gradient
    areas and sparsely in homogeneous ones, minimising spatial
    correlation while preserving edge information.

    Parameters
    ----------
    raster_path : str
        Path to the magnetic raster.
    extent : list
        [min_x, max_x, min_y, max_y, min_z, max_z].
    dem_path : str, optional
        Path to DEM for Z elevation.
    high_grad_step : int
        Grid step (pixels) in high-gradient areas.
    low_grad_step : int
        Grid step (pixels) in low-gradient areas.
    grad_threshold_percentile : float
        Percentile of gradient magnitude used to separate high/low zones.
    nodata : float
        Nodata value in raster.
    margin : float
        Fraction of rows/cols to crop from edges.

    Returns
    -------
    xyz : np.ndarray
        Columns: X, Y, Z.
    mag_values : np.ndarray
        Magnetic values at each point.
    """
    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        data = src.read(1, window=window)
        transform = src.window_transform(window)

    rows, cols = data.shape

    # Crop margin
    row_margin = int(rows * margin)
    col_margin = int(cols * margin)
    crop_mask = np.zeros_like(data, dtype=bool)
    crop_mask[row_margin:rows - row_margin, col_margin:cols - col_margin] = True

    # Valid-data mask (nodata excluded)
    valid_mask = (data != nodata) & ~np.isnan(data) & crop_mask

    # Compute gradient magnitude via Sobel
    data_filled = np.nan_to_num(data, nan=0.0)
    data_filled[~valid_mask] = 0.0
    grad_x = ndimage.sobel(data_filled, axis=1)
    grad_y = ndimage.sobel(data_filled, axis=0)
    grad_mag = np.sqrt(grad_x ** 2 + grad_y ** 2)

    # Threshold to split high / low gradient
    valid_grad = grad_mag[valid_mask]
    if len(valid_grad) == 0:
        return np.zeros((0, 3)), np.zeros(0)
    threshold = np.percentile(valid_grad, grad_threshold_percentile)

    high_grad_mask = (grad_mag >= threshold) & valid_mask
    low_grad_mask = (grad_mag < threshold) & valid_mask

    # Dense grid in high-gradient areas
    ii_h, jj_h = np.meshgrid(
        np.arange(0, rows, high_grad_step),
        np.arange(0, cols, high_grad_step),
        indexing='ij',
    )
    ii_h, jj_h = ii_h.flatten(), jj_h.flatten()
    keep_h = high_grad_mask[ii_h, jj_h]
    ii_h, jj_h = ii_h[keep_h], jj_h[keep_h]

    # Sparse grid in low-gradient areas
    ii_l, jj_l = np.meshgrid(
        np.arange(0, rows, low_grad_step),
        np.arange(0, cols, low_grad_step),
        indexing='ij',
    )
    ii_l, jj_l = ii_l.flatten(), jj_l.flatten()
    keep_l = low_grad_mask[ii_l, jj_l]
    ii_l, jj_l = ii_l[keep_l], jj_l[keep_l]

    ii = np.concatenate([ii_h, ii_l])
    jj = np.concatenate([jj_h, jj_l])

    if len(ii) == 0:
        return np.zeros((0, 3)), np.zeros(0)

    mag_values = data[ii, jj]

    # XY coordinates
    xs, ys = rasterio.transform.xy(transform, ii.tolist(), jj.tolist())
    xs = np.array(xs)
    ys = np.array(ys)

    # Z from DEM
    zs = _sample_dem(dem_path, xs, ys, extent)

    xyz = np.column_stack((xs, ys, zs))
    return xyz, mag_values


def _extract_magnetic_points_gradient_peaks(
    raster_path,
    extent,
    dem_path=None,
    min_distance=50,
    n_strongest=None,
    nodata=-999999,
    margin=0.02,
):
    """Extract magnetic points at local maxima of gradient magnitude.

    Uses peak_local_max on the Sobel gradient magnitude to identify
    the most informative locations — magnetic contacts, anomaly edges,
    and structural boundaries.  The ``min_distance`` parameter
    directly controls the minimum spatial separation between selected
    points, guaranteeing low spatial correlation.

    Parameters
    ----------
    raster_path : str
        Path to the magnetic raster.
    extent : list
        [min_x, max_x, min_y, max_y, min_z, max_z].
    dem_path : str, optional
        Path to DEM for Z elevation.
    min_distance : int
        Minimum pixel distance between selected peaks.
    n_strongest : int, optional
        Keep only the *n* strongest peaks (by gradient magnitude).
    nodata : float
        Nodata value in raster.
    margin : float
        Fraction of rows/cols to crop from edges.

    Returns
    -------
    xyz : np.ndarray
        Columns: X, Y, Z.
    mag_values : np.ndarray
        Magnetic values at each point.
    """
    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        data = src.read(1, window=window)
        transform = src.window_transform(window)

    rows, cols = data.shape

    # Margin crop mask (applied after peak finding)
    row_margin = int(rows * margin)
    col_margin = int(cols * margin)

    # Valid-data mask
    valid_mask = (data != nodata) & ~np.isnan(data)

    # Compute gradient magnitude
    data_filled = np.nan_to_num(data, nan=0.0)
    data_filled[~valid_mask] = 0.0
    grad_x = ndimage.sobel(data_filled, axis=1)
    grad_y = ndimage.sobel(data_filled, axis=0)
    grad_mag = np.sqrt(grad_x ** 2 + grad_y ** 2)

    # Mask out nodata / marginal regions for peak finding
    peak_mask = valid_mask.copy()
    if margin > 0:
        peak_mask[:row_margin, :] = False
        peak_mask[rows - row_margin:, :] = False
        peak_mask[:, :col_margin] = False
        peak_mask[:, cols - col_margin:] = False

    # Remove pixels with zero gradient (perfectly homogeneous — uninformative)
    peak_mask &= (grad_mag > 0)

    peaks = peak_local_max(
        grad_mag,
        min_distance=min_distance,
        labels=peak_mask,
        exclude_border=False,
    )

    if len(peaks) == 0:
        return np.zeros((0, 3)), np.zeros(0)

    # Optionally keep only the strongest peaks
    if n_strongest is not None and len(peaks) > n_strongest:
        peak_strength = grad_mag[tuple(peaks.T)]
        idx = np.argsort(peak_strength)[::-1][:n_strongest]
        peaks = peaks[idx]

    ii, jj = peaks[:, 0], peaks[:, 1]

    mag_values = data[ii, jj]

    # XY coordinates
    xs, ys = rasterio.transform.xy(transform, ii.tolist(), jj.tolist())
    xs = np.array(xs)
    ys = np.array(ys)

    # Z from DEM
    zs = _sample_dem(dem_path, xs, ys, extent)

    xyz = np.column_stack((xs, ys, zs))
    return xyz, mag_values


def _sample_dem(dem_path, xs, ys, extent):
    """Sample Z from a DEM raster, falling back to model top."""
    if dem_path and os.path.exists(dem_path):
        with rasterio.open(dem_path) as topo_src:
            coords = list(zip(xs, ys))
            zs = np.array([val[0] for val in topo_src.sample(coords)])
        z_valid = zs > 0
        if not np.all(z_valid):
            zs[~z_valid] = extent[5]
    else:
        zs = np.full_like(xs, extent[5], dtype=float)
    return zs


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def test_extract_magnetic_points(data_dir):
    """Extract magnetic point observations from the Soricom raster.

    Extracts regularly-spaced points from the merged B1B2 magnetic raster
    within the Soricom model extent, with Z elevation from the DEM.
    Plots the extracted points for visual verification.
    """
    raster_path = os.path.join(
        data_dir, 'Geophysical_Raw_Data',
        'magnetics_model2.ers',
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


def test_spatial_correlation_reduction(data_dir):
    """Compare gradient-adaptive sampling against standard regular grid.

    The adaptive strategy keeps more points where the magnetic field
    changes rapidly (magnetic contacts / structural boundaries) and
    fewer in homogeneous areas, reducing total point count while
    preserving the information content at edges.
    """
    raster_path = os.path.join(
        data_dir, 'Geophysical_Raw_Data',
        'magnetics_model2.ers',
    )
    dem_path = os.path.join(data_dir, 'Topographic_Data', 'soricom_DEM10m.tif')

    assert os.path.exists(raster_path), f"Raster not found: {raster_path}"
    assert os.path.exists(dem_path), f"DEM not found: {dem_path}"

    extent = SoricomModelConfig.EXTENT

    # Standard regular grid (baseline)
    xyz_std, mag_std = _extract_magnetic_points_from_raster(
        raster_path, extent, dem_path=dem_path, step=5,
    )

    # Gradient-adaptive sampling
    xyz_adapt, mag_adapt = _extract_magnetic_points_spatially_reduced(
        raster_path, extent, dem_path=dem_path,
        high_grad_step=10, low_grad_step=50,
        grad_threshold_percentile=80,
    )

    print(f"\nComparison — Standard vs Gradient-Adaptive:")
    print(f"  Standard grid (step=5):          {len(xyz_std):>6} points")
    print(f"  Gradient-adaptive:              {len(xyz_adapt):>6} points")
    print(f"  Reduction factor:               {len(xyz_std) / max(len(xyz_adapt), 1):.2f}x")
    print(f"  Value range — standard:          {mag_std.min():.1f} – {mag_std.max():.1f} nT")
    print(f"  Value range — adaptive:         {mag_adapt.min():.1f} – {mag_adapt.max():.1f} nT")

    # Plot comparison
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    with rasterio.open(raster_path) as src:
        full_data = src.read(1)
        full_valid = full_data != -999999
        full_masked = np.where(full_valid, full_data, np.nan)
        raster_extent = (src.bounds.left, src.bounds.right,
                         src.bounds.bottom, src.bounds.top)

    for ax, xyz, mag, title in zip(
        axes,
        [xyz_std, xyz_adapt],
        [mag_std, mag_adapt],
        [f'Standard Grid (n={len(xyz_std)})',
         f'Gradient-Adaptive (n={len(xyz_adapt)})'],
    ):
        ax.imshow(full_masked, extent=raster_extent, cmap='RdYlBu_r',
                  aspect='equal', alpha=0.4)
        sc = ax.scatter(xyz[:, 0], xyz[:, 1], c=mag, cmap='RdYlBu_r',
                        s=6, alpha=0.85, edgecolors='black', linewidth=0.3)
        ax.plot(
            [extent[0], extent[1], extent[1], extent[0], extent[0]],
            [extent[2], extent[2], extent[3], extent[3], extent[2]],
            'r--', linewidth=1,
        )
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        fig.colorbar(sc, ax=ax, shrink=0.7)

    # Histogram of gradient magnitudes in the selected regions
    ax2 = axes[2]
    # Read raster and compute gradient for histogram
    with rasterio.open(raster_path) as src:
        window = from_bounds(extent[0], extent[2], extent[1], extent[3], src.transform)
        data_full = src.read(1, window=window)
    valid = (data_full != -999999) & ~np.isnan(data_full)
    data_filled = np.nan_to_num(data_full, nan=0.0)
    data_filled[~valid] = 0.0
    gx = ndimage.sobel(data_filled, axis=1)
    gy = ndimage.sobel(data_filled, axis=0)
    gm = np.sqrt(gx ** 2 + gy ** 2)
    ax2.hist(gm[valid].ravel(), bins=50, color='gray', alpha=0.6, label='All pixels')
    ax2.axvline(np.percentile(gm[valid], 80), color='red', linestyle='--',
                label='80th percentile (threshold)')
    ax2.set_xlabel('Gradient magnitude')
    ax2.set_ylabel('Count')
    ax2.set_title('Gradient Distribution & Sampling Threshold')
    ax2.legend()

    fig.tight_layout()
    plt.show()

    # Assertions
    assert len(xyz_adapt) > 0, "Should extract points"
    assert len(xyz_adapt) < len(xyz_std), \
        f"Adaptive ({len(xyz_adapt)}) should reduce point count vs standard ({len(xyz_std)})"
    assert abs(mag_adapt.mean() - mag_std.mean()) / abs(mag_std.mean()) < 0.15, \
        "Mean should be preserved within 15% after adaptive sampling"

    # Save
    out_dir = os.path.join(os.path.dirname(__file__))
    np.save(os.path.join(out_dir, 'soricom_magnetic_xyz_adaptive.npy'), xyz_adapt)
    np.save(os.path.join(out_dir, 'soricom_magnetic_values_adaptive.npy'), mag_adapt)
    print(f"\nSaved gradient-adaptive points to soricom_magnetic_*_adaptive.npy")


def test_gradient_peak_extraction(data_dir):
    """Extract minimal magnetic points at local gradient maxima.

    Uses ``peak_local_max`` on the Sobel gradient magnitude with a
    ``min_distance`` constraint to guarantee a minimum spatial
    separation between selected points, directly controlling spatial
    correlation.
    """
    raster_path = os.path.join(
        data_dir, 'Geophysical_Raw_Data',
        'magnetics_model2.ers',
    )
    dem_path = os.path.join(data_dir, 'Topographic_Data', 'soricom_DEM10m.tif')

    assert os.path.exists(raster_path), f"Raster not found: {raster_path}"
    assert os.path.exists(dem_path), f"DEM not found: {dem_path}"

    extent = SoricomModelConfig.EXTENT

    # 1. Standard grid (reference)
    xyz_std, mag_std = _extract_magnetic_points_from_raster(
        raster_path, extent, dem_path=dem_path, step=5,
    )

    # 2. Gradient peaks with moderate separation
    xyz_peaks, mag_peaks = _extract_magnetic_points_gradient_peaks(
        raster_path, extent, dem_path=dem_path,
        min_distance=50, n_strongest=None,
    )

    # 3. Gradient peaks with stronger separation (fewer points)
    xyz_peaks_sparse, mag_peaks_sparse = _extract_magnetic_points_gradient_peaks(
        raster_path, extent, dem_path=dem_path,
        min_distance=100, n_strongest=None,
    )

    print(f"\nComparison — Standard vs Gradient Peaks:")
    print(f"  Standard grid (step=5):         {len(xyz_std):>6} points")
    print(f"  Gradient peaks (d=50):           {len(xyz_peaks):>6} points "
          f"({len(xyz_std) / max(len(xyz_peaks), 1):.1f}x reduction)")
    print(f"  Gradient peaks (d=100):          {len(xyz_peaks_sparse):>6} points "
          f"({len(xyz_std) / max(len(xyz_peaks_sparse), 1):.1f}x reduction)")

    print(f"  Value range — standard:         {mag_std.min():.1f} – {mag_std.max():.1f} nT")
    print(f"  Value range — peaks d=50:       {mag_peaks.min():.1f} – {mag_peaks.max():.1f} nT")
    print(f"  Value range — peaks d=100:      {mag_peaks_sparse.min():.1f} – {mag_peaks_sparse.max():.1f} nT")

    # Plot
    fig, axes = plt.subplots(1, 3, figsize=(22, 7))

    with rasterio.open(raster_path) as src:
        full_data = src.read(1)
        full_valid = full_data != -999999
        full_masked = np.where(full_valid, full_data, np.nan)
        raster_extent = (src.bounds.left, src.bounds.right,
                         src.bounds.bottom, src.bounds.top)

    datasets = [
        (xyz_std, mag_std, f'Standard Grid (n={len(xyz_std)})'),
        (xyz_peaks, mag_peaks, f'Gradient Peaks d=50 (n={len(xyz_peaks)})'),
        (xyz_peaks_sparse, mag_peaks_sparse,
         f'Gradient Peaks d=100 (n={len(xyz_peaks_sparse)})'),
    ]
    for ax, (xyz, mag, title) in zip(axes, datasets):
        ax.imshow(full_masked, extent=raster_extent, cmap='RdYlBu_r',
                  aspect='equal', alpha=0.4)
        sc = ax.scatter(xyz[:, 0], xyz[:, 1], c=mag, cmap='RdYlBu_r',
                        s=12, alpha=0.9, edgecolors='black', linewidth=0.4)
        ax.plot(
            [extent[0], extent[1], extent[1], extent[0], extent[0]],
            [extent[2], extent[2], extent[3], extent[3], extent[2]],
            'r--', linewidth=1,
        )
        ax.set_title(title, fontweight='bold')
        ax.set_xlabel('Easting (m)')
        ax.set_ylabel('Northing (m)')
        fig.colorbar(sc, ax=ax, shrink=0.7)

    fig.tight_layout()
    plt.show()

    # Assertions
    assert len(xyz_peaks) > 0, "Should extract peak points"
    assert len(xyz_peaks) < len(xyz_std), \
        f"Peaks ({len(xyz_peaks)}) should be fewer than standard ({len(xyz_std)})"
    assert len(xyz_peaks_sparse) < len(xyz_peaks), \
        "Larger min_distance should yield fewer points"
    assert abs(mag_peaks.mean() - mag_std.mean()) / abs(mag_std.mean()) < 0.20, \
        "Peak-based mean should be within 20% of standard mean"

    # Save
    out_dir = os.path.join(os.path.dirname(__file__))
    np.save(os.path.join(out_dir, 'soricom_magnetic_xyz_peaks.npy'), xyz_peaks)
    np.save(os.path.join(out_dir, 'soricom_magnetic_values_peaks.npy'), mag_peaks)
    print(f"\nSaved gradient peak points to soricom_magnetic_*_peaks.npy")