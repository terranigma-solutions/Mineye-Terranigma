import os

import numpy as np
import pytest
import rasterio
from matplotlib import pyplot as plt
from rasterio.windows import from_bounds

if rasterio is None:
    pytest.skip("rasterio is required for reading EnMap files")


def test_read_EnMap(base_dir, model_extent, simple_geo_model):
    """Test reading EnMap files and plotting them.

    This test reads the EnMap segmentation results (which are at the same
    resolution as the GemPy model) and creates visualizations.
    """

    # Path to EnMap segmentation results
    enmap_dir = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap')

    # Available EnMap files (segmentation results)
    enmap_files = {
            'result_n4' : os.path.join(enmap_dir, 'EPSG3857_EnMap_result_n4_betajump0.1.tif'),
            'result_n6' : os.path.join(enmap_dir, 'EnMap_result_n6_betajump0.1.tif'),
            'result_n8' : os.path.join(enmap_dir, 'EnMap_result_n8_betajump0.1.tif'),
            'entropy_n4': os.path.join(enmap_dir, 'EnMap_entropy_n4_betajump0.1.tif'),
    }

    # Check which files exist
    available_files = {k: v for k, v in enmap_files.items() if os.path.exists(v)}

    if not available_files:
        print(f"\n⚠️  No EnMap files found in {enmap_dir}")
        print("    Please add EnMap data to proceed with this test.")
        pytest.skip("No EnMap data files available")

    print(f"\n✓ Found {len(available_files)} EnMap file(s)")

    # Read and plot the first available file
    file_key = list(available_files.keys())[0]
    file_path = available_files[file_key]

    print(f"\n📖 Reading: {os.path.basename(file_path)}")

    with rasterio.open(file_path) as src:
        # Read the data
        data = src.read(1)  # Read first band
        transform = src.transform
        bounds = src.bounds
        crs = src.crs

        print(f"   Shape: {data.shape}")
        print(f"   Bounds: {bounds}")
        print(f"   CRS: {crs}")
        print(f"   Resolution: {transform.a:.2f} x {abs(transform.e):.2f} m")
        print(f"   Data range: [{np.nanmin(data):.2f}, {np.nanmax(data):.2f}]")

        # Create coordinates for plotting
        rows, cols = data.shape
        x = np.linspace(bounds.left, bounds.right, cols)
        y = np.linspace(bounds.bottom, bounds.top, rows)
        X, Y = np.meshgrid(x, y)

        # Check if coordinates align with model extent
        print(f"\n   Model extent: {model_extent[:4]}")
        print(f"   EnMap bounds: [{bounds.left:.0f}, {bounds.right:.0f}, {bounds.bottom:.0f}, {bounds.top:.0f}]")

        # Create visualization
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        # Plot 1: The EnMap data
        if 'result' in file_key:
            # Segmentation result - use discrete colormap
            im1 = axes[0].imshow(data, cmap='tab10', interpolation='nearest',
                                 extent=[bounds.left, bounds.right, bounds.bottom, bounds.top])
            axes[0].set_title(f'EnMap Segmentation Result ({file_key})', fontsize=12, fontweight='bold')
            cbar1 = plt.colorbar(im1, ax=axes[0])
            cbar1.set_label('Class ID', fontsize=10)
        else:
            # Entropy - use continuous colormap
            im1 = axes[0].imshow(data, cmap='viridis', interpolation='bilinear',
                                 extent=[bounds.left, bounds.right, bounds.bottom, bounds.top])
            axes[0].set_title(f'EnMap Entropy ({file_key})', fontsize=12, fontweight='bold')
            cbar1 = plt.colorbar(im1, ax=axes[0])
            cbar1.set_label('Entropy', fontsize=10)

        axes[0].set_xlabel('X (m)', fontsize=10)
        axes[0].set_ylabel('Y (m)', fontsize=10)
        axes[0].grid(True, alpha=0.3)
        axes[0].scatter(
            simple_geo_model.surface_points_copy.xyz[:, 0],
            simple_geo_model.surface_points_copy.xyz[:, 1], c='red', s=10,
            zorder=5, label='Model Surface Points', edgecolors='black', linewidth=0.5)

        # Plot 2: Histogram
        valid_data = data[~np.isnan(data)]
        axes[1].hist(valid_data.flatten(), bins=50, color='steelblue', alpha=0.7, edgecolor='black')
        axes[1].set_xlabel('Value', fontsize=10)
        axes[1].set_ylabel('Frequency', fontsize=10)
        axes[1].set_title('Data Distribution', fontsize=12, fontweight='bold')
        axes[1].grid(True, alpha=0.3)

        plt.tight_layout()
        plt.show()

        # Basic assertions
        assert data.shape[0] > 0 and data.shape[1] > 0, "Data should not be empty"
        assert not np.all(np.isnan(data)), "Data should not be all NaN"

        print(f"\n✓ Successfully read and plotted EnMap data")
        print(f"   Valid pixels: {np.sum(~np.isnan(data))} / {data.size} ({100 * np.sum(~np.isnan(data)) / data.size:.1f}%)")


def _extract_points_from_raster(raster_path, extent, step=10, z_value=None, topo_path=None, margin=0.0):
    """
    Private method to extract points from a raster within a given extent.
    """
    
    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        
        # Read the data within the window
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        
        # 0. Crop points near borders
        rows, cols = data.shape
        row_margin = int(rows * margin)
        col_margin = int(cols * margin)
        
        crop_mask = np.zeros_like(data, dtype=bool)
        crop_mask[row_margin:rows-row_margin, col_margin:cols-col_margin] = True
        
        rows, cols = data.shape
        row_indices = np.arange(0, rows, step)
        col_indices = np.arange(0, cols, step)
        
        ii, jj = np.meshgrid(row_indices, col_indices, indexing='ij')
        
        # Flatten indices
        ii_flat = ii.flatten()
        jj_flat = jj.flatten()

        # Filter by crop mask
        in_crop = crop_mask[ii_flat, jj_flat]
        ii_flat = ii_flat[in_crop]
        jj_flat = jj_flat[in_crop]
        
        # Get labels
        labels = data[ii_flat, jj_flat]
        
        # Filter out NaNs
        valid_mask = ~np.isnan(labels)
        
        # Apply specific label logic:
        # 1. Ignore label 1
        # 2. Combine label 3 and 0 (let's map 3 to 0)
        
        # Update valid_mask to ignore label 1
        valid_mask &= (labels != 1)
        
        ii_valid = ii_flat[valid_mask]
        jj_valid = jj_flat[valid_mask]
        labels_valid = labels[valid_mask]
        
        # Combine label 3 and 0 -> set all 3s to 0
        labels_valid[labels_valid == 3] = 0
        
        # Get xy coordinates
        xs, ys = rasterio.transform.xy(transform, ii_valid.tolist(), jj_valid.tolist())
        xs = np.array(xs)
        ys = np.array(ys)
        
        if topo_path:
            with rasterio.open(topo_path) as topo_src:
                coords = zip(xs, ys)
                zs = np.array([val[0] for val in topo_src.sample(coords)])
        else:
            if z_value is None:
                z_value = extent[5]
            zs = np.full_like(xs, z_value)
        
        xyz = np.column_stack((xs, ys, zs))
        
        return xyz, labels_valid, data, (left, right, bottom, top)


def _extract_points_spatially_reduced(raster_path, extent, step_boundary=2, step_inner=20, kernel_size=5, z_value=None, topo_path=None, margin=0.0):
    """
    Extract points using a spatial kernel to identify boundaries and reduce points in homogeneous areas.
    """
    from skimage.segmentation import find_boundaries

    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        
        # Read the data within the window at full resolution
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        
        # 0. Crop points near borders
        rows, cols = data.shape
        row_margin = int(rows * margin)
        col_margin = int(cols * margin)
        
        crop_mask = np.zeros_like(data, dtype=bool)
        crop_mask[row_margin:rows-row_margin, col_margin:cols-col_margin] = True
        
        # Clean data (handle NaNs and label mapping as in previous step)
        # We need to fill NaNs for boundary detection, but we'll mask them out later
        data_clean = data.copy()
        mask_nan = np.isnan(data)
        data_clean[mask_nan] = 255 # Temporary label for NaNs (uint8 safe)
        
        # Map labels (label 3 to 0, label 1 to be ignored)
        # Note: We do this on the full resolution data to find boundaries correctly
        data_mapped = data_clean.copy()
        data_mapped[data_mapped == 3] = 0
        
        # Identify boundaries using a spatial kernel (find_boundaries uses a neighborhood check)
        boundaries = find_boundaries(data_mapped, mode='thick')
        
        # Create a sampling mask
        # 1. Sample boundary points at high frequency
        boundary_mask = np.zeros_like(boundaries, dtype=bool)
        boundary_mask[::step_boundary, ::step_boundary] = True
        boundary_mask &= boundaries
        
        # 2. Sample inner points at low frequency
        inner_mask = np.zeros_like(boundaries, dtype=bool)
        inner_mask[::step_inner, ::step_inner] = True
        inner_mask &= ~boundaries
        
        # Combine masks
        sampling_mask = (boundary_mask | inner_mask) & crop_mask
        
        # Mask out NaNs and ignored labels (label 1)
        sampling_mask &= ~mask_nan
        sampling_mask &= (data_mapped != 1)
        
        # Get indices
        ii, jj = np.where(sampling_mask)
        labels_valid = data_mapped[ii, jj]
        
        # Get xy coordinates
        xs, ys = rasterio.transform.xy(transform, ii.tolist(), jj.tolist())
        xs = np.array(xs)
        ys = np.array(ys)
        
        if topo_path:
            with rasterio.open(topo_path) as topo_src:
                coords = zip(xs, ys)
                zs = np.array([val[0] for val in topo_src.sample(coords)])
        else:
            if z_value is None:
                z_value = extent[5]
            zs = np.full_like(xs, z_value)
        
        xyz = np.column_stack((xs, ys, zs))
        
        return xyz, labels_valid, data, (left, right, bottom, top)


def _extract_points_central_reduced(raster_path, extent, min_distance=25, z_value=None, topo_path=None, balance_patches=True, margin=0.05):
    """
    Extract points from the center of geological bodies using distance transform.
    Prioritizes points furthest from boundaries and ensures they are spatially separated.
    
    If balance_patches is True, it tries to limit the number of points from small patches
    and potentially increase the number of points from large patches to avoid over-representation.
    """
    from skimage.segmentation import find_boundaries
    from skimage.feature import peak_local_max
    from scipy import ndimage

    with rasterio.open(raster_path) as src:
        left, right, bottom, top = extent[0], extent[1], extent[2], extent[3]
        window = from_bounds(left, bottom, right, top, src.transform)
        
        # Read the data within the window at full resolution
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        
        # 0. Crop points near borders
        rows, cols = data.shape
        row_margin = int(rows * margin)
        col_margin = int(cols * margin)
        
        crop_mask = np.zeros_like(data, dtype=bool)
        crop_mask[row_margin:rows-row_margin, col_margin:cols-col_margin] = True
        
        # Clean data (handle NaNs and label mapping)
        data_mapped = data.copy()
        mask_nan = np.isnan(data)
        data_mapped[data_mapped == 3] = 0
        
        # 1. Identify boundaries
        # We need to fill NaNs for boundary detection
        data_temp = data_mapped.copy()
        data_temp[mask_nan] = 255
        boundaries = find_boundaries(data_temp, mode='thick')
        
        # 2. Distance transform: distance to nearest boundary or NaN
        # We want to be far from boundaries AND far from NaN areas (which are outside the domain)
        dist_mask = ~boundaries & ~mask_nan & crop_mask
        dist_transform = ndimage.distance_transform_edt(dist_mask)
        
        # 3. Extract peaks for each label
        unique_labels = np.unique(data_mapped)
        unique_labels = unique_labels[~np.isnan(unique_labels)]
        unique_labels = unique_labels[unique_labels != 1]
        
        all_ii = []
        all_jj = []
        all_labels = []

        # Area calculation for balancing
        if balance_patches:
            label_areas = {}
            for label_val in unique_labels:
                label_areas[label_val] = np.sum(data_mapped == label_val)
            
            max_area = max(label_areas.values()) if label_areas else 1
            min_area = min(label_areas.values()) if label_areas else 1
            
            print(f"\n   Area balancing enabled:")
            print(f"   Area range: {min_area} to {max_area} pixels")
        
        for label_val in unique_labels:
            mask = (data_mapped == label_val)
            
            # Dynamic min_distance based on area? 
            # Or just limit the number of points?
            
            current_min_dist = min_distance
            if balance_patches:
                # If area is small, we might want to increase min_distance to get fewer points
                # or just use num_peaks in peak_local_max
                area_ratio = label_areas[label_val] / max_area
                # If area is only 1% of the max area, we probably don't want many points.
                # However, we want at least one point if possible.
                
                # Heuristic: cap number of points proportional to area
                # total_points_budget = (data_mapped.size / (min_distance**2))
                # expected_points = total_points_budget * (label_areas[label_val] / data_mapped.size)
                
                # Another approach: adjust min_distance. 
                # Larger area -> smaller min_distance (more dense)
                # Smaller area -> larger min_distance (less dense)
                # This is counter-intuitive if we want to AVOID over-representing small patches.
                # Actually, small patches get over-represented because they might only fit ONE point,
                # but that one point represents a tiny area, while a large patch has many points,
                # but maybe fewer per unit area than the small patch.
                
                # Actually, the user says "By removing redundant points on big areas it seems that we are over representing small patches."
                # This means the big areas are being thinned out TOO MUCH compared to small patches.
                # So we should either:
                # 1. Thin out small patches more.
                # 2. Thin out big areas less.
                
                # Let's try to limit the number of peaks for small patches.
                pass

            # Use peak_local_max to find central points
            peaks = peak_local_max(
                dist_transform, 
                min_distance=current_min_dist, 
                labels=mask,
                exclude_border=False
            )
            
            if balance_patches:
                # Limit points for small patches
                # Say we want point density to be somewhat uniform.
                # area / num_peaks should be roughly constant.
                
                # We use the distance transform to estimate a "natural" number of peaks 
                # given the min_distance. But for very small patches, even 1 peak might 
                # be "over-representing" them relative to their area if we have thinned 
                # out big areas significantly.
                
                target_density = 1 / (min_distance**2) # points per pixel
                # Heuristic: Target number of points proportional to area.
                # We want a more uniform distribution.
                # Actually, if we use target_density = 1 / (min_distance**2), 
                # peak_local_max ALREADY tries to achieve this density.
                # To actually thinning it out MORE than peak_local_max, we should use a stricter density.
                target_num_peaks = max(1, int(label_areas[label_val] * target_density * 0.2)) 
                
                # Further correction: if a patch is very small, we might want to skip it 
                # or strictly limit it to 1 point if it's below a threshold.
                min_patch_area = (min_distance ** 2)
                if label_areas[label_val] < min_patch_area:
                     # Even if it's small, we keep at most 1 point if it has any valid pixel
                     target_num_peaks = 1 
                
                if len(peaks) > target_num_peaks:
                     print(f"      Label {label_val}: Area={label_areas[label_val]}, Peaks={len(peaks)} -> Reduced to {target_num_peaks}")
                     # Sort peaks by distance transform value (most central first)
                     peak_vals = dist_transform[tuple(peaks.T)]
                     idx = np.argsort(peak_vals)[::-1]
                     peaks = peaks[idx[:target_num_peaks]]

            if len(peaks) > 0:
                all_ii.extend(peaks[:, 0])
                all_jj.extend(peaks[:, 1])
                all_labels.extend([label_val] * len(peaks))
        
        ii = np.array(all_ii)
        jj = np.array(all_jj)
        labels_valid = np.array(all_labels)
        
        if len(ii) == 0:
            return np.zeros((0, 3)), np.zeros(0), data, (left, right, bottom, top)

        # Get xy coordinates
        xs, ys = rasterio.transform.xy(transform, ii.tolist(), jj.tolist())
        xs = np.array(xs)
        ys = np.array(ys)
        
        if topo_path:
            with rasterio.open(topo_path) as topo_src:
                coords = zip(xs, ys)
                zs = np.array([val[0] for val in topo_src.sample(coords)])
        else:
            if z_value is None:
                z_value = extent[5]
            zs = np.full_like(xs, z_value)
        
        xyz = np.column_stack((xs, ys, zs))
        
        return xyz, labels_valid, data, (left, right, bottom, top)


def test_spatial_correlation_reduction(base_dir, model_extent, topography_dir):
    """
    Test point reduction using spatial correlation (boundary detection).
    """
    enmap_path = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap', 'EPSG3857_EnMap_result_n4_betajump0.1.tif')
    topo_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')

    if not os.path.exists(enmap_path):
        pytest.skip(f"EnMap file not found at {enmap_path}")

    # 1. Standard extraction (for comparison)
    step_standard = 10
    xyz_std, labels_std, _, _ = _extract_points_from_raster(enmap_path, model_extent, step=step_standard, topo_path=topo_path)
    
    # 2. Spatially reduced extraction
    # We want to keep boundaries dense but interior sparse
    xyz_red, labels_red, data, bounds = _extract_points_spatially_reduced(
        enmap_path, model_extent, step_boundary=50, step_inner=50, topo_path=topo_path
    )
    
    print(f"\n📊 Reduction Comparison:")
    print(f"   Standard (step={step_standard}): {len(xyz_std)} points")
    print(f"   Spatially Reduced: {len(xyz_red)} points")
    print(f"   Reduction factor: {len(xyz_std)/len(xyz_red):.2f}x")
    
    # 3. Visualization
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), sharex=True, sharey=True)
    
    # Standard plot
    axes[0].imshow(data, extent=bounds, cmap='tab10', interpolation='nearest', alpha=0.3)
    axes[0].scatter(xyz_std[:, 0], xyz_std[:, 1], c=labels_std, cmap='tab10', s=1, alpha=0.8)
    axes[0].set_title(f'Standard Regular Grid (n={len(xyz_std)})')
    
    # Reduced plot
    axes[1].imshow(data, extent=bounds, cmap='tab10', interpolation='nearest', alpha=0.3)
    axes[1].scatter(xyz_red[:, 0], xyz_red[:, 1], c=labels_red, cmap='tab10', s=1, alpha=0.8)
    axes[1].set_title(f'Spatially Reduced (n={len(xyz_red)})\nDense at boundaries, sparse in interior')
    
    for ax in axes:
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')

    plt.tight_layout()
    plt.show()
    
    # Assertions
    assert len(xyz_red) < len(xyz_std), "Reduced set should have fewer points than standard if steps are chosen correctly"
    assert len(xyz_red) > 0, "Should have extracted some points"
    assert 1 not in labels_red, "Label 1 should have been ignored"
    assert 3 not in labels_red, "Label 3 should have been combined with 0"
    
    # Check that we have a mix of boundary and inner points
    # (This is a bit heuristic but helps ensure the logic is working)
    with rasterio.open(enmap_path) as src:
        window = from_bounds(model_extent[0], model_extent[2], model_extent[1], model_extent[3], src.transform)
        data = src.read(1, window=window)
        data_mapped = data.copy()
        data_mapped[data_mapped == 3] = 0
        from skimage.segmentation import find_boundaries
        boundaries = find_boundaries(data_mapped, mode='thick')
        
        # Transform xyz back to pixel coordinates to check if they are on boundaries
        inv_transform = ~src.window_transform(window)
        cols, rows = inv_transform * (xyz_red[:, 0], xyz_red[:, 1])
        cols = np.round(cols).astype(int)
        rows = np.round(rows).astype(int)
        
        on_boundary = boundaries[rows, cols]
        n_boundary = np.sum(on_boundary)
        n_inner = len(xyz_red) - n_boundary
        
        print(f"   Points on boundaries: {n_boundary}")
        print(f"   Points in interior: {n_inner}")
        
        assert n_boundary > 0, "Should have some points on boundaries"
        assert n_inner > 0, "Should have some points in interior"
    
    # Save reduced set
    np.save(os.path.join(base_dir, 'reduced_xyz.npy'), xyz_red)
    np.save(os.path.join(base_dir, 'reduced_labels.npy'), labels_red)
    
    print(f"✓ Spatially reduced points saved to reduced_*.npy")


def test_central_body_extraction(base_dir, model_extent, topography_dir):
    """
    Test extraction of points from the center of bodies to minimize spatial correlation.
    """
    enmap_path = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap', 'EPSG3857_EnMap_result_n4_betajump0.1.tif')
    topo_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')

    if not os.path.exists(enmap_path):
        pytest.skip(f"EnMap file not found at {enmap_path}")

    # 1. Spatially reduced extraction (previous approach for comparison)
    xyz_red, labels_red, _, _ = _extract_points_spatially_reduced(
        enmap_path, model_extent, step_boundary=50, step_inner=50, topo_path=topo_path
    )
    
    # 2. Central extraction (new approach)
    # min_distance ensures points are not correlated
    xyz_central, labels_central, data, bounds = _extract_points_central_reduced(
        enmap_path, model_extent, min_distance=25, topo_path=topo_path
    )
    
    print(f"\n📊 Strategy Comparison:")
    print(f"   Boundary-focused reduction: {len(xyz_red)} points")
    print(f"   Central-focused reduction: {len(xyz_central)} points")
    
    # 3. Visualization
    fig, axes = plt.subplots(1, 2, figsize=(18, 8), sharex=True, sharey=True)
    
    # Boundary-focused plot
    axes[0].imshow(data, extent=bounds, cmap='tab10', interpolation='nearest', alpha=0.3)
    axes[0].scatter(xyz_red[:, 0], xyz_red[:, 1], c=labels_red, cmap='tab10', s=10, edgecolors='black', linewidth=0.5)
    axes[0].set_title(f'Boundary-Focused (n={len(xyz_red)})')
    
    # Central plot
    axes[1].imshow(data, extent=bounds, cmap='tab10', interpolation='nearest', alpha=0.3)
    axes[1].scatter(xyz_central[:, 0], xyz_central[:, 1], c=labels_central, cmap='tab10', s=10, edgecolors='black', linewidth=0.5)
    axes[1].set_title(f'Central-Focused (n={len(xyz_central)})\nPoints at local maxima of distance transform')
    
    for ax in axes:
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')

    plt.tight_layout()
    plt.show()
    
    # Assertions
    assert len(xyz_central) > 0, "Should extract at least some points"
    assert not np.any(labels_central == 1), "Label 1 should be excluded"
    assert not np.any(labels_central == 3), "Label 3 should be mapped to 0"
    
    # Final output
    np.save(os.path.join(base_dir, 'central_xyz.npy'), xyz_central)
    np.save(os.path.join(base_dir, 'central_labels.npy'), labels_central)
    print(f"\n✅ Central points saved to 'central_xyz.npy' and 'central_labels.npy'")


def test_extract_reference_points(base_dir, model_extent, simple_geo_model, topography_dir):
    enmap_path = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap', 'EPSG3857_EnMap_result_n4_betajump0.1.tif')
    topo_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')

    if not os.path.exists(enmap_path):
        pytest.skip(f"EnMap file not found at {enmap_path}")

    # Extract points from main model extent
    xyz, labels, data, bounds = _extract_points_from_raster(enmap_path, model_extent, topo_path=topo_path)
    
    print(f"   Extracted {len(xyz)} points from main extent")
    
    # 3. Plotting for verification
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(data, extent=bounds, cmap='tab10', interpolation='nearest')
    plt.colorbar(im, label='Class ID')
    ax.scatter(xyz[:, 0], xyz[:, 1], c='white', s=2, alpha=0.5, label='Extracted Points')
    ax.set_title('Cropped EnMap with Extracted Points (Label 1 ignored, 3&0 combined)')
    ax.set_xlabel('X (m)')
    ax.set_ylabel('Y (m)')
    plt.legend()
    plt.show()

    # Assertions
    assert len(xyz) > 0, "Should have extracted some points"
    assert 1 not in labels, "Label 1 should have been ignored"
    assert 3 not in labels, "Label 3 should have been combined with 0"
    
    # Save results for optimization
    np.save(os.path.join(base_dir, 'extracted_xyz.npy'), xyz)
    np.save(os.path.join(base_dir, 'extracted_labels.npy'), labels)
    
    print(f"✓ Extracted {len(xyz)} points saved to .npy files")

    # Extra extent demonstration
    extra_extent = [
        model_extent[0], 
        model_extent[0] + (model_extent[1] - model_extent[0]) / 2,
        model_extent[2],
        model_extent[2] + (model_extent[3] - model_extent[2]) / 2,
        model_extent[4],
        model_extent[5]
    ]
    
    xyz_extra, labels_extra, _, _ = _extract_points_from_raster(enmap_path, extra_extent, step=5, topo_path=topo_path)
    print(f"✓ Extracted {len(xyz_extra)} points from extra extent (quarter of model)")
