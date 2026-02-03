import os
import time

import numpy as np
import pytest
import rasterio
from matplotlib import pyplot as plt

import gempy as gp
import gempy_viewer as gpv

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


def test_extract_reference_points(base_dir, model_extent, simple_geo_model):
    enmap_path = os.path.join(base_dir, 'examples', 'Data', 'Segmentation_Input_Data', 'Enmap', 'EPSG3857_EnMap_result_n4_betajump0.1.tif')

    if not os.path.exists(enmap_path):
        pytest.skip(f"EnMap file not found at {enmap_path}")

    with rasterio.open(enmap_path) as src:
        # 1. Define the cropping window based on model_extent
        # model_extent: [min_x, max_x, min_y, max_y, min_z, max_z]
        left, right, bottom, top = model_extent[0], model_extent[1], model_extent[2], model_extent[3]
        
        from rasterio.windows import from_bounds
        window = from_bounds(left, bottom, right, top, src.transform)
        
        # Read the data within the window
        data = src.read(1, window=window)
        transform = src.window_transform(window)
        
        print(f"   Cropped Shape: {data.shape}")
        print(f"   Cropped Bounds: {left}, {right}, {bottom}, {top}")

        # 2. Extract well distributed amount of points
        # We can use a regular grid or random sampling. 
        # Let's use a step-based approach to get a well-distributed sample.
        step = 10  # Adjust step for density
        rows, cols = data.shape
        
        row_indices = np.arange(0, rows, step)
        col_indices = np.arange(0, cols, step)
        
        ii, jj = np.meshgrid(row_indices, col_indices, indexing='ij')
        
        # Flatten indices
        ii_flat = ii.flatten()
        jj_flat = jj.flatten()
        
        # Get labels
        labels = data[ii_flat, jj_flat]
        
        # Filter out NaNs if any
        valid_mask = ~np.isnan(labels)
        ii_valid = ii_flat[valid_mask]
        jj_valid = jj_flat[valid_mask]
        labels_valid = labels[valid_mask]
        
        # Get xy coordinates
        # rasterio.transform.xy(transform, rows, cols, offset='center')
        xs, ys = rasterio.transform.xy(transform, ii_valid.tolist(), jj_valid.tolist())
        xs = np.array(xs)
        ys = np.array(ys)
        
        # For Z, we might want to use topography or a fixed depth. 
        # Given it's a surface map, Z should probably come from topography or be at the top.
        # If we have the model, we can try to get topography at these points.
        # For now, let's assume z=0 or top of model if not specified.
        zs = np.full_like(xs, model_extent[5]) # Using max_z
        
        xyz = np.column_stack((xs, ys, zs))
        
        print(f"   Extracted {len(xyz)} points")
        
        # 3. Plotting for verification
        fig, ax = plt.subplots(figsize=(10, 8))
        im = ax.imshow(data, extent=(left, right, bottom, top), cmap='tab10', interpolation='nearest')
        plt.colorbar(im, label='Class ID')
        ax.scatter(xyz[:, 0], xyz[:, 1], c='white', s=2, alpha=0.5, label='Extracted Points')
        ax.set_title('Cropped EnMap with Extracted Points')
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        plt.legend()
        plt.show()

        # Assertions
        assert len(xyz) > 0, "Should have extracted some points"
        assert xyz.shape[1] == 3, "xyz should have 3 columns"
        assert len(labels_valid) == len(xyz), "Labels and xyz should have same length"
        
        # Save results for optimization (optional, but requested format was array or xyz and labels)
        np.save(os.path.join(base_dir, 'extracted_xyz.npy'), xyz)
        np.save(os.path.join(base_dir, 'extracted_labels.npy'), labels_valid)
        
        print(f"✓ Extracted {len(xyz)} points saved to .npy files")
