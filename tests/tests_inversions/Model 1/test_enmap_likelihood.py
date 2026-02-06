import os
import time

import numpy as np
import pytest
import rasterio
from matplotlib import pyplot as plt

import gempy as gp
import gempy_viewer as gpv


def test_simple_model_with_topography(simple_geo_model, topography_dir):
    """Test reading and computing a geological model with topography."""
    topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')
    gp.set_topography_from_file(
        grid=simple_geo_model.grid,
        filepath=topography_path,
        crop_to_extent=[
                simple_geo_model.grid.extent[0],
                simple_geo_model.grid.extent[2],
                simple_geo_model.grid.extent[1],
                simple_geo_model.grid.extent[3]
        ]
    )

    start_time = time.time()

    simple_geo_model.interpolation_options.evaluation_options.number_octree_levels_surface = 3
    gp.compute_model(simple_geo_model)
    elapsed_time = time.time() - start_time
    print(f"\n⏱️  Model computation time: {elapsed_time:.2f} seconds")

    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None

    gpv.plot_3d(simple_geo_model, ve=5, image=True,
                kwargs_pyvista_bounds={
                        'show_xlabels': False,
                        'show_ylabels': False,
                        'show_zlabels': False
                }
                )
    gpv.plot_2d(model=simple_geo_model, section_names=['topography'], show_topography=True)


def test_read_EnMap(base_dir, model_extent, simple_geo_model):
    """Test reading EnMap files and plotting them.

    This test reads the EnMap segmentation results (which are at the same
    resolution as the GemPy model) and creates visualizations.
    """
    if rasterio is None:
        pytest.skip("rasterio is required for reading EnMap files")

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
