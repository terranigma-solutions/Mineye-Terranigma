import os
import time

import gempy as gp
import gempy_viewer as gpv


def test_simple_model(simple_geo_model):
    """Test reading and computing a geological model."""
    start_time = time.time()
    gp.compute_model(simple_geo_model)
    elapsed_time = time.time() - start_time

    print(f"\n⏱️  Model computation time: {elapsed_time:.2f} seconds")

    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None

    gpv.plot_3d(simple_geo_model, ve=5, image=True)


def test_simple_model_with_topography(simple_geo_model, topography_dir):
    """Test reading and computing a geological model with topography."""
    topography_path = os.path.join(topography_dir,  'topo_reduced_sf0.1.tif')
    gp.set_topography_from_file(
        grid=simple_geo_model.grid,
        filepath=topography_path,
        crop_to_extent=[-695558, simple_geo_model.grid.extent[2],
                        simple_geo_model.grid.extent[1], simple_geo_model.grid.extent[3]]
    )

    start_time = time.time()
    gp.compute_model(simple_geo_model)
    elapsed_time = time.time() - start_time
    print(f"\n⏱️  Model computation time: {elapsed_time:.2f} seconds")

    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None

    gpv.plot_3d(simple_geo_model, ve=5, image=True)
