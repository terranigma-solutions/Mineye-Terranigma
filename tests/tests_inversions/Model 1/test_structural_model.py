import gempy as gp
import gempy_viewer as gpv


def test_read_data(simple_geo_model, geomodel_dir):
    """Test reading and computing a geological model."""
    # Note: topography_path needs to be defined - uncomment and fix path if topography file exists
    # topography_path = os.path.join(geomodel_dir, 'path', 'to', 'topography.csv')
    # gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topography_path)

    # gp.set_topography_from_file(grid=simple_geo_model.grid, filepath=topography_path)
    gp.compute_model(simple_geo_model)

    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None
     
    gpv.plot_3d(simple_geo_model, ve=5, image=False)
    
    