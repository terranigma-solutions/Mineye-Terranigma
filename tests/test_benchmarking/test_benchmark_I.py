import gempy as gp
import gempy_viewer as gpv


def test_simple_model(simple_geo_model):
    """Test reading and computing a geological model."""
    gp.compute_model(simple_geo_model)
    # Add assertions here to verify the model is computed correctly
    assert simple_geo_model is not None

    if False:
        gpv.plot_3d(simple_geo_model, ve=5, image=False)
