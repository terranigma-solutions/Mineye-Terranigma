import gempy_viewer as gpv

import gempy as gp
from gempy_engine.core.backend_tensor import BackendTensor

import pyro.distributions as dist
import torch

from pyro.infer import Predictive
import arviz as az
import matplotlib.pyplot as plt

from pyro.distributions import Distribution
import gempy_probability as gpp
from gempy_engine.core.data.interpolation_input import InterpolationInput


def modify_z_for_surface_point1(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
) -> InterpolationInput:
    # TODO: We can make a factory for this type of functions
    prior_key = r'$\mu_{top}$'

    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame
    interp_input = interpolation_input_from_structural_frame(geo_model)
    new_tensor: torch.Tensor = torch.index_put(
        input=interp_input.surface_points.sp_coords,
        indices=(torch.tensor([0]), torch.tensor([2])),  # * This has to be Tensors
        values=(samples[prior_key])
    )
    interp_input.surface_points.sp_coords = new_tensor
    return interp_input


def test_error_propagation(simple_geo_model):
    """Test reading and computing a geological model with topography."""
    # TODO: Code to define the distributions of dips
    geo_model = simple_geo_model
    BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

    normal = dist.Normal(
        loc=(geo_model.surface_points_copy_transformed.xyz[0, 2]),
        scale=torch.tensor(0.1, dtype=torch.float64)
    )

    model_priors = {
            r'$\mu_{top}$': dist.Normal(
                loc=geo_model.surface_points_copy_transformed.xyz[0, 2],
                scale=torch.tensor(0.1, dtype=torch.float64)
            )
    }

    prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model(
        priors=model_priors,
        set_interp_input_fn=modify_z_for_surface_point1,
        likelihood_fn=None,
        obs_name=None
    )

    # TODO: Code to define the method and all of that

    prior_inference_data: az.InferenceData = gpp.run_predictive(
        prob_model=prob_model,
        geo_model=geo_model,
        y_obs_list=[],
        n_samples=50,
        plot_trace=True
    )


    # TODO: Code to visualize the results

    az.plot_trace(prior_inference_data.prior)
    plt.show()
