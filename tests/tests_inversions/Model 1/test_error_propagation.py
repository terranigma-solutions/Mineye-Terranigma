import gempy_viewer as gpv

import gempy as gp
import numpy as np

from gempy_engine.core.backend_tensor import BackendTensor

import pyro.distributions as dist
import torch

from pyro.infer import Predictive
import arviz as az
import matplotlib.pyplot as plt

import pyro
from pyro.distributions import Distribution
import gempy_probability as gpp
from gempy_engine.core.data.interpolation_input import InterpolationInput


seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)

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

    model_priors = {
            r'$\mu_{top}$': dist.Normal(
                loc=geo_model.surface_points_copy_transformed.xyz[0, 2],
                scale=torch.tensor(0.001, dtype=torch.float64)
            )
    }

    prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model(
        priors=model_priors,
        set_interp_input_fn=modify_z_for_surface_point1,
        likelihood_fn=None,
        obs_name=None
    )

    # TODO: Code to define the method and all of that

    n_samples = 50
    prior_inference_data: az.InferenceData = gpp.run_predictive(
        prob_model=prob_model,
        geo_model=geo_model,
        y_obs_list=[],
        n_samples=n_samples,
        plot_trace=True
    )

    p2d = gpv.plot_2d(
        model=geo_model,
        show_topography=False,
        legend=False,
        show_lith=False,
        show_data=False,
        show=False
    )

    posterior_top_mean_z: np.ndarray = (prior_inference_data.prior[r'$\mu_{top}$'].values[0, :])
    xyz = np.zeros((posterior_top_mean_z.shape[0], 3))
    xyz[:, 2] = posterior_top_mean_z
    world_coord = geo_model.input_transform.apply_inverse(xyz)
    
    for i in np.linspace(0, n_samples-1, n_samples,).astype(int):
        gp.modify_surface_points(
            geo_model=geo_model,
            slice=0,
            Z=world_coord[i, 2]
        )
        gp.compute_model(gempy_model=geo_model)

        from gempy_viewer.API._plot_2d_sections_api import plot_sections
        from gempy_viewer.core.data_to_show import DataToShow
        plot_sections(
            gempy_model=geo_model,
            sections_data=p2d.section_data_list,
            data_to_show=DataToShow(
                n_axis=1,
                show_data=True,
                show_surfaces=True,
                show_lith=False
            ),
            kwargs_boundaries={
                    "linewidth": 0.5,
                    "alpha"    : 0.1,
            },
            kwargs_surface_points={
                    'alpha': 0.1  # semi-transparent
            },
        )

    p2d.fig.show()


    
    
