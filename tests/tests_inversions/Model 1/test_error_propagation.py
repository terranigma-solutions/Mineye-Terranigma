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
from gempy_probability.modules.plot.plot_gempy import plot_gempy

seed = 4003
pyro.set_rng_seed(seed)
torch.manual_seed(seed)


def test_error_propagation(simple_geo_model):
    """Test reading and computing a geological model with topography."""

    def _modify_z_for_surface_point1(
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

    def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: float, sample_idx: int):
        xyz = np.zeros((1, 3))
        xyz[0, 2] = sample_value
        world_coord = geo_model.input_transform.apply_inverse(xyz)

        # Modify the surface point
        gp.modify_surface_points(
            geo_model=geo_model,
            slice=0,
            Z=world_coord[0, 2]
        )

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
        set_interp_input_fn=_modify_z_for_surface_point1,
        likelihood_fn=None,
        obs_name=None
    )

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

    plot_gempy(
        geo_model=geo_model,
        n_samples=n_samples,
        samples=(prior_inference_data.prior[r'$\mu_{top}$'].values[0, :]),
        update_model_fn=-update_model_for_plotting,
        gempy_plot=p2d
    )
