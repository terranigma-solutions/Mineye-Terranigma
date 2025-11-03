from typing import Any

import arviz as az
import numpy as np
import torch
from pyro import distributions as dist
from pyro.distributions import Distribution

import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv
from gempy_engine.core.backend_tensor import BackendTensor
from gempy_engine.core.data.interpolation_input import InterpolationInput
from gempy_probability.modules.plot.plot_gempy import plot_gempy
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_1_aux import setup_geomodel, baseline, normalize, plot

# noinspection PyUnusedImports
from tests import conftest


class TestErrorPropagationDips:
    prior_key = r'dips'

    @staticmethod
    def _modify_dips_for_orientations(
            samples: dict[str, Distribution],
            geo_model: gp.data.GeoModel,
    ) -> InterpolationInput:
        # # TODO: We can make a factory for this type of functions

        from gempy.modules.data_manipulation import interpolation_input_from_structural_frame
        from gempy.modules.data_manipulation.manipulate_points import compute_adp_from_gradients, convert_orientation_to_pole_vector

        interp_input: InterpolationInput = interpolation_input_from_structural_frame(geo_model)
        smaples = samples[TestErrorPropagationDips.prior_key]

        og_gradients = interp_input.orientations.dip_gradients
        azimuth, dip, polarity = compute_adp_from_gradients(
            G_x=og_gradients[:, 0],
            G_y=og_gradients[:, 1],
            G_z=og_gradients[:, 2],
        )

        gradients = convert_orientation_to_pole_vector(
            azimuth=azimuth,
            dip=smaples,
            polarity=polarity
        )

        interp_input.orientations.dip_gradients = gradients
        return interp_input

    @staticmethod
    def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: float, sample_idx: int):
        # # Modify the surface point
        gp.modify_orientations(
            geo_model=geo_model,
            dip=sample_value,
        )

    def test_error_propagation(self, simple_geo_model, n_samples=50):
        """Test reading and computing a geological model with topography."""

        geo_model = simple_geo_model
        BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

        mean_orientations = torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10
        model_priors = {
                r'dips': dist.Normal(
                    loc=mean_orientations,
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                )
        }

        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=self._modify_dips_for_orientations,
            likelihood_fn=None,
            obs_name=None
        )

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
            n_samples=20,
            samples=(prior_inference_data.prior[r'dips'].values[0, :]),
            update_model_fn=self._update_model_for_plotting,
            gempy_plot=p2d
        )

    def test_gravity_error_propagation(self, simple_geo_model, geophysical_dir, n_samples=50):
        import geopandas as gpd
        import os

        # Model parameters
        # Use actual gravity measurement device locations
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model, )

        # ============ COMPUTE BASELINE FORWARD MODEL ============
        # CRITICAL: Compute forward model with INITIAL/MEAN parameters before inference
        # This provides the baseline statistics needed to preserve prior variability
        baseline_fw_gravity_np = baseline(geo_model)

        observed_gravity_ugal = observed_gravity * 1000
        norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)

        # region PYRO MODEL SETUP ============
        mean_orientations = torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10
        model_priors = {
                r'dips': dist.Normal(
                    loc=mean_orientations,
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                )
        }

        # Pre-forward deterministics (parameter transformations)
        pre_forward_dets = {
                "dips_degrees": lambda samples, gm: samples["dips"],  # Just pass through
        }

        # Post-forward deterministics (extract model outputs)
        # Normalization is applied HERE during inference using the pre-computed parameters

        post_forward_dets = {
                "gravity_response_raw": lambda samples, gm, sol: sol.gravity,  # Store raw gravity
                "gravity_response"    : lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),  # Normalized!
                "mean_gravity"        : lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
                "max_gravity"         : lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
        }

        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=self._modify_dips_for_orientations,
            likelihood_fn=None,
            pre_forward_deterministics=pre_forward_dets,
            post_forward_deterministics=post_forward_dets,
            obs_name=None
        )

        prior_inference_data: az.InferenceData = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=[],
            n_samples=n_samples,
            plot_trace=True
        )

        # endregion

        # region EXTRACT NORMALIZED SAMPLES ============
        gravity_samples_norm = prior_inference_data.prior[r'gravity_response'].values[0, :]  # (n_samples, n_devices)
        gravity_samples_norm, unit_label = plot(
            gravity_samples_norm=gravity_samples_norm,
            observed_gravity_ugal=observed_gravity_ugal,
            xy_ravel=xy_ravel
        )

        # Print final summary
        print(f"\n{'=' * 60}")
        print(f"INFERENCE COMPLETE - NORMALIZED RESULTS")
        print(f"{'=' * 60}")
        print(f"Method: {norm_params['method']}")
        print(f"Parameters: {norm_params}")
        print(f"Unit label: {unit_label}")
        print(f"Samples shape: {gravity_samples_norm.shape}")
        print(f"✓ Normalization applied DURING inference (in post_forward_deterministics)")
        print(f"✓ All {gravity_samples_norm.shape[0]} samples use consistent parameters")
        print(f"✓ Parameters computed from observed data before inference")
        print(f"{'=' * 60}\n")

        # endregion

        # region gempy representation
        gp.set_active_grid(
            grid=geo_model.grid,
            grid_type=[geo_model.grid.GridTypes.OCTREE],
            reset=True
        )

        geo_model.geophysics_input = None

        gp.compute_model(gempy_model=geo_model)

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
            n_samples=20,
            samples=(prior_inference_data.prior[r'dips'].values[0, :]),
            update_model_fn=self._update_model_for_plotting,
            gempy_plot=p2d
        )

        # endregion
   