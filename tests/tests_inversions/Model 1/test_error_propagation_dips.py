import arviz as az
import torch
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
import gempy_viewer as gpv
from gempy_engine.core.backend_tensor import BackendTensor
from gempy_probability.modules.plot.plot_gempy import plot_gempy
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.model_setup import baseline, setup_geomodel
from mineye.GeoModel.model_one.probabilistic_model import normalize, create_orientation_modifier
from mineye.GeoModel.model_one.visualization import plot, gempy_viz
# noinspection PyUnusedImports
from tests import conftest


class TestErrorPropagationDips:
    prior_key = r'dips'


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
            set_interp_input_fn=create_orientation_modifier(key=TestErrorPropagationDips.prior_key),
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
        geo_model.interpolation_options.sigmoid_slope = 100

        observed_gravity_ugal = observed_gravity * 1000
        norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)

        # region PYRO MODEL SETUP ============
        mean_orientations = torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10
        model_priors = {
                r'dips': dist.Normal(
                    loc=mean_orientations,
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                ).to_event(1)
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
            set_interp_input_fn=create_orientation_modifier(key=TestErrorPropagationDips.prior_key),
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

        gempy_viz(geo_model, prior_inference_data)

  