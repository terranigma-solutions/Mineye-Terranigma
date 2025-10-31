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

    def test_error_propagation(self, simple_geo_model):
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
            samples=(prior_inference_data.prior[r'dips'].values[0, :]),
            update_model_fn=self._update_model_for_plotting,
            gempy_plot=p2d
        )

    def _gravity_precomputations(self, density_plutonites: float, density_sedimentary_host: float,
                                 xy_ravel: np.ndarray, simple_geo_model: gp.data.GeoModel):
        print("Using actual gravity measurement locations...")
        print(f"Using {len(xy_ravel)} actual measurement points")

        print("Computing forward gravity model...")

        # Step 1: Set centered grid
        print("Setting up centered grid...")
        gp.set_centered_grid(
            grid=simple_geo_model.grid,
            centers=xy_ravel,
            resolution=np.array([10, 10, 15]),
            radius=np.array([5000, 5000, 5000])
        )

        # Step 2: Calculate gravity gradient (tz component)
        print("Calculating gravity gradient...")
        gravity_gradient = gp.calculate_gravity_gradient(simple_geo_model.grid.centered_grid)

        # Step 3: Configure geophysics input
        print("Configuring geophysics input...")
        simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
            tz=gravity_gradient,
            densities=np.array([density_sedimentary_host, density_plutonites])  # kg/m³ for different formations,
        )

    def test_gravity_error_propagation(self, simple_geo_model, geophysical_dir):
        import geopandas as gpd
        import os

        # Model parameters
        # Use actual gravity measurement device locations
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

        geo_model, xy_ravel = self._setup_geomodel(gravity_data, simple_geo_model)

        # ============ COMPUTE BASELINE FORWARD MODEL ============
        # CRITICAL: Compute forward model with INITIAL/MEAN parameters before inference
        # This provides the baseline statistics needed to preserve prior variability
        baseline_fw_gravity_np = self._baseline(geo_model)

        norm_params, normalization_method, observed_gravity_ugal = self._normalize(baseline_fw_gravity_np, observed_gravity)

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

        from mineye.GeoModel.geophysics import apply_normalization_torch
        # Post-forward deterministics (extract model outputs)
        # Normalization is applied HERE during inference using the pre-computed parameters
        post_forward_dets = {
                "gravity_response_raw": lambda samples, gm, sol: sol.gravity,  # Store raw gravity
                "gravity_response"    : lambda samples, gm, sol: apply_normalization_torch(-sol.gravity, norm_params),  # Normalized!
                "mean_gravity"        : lambda samples, gm, sol: torch.mean(apply_normalization_torch(-sol.gravity, norm_params)),
                "max_gravity"         : lambda samples, gm, sol: torch.max(apply_normalization_torch(-sol.gravity, norm_params), 0),
        }

        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=self._modify_dips_for_orientations,
            likelihood_fn=None,
            pre_forward_deterministics=pre_forward_dets,
            post_forward_deterministics=post_forward_dets,
            obs_name=None
        )

        n_samples = 10
        prior_inference_data: az.InferenceData = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=[],
            n_samples=n_samples,
            plot_trace=True
        )
        
        # endregion

        # region EXTRACT NORMALIZED SAMPLES ============
        # Extract gravity samples: shape (n_chains, n_samples, n_devices)
        # We use chain 0 for visualization
        # IMPORTANT: These are ALREADY NORMALIZED (done during inference)
        gravity_samples_norm = prior_inference_data.prior[r'gravity_response'].values[0, :]  # (n_samples, n_devices)

        gravity_samples_norm, unit_label = self._plot(gravity_samples_norm, observed_gravity_ugal, xy_ravel)

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
        
        #endreion

    def _plot(self, gravity_samples_norm, observed_gravity_ugal, xy_ravel) -> tuple[str, Any]:
        # Import visualization functions
        from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_map_interpolated
        from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_uncertainty_profiles
        from mineye.GeoModel.plotting.probabilistic_analysis import plot_gravity_with_uncertainty

        # Convert PyTorch tensor to numpy if needed
        if hasattr(gravity_samples_norm, 'numpy'):
            gravity_samples_norm = gravity_samples_norm.numpy()
            
        unit_label = 'Aligned Gravity (μGal)'

        print(f"\n{'=' * 60}")
        print(f"EXTRACTED NORMALIZED SAMPLES")
        print(f"{'=' * 60}")
        print(f"Number of samples: {gravity_samples_norm.shape[0]}")
        print(f"Number of devices: {gravity_samples_norm.shape[1]}")
        print(f"Normalization was applied DURING inference (not post-processing)")
        print(f"All samples use consistent normalization parameters from observed data")
        print(f"{'=' * 60}\n")

        # 1. Comprehensive uncertainty visualization (with normalized data)
        plot_gravity_with_uncertainty(
            gravity_samples=gravity_samples_norm,
            xy_coords=xy_ravel,
            observed_data=observed_gravity_ugal,
            confidence_level=0.95,
            title="Gravity Uncertainty Propagation from Dip Uncertainty"
        )

        # 2. Profile plots with confidence bands (with normalized data)
        plot_gravity_uncertainty_profiles(
            gravity_samples=gravity_samples_norm,
            xy_coords=xy_ravel,
            observed_data=(observed_gravity_ugal),
            n_profiles=4,
            confidence_level=0.95
        )

        # 3. Interpolated uncertainty maps (smoother visualization with normalized data)
        plot_gravity_uncertainty_map_interpolated(
            gravity_samples=gravity_samples_norm,
            xy_coords=xy_ravel,
            observed_data=(observed_gravity_ugal),
            grid_resolution=100
        )
        return gravity_samples_norm, unit_label

    def _baseline(self, geo_model) -> Any:
        print("\n" + "=" * 60)
        print("COMPUTING BASELINE FORWARD MODEL")
        print("=" * 60)
        print("Computing gravity with mean/initial prior parameters...")

        baseline_fw_gravity = geo_model.solutions.gravity

        # Convert to numpy for normalization parameter computation
        if hasattr(baseline_fw_gravity, 'numpy'):
            baseline_fw_gravity_np = baseline_fw_gravity.numpy()
        else:
            baseline_fw_gravity_np = np.array(baseline_fw_gravity)

        print(f"Baseline forward model statistics:")
        print(f"  Mean: {np.mean(baseline_fw_gravity_np):.2f} μGal")
        print(f"  Std:  {np.std(baseline_fw_gravity_np):.2f} μGal")
        print(f"  Min:  {np.min(baseline_fw_gravity_np):.2f} μGal")
        print(f"  Max:  {np.max(baseline_fw_gravity_np):.2f} μGal")
        print("=" * 60 + "\n")
        return -baseline_fw_gravity_np

    def _normalize(self, baseline_fw_gravity_np , observed_gravity):
        from mineye.GeoModel.geophysics import compute_normalization_params

        # Convert observed gravity from mGal to μGal for comparison
        observed_gravity_ugal = observed_gravity * 1000

        # Compute normalization parameters once from observed data AND baseline forward model
        # CRITICAL: Pass baseline_forward_model to preserve prior variability
        normalization_method = 'align_to_reference'
        norm_params = compute_normalization_params(
            reference_data=observed_gravity_ugal,
            baseline_forward_model=baseline_fw_gravity_np,  # CRITICAL: Baseline model
            method=normalization_method,
            verbose=True
        )


        print(f"\n{'=' * 60}")
        print(f"NORMALIZATION PARAMETERS (FIXED for all samples)")
        print(f"{'=' * 60}")
        print(f"Method: {norm_params['method']}")
        print(f"Reference (observed) mean: {norm_params['reference_mean']:.2f} μGal")
        print(f"Reference (observed) std:  {norm_params['reference_std']:.2f} μGal")
        print(f"Baseline forward mean:     {norm_params['baseline_forward_mean']:.2f} μGal")
        print(f"Baseline forward std:      {norm_params['baseline_forward_std']:.2f} μGal")
        print(f"\n✓ Baseline statistics are FIXED (from initial model)")
        print(f"✓ All samples normalized using these fixed parameters")
        print(f"✓ This PRESERVES variability from prior uncertainty")
        print(f"{'=' * 60}\n")
        return norm_params, normalization_method, observed_gravity_ugal

    def _setup_geomodel(self, gravity_data, simple_geo_model: gp.data.GeoModel): 
        geo_model: gp.data.GeoModel = simple_geo_model
        BackendTensor.change_backend_gempy(engine_backend=gp.data.AvailableBackends.PYTORCH)

        xy_ravel = np.column_stack([
                np.array(gravity_data.geometry.x.values),
                np.array(gravity_data.geometry.y.values),
                np.full(len(gravity_data), 0)  # Set Z to surface elevation
        ])

        self._gravity_precomputations(
            density_plutonites=2.9,  # kg/m³
            density_sedimentary_host=2.3,  # kg/m³
            xy_ravel=xy_ravel,
            simple_geo_model=geo_model
        )

        geo_model.interpolation_options.mesh_extraction = False

        gp.set_active_grid(
            grid=geo_model.grid,
            grid_type=[geo_model.grid.GridTypes.CENTERED],
            reset=True
        )
        gp.compute_model(geo_model)
        return geo_model, xy_ravel
