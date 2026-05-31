"""Test magnetic inversion using the Soricom fault structural model.

Uses the Soricom fault model from ``examples/01_basic_examples/04_soricom_fault_model.py``
and magnetic point observations extracted from the B1B2 merged raster
(see ``test_extract_magnetic_points.py`` for extraction).
"""

import os

import arviz as az
import numpy as np
import pytest
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
from gempy_engine.core.data.geophysics_input import MagneticsInput
from gempy_engine.modules.geophysics.magnetic_gradient import calculate_magnetic_gradient_tensor
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multimagnetic_likelihood_fixed_std
from mineye.GeoModel.model_one.visualization import plot_many_observed_vs_forward, probability_fields_for
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig

nc_path ="/home/leguark/Projects/gempy_project/Mineye-Terranigma/tests/tests_inversions/model2/arviz_data_magnetic_soricom_z_1779986757.nc" 

RUN_SIMULATION = False

_original_z_coords = None


def _update_model_for_plotting(geo_model: gp.data.GeoModel, sample_value: np.ndarray, sample_idx: int):
    global _original_z_coords
    if _original_z_coords is None:
        _original_z_coords = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)

    scale_z = geo_model.input_transform.scale[2]
    shifts_m = sample_value / scale_z  # Convert relative shifts in transformed coordinates to meters

    new_z = _original_z_coords.copy()
    # Point 0: Main_Fault (no change)
    # Points 1-12: host_rock -> shifted by shifts_m[0] (layer wide)
    new_z[1:13] = _original_z_coords[1:13] + shifts_m[0]
    # Points 13-21: chromite lense -> shifted by shifts_m[1:] (independent)
    new_z[13:22] = _original_z_coords[13:22] + shifts_m[1:]

    gp.modify_surface_points(
        geo_model=geo_model,
        Z=new_z
    )


def gempy_viz(geo_model: gp.data.GeoModel, prior_inference_data: az.InferenceData,
              n_samples=20, ve=3):
    from gempy_probability.modules.plot.plot_gempy import plot_gempy
    import gempy_viewer as gpv

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
        show=False,
        ve=ve
    )

    # Cache original Z *before* the prior loop mutates the model in-place
    original_z = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)
    global _original_z_coords
    _original_z_coords = None

    plot_gempy(
        geo_model=geo_model,
        n_samples=n_samples,
        samples=(prior_inference_data.prior['surface_points_z'].values[0, :]),
        update_model_fn=_update_model_for_plotting,
        gempy_plot=p2d,
    )

    if hasattr(prior_inference_data, 'posterior'):
        # Restore model to original Z so the posterior overlay is correctly anchored
        gp.modify_surface_points(geo_model=geo_model, Z=original_z)
        _original_z_coords = None
        n_surfaces = len(geo_model.structural_frame.elements_colors_contacts)
        plot_gempy(
            geo_model=geo_model,
            n_samples=n_samples,
            samples=(prior_inference_data.posterior['surface_points_z'].values[0, :]),
            update_model_fn=_update_model_for_plotting,
            gempy_plot=p2d,
            contour_colors=[default_red] * n_surfaces,
        )

    return p2d


def _create_soricom_geomodel():
    """Build the Soricom fault geomodel (no topography, as requested)."""
    geo_model = gp.create_geomodel(
        project_name=SoricomModelConfig.PROJECT_NAME,
        extent=SoricomModelConfig.EXTENT,
        refinement=SoricomModelConfig.REFINEMENT,
        importer_helper=gp.data.ImporterHelper(
            path_to_orientations=paths.get_soricom_orientations(),
            path_to_surface_points=paths.get_soricom_formation_points(),
        ),
    )

    geo_model.grid = geo_model.grid.init_octree_grid(
        extent=SoricomModelConfig.EXTENT,
        octree_levels=SoricomModelConfig.REFINEMENT,
    )
    geo_model.interpolation_options.number_octree_levels_surface = 2

    gp.map_stack_to_surfaces(
        gempy_model=geo_model,
        mapping_object=SoricomModelConfig.SURFACE_MAPPING,
    )

    geo_model.structural_frame.structural_groups[
        SoricomModelConfig.FAULT_GROUP_INDEX
    ].structural_relation = gp.data.StackRelationType.FAULT
    geo_model.structural_frame.fault_relations = SoricomModelConfig.FAULT_RELATIONS_MATRIX

    return geo_model


class TestMagneticInversion:
    prior_key_susceptibility = r'susceptibility'
    prior_key_surface_points_z = r'surface_points_z'

    def test_run_predictive(self, geophysical_dir, n_samples=50):
        if not RUN_SIMULATION:
            pytest.skip("Skipping expensive predictive simulation.")
        import time
        import gempy_viewer as gpv

        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir=geophysical_dir,
        )

        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)

        start_time = time.time()
        prior_inference_data: az.InferenceData = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=magnetic_observations_tensor,
            n_samples=100,
            plot_trace=True
        )
        elapsed_time = time.time() - start_time
        print(f"Prior predictive sampling completed in {elapsed_time:.2f} seconds ({elapsed_time / 60:.2f} minutes)")

        if False:
            geo_model.interpolation_options.mesh_extraction = True
            gp.set_active_grid(
                grid=geo_model.grid,
                grid_type=[geo_model.grid.GridTypes.OCTREE],
                reset=True
            )

            geo_model.geophysics_input = None
            draw = 15
            global _original_z_coords
            _original_z_coords = None
            _update_model_for_plotting(
                geo_model=geo_model,
                sample_value=(prior_inference_data.prior[self.prior_key_surface_points_z].values[0, :][draw]),
                sample_idx=draw
            )
            gp.compute_model(gempy_model=geo_model)
            gpv.plot_3d(geo_model)

            draw = 10
            _update_model_for_plotting(
                geo_model=geo_model,
                sample_value=(prior_inference_data.prior[self.prior_key_surface_points_z].values[0, :][draw]),
                sample_idx=draw
            )
            gp.compute_model(gempy_model=geo_model)
            gpv.plot_3d(geo_model)

    def test_magnetic_inversion(
            self, geophysical_dir, n_samples=50,
            arviz_data_filename="arviz_data_magnetic_soricom_z.nc",
    ):
        """Test magnetic inversion using the Soricom fault model and NUTS sampler with Z position priors."""
        if not RUN_SIMULATION:
            pytest.skip("Skipping expensive magnetic inversion simulation.")
        print("Soricom magnetic inversion (Z-priors) with NUTS...")
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir,
        )
        import time
        arviz_data_filename = f"arviz_data_magnetic_soricom_z_{int(time.time())}.nc"
        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)

        # Prior predictive
        if True:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=magnetic_observations_tensor,
                n_samples=100,
                plot_trace=True,
            )

        # NUTS inference
        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=magnetic_observations_tensor,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                num_samples=200,
                warmup_steps=200,
            ),
            plot_trace=True,
            run_posterior_predictive=True,
        )

        data.extend(prior_inference_data)
        data.to_netcdf(os.path.join(os.path.dirname(__file__), arviz_data_filename))

    def _create_probabilistic_model(self, geophysical_dir):
        # 1) Read magnetic point data
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)

        # 2) Build Soricom geomodel + setup magnetic grid
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        geo_model.interpolation_options.sigmoid_slope = 100
        baseline_fw_magnetics_np = self._baseline_magnetics(geo_model)

        norm_params = normalize(
            baseline_fw_gravity_np=baseline_fw_magnetics_np,
            observed_gravity=observed_magnetics_nt,
            method="align_to_reference",
            extrapolation_buffer=0.3,
        )
        device= torch.device("cuda" if torch.cuda.is_available() and True else "cpu")
        float_ = torch.float32
        # 3) Priors
        model_priors = {
            self.prior_key_surface_points_z: dist.Normal(
                loc=torch.tensor([0.0] * 10, dtype=float_, device=device),
                scale=torch.tensor([15.0] * 10, dtype=float_, device=device) * geo_model.input_transform.scale[2],
                validate_args=True,
            ).to_event(1),
            self.prior_key_susceptibility: dist.LogNormal(
                loc=torch.tensor([-9.21, -9.21, -0.69, -9.21], dtype=float_, device=device),
                scale=torch.tensor([0.1, 0.1, 0.2, 0.1], dtype=float_, device=device),
            ).to_event(1),
        }

        # 4) Deterministics
        post_forward_dets = {
                "magnetic_response_raw": lambda samples, gm, sol: sol.magnetics,
                "magnetic_response"    : lambda samples, gm, sol: align_forward_to_observed(
                    sol.magnetics, norm_params,
                ),
                "mean_magnetics"       : lambda samples, gm, sol: torch.mean(
                    align_forward_to_observed(sol.magnetics, norm_params),
                ),
                "max_magnetics"        : lambda samples, gm, sol: torch.max(
                    align_forward_to_observed(sol.magnetics, norm_params), 0,
                ),
        }

        # 5) Likelihood — fixed std to let NUTS explore parameters freely
        likelihood_fn = generate_multimagnetic_likelihood_fixed_std(
            norm_params=norm_params,
            sigma_value=150.0,  # ~2x baseline residual std, reasonable nT noise floor
        )

        # 6) Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=self._set_magnetic_priors,
            likelihood_fn=likelihood_fn,
            obs_name="Magnetic Measurement (Soricom)",
        )

        return geo_model, observed_magnetics_nt, prob_model

    @staticmethod
    def _read_magnetics(geophysical_dir):
        import os as _os
        xyz_path = _os.path.join(_os.path.dirname(__file__), 'soricom_magnetic_xyz_adaptive.npy')
        mag_path = _os.path.join(_os.path.dirname(__file__), 'soricom_magnetic_values_adaptive.npy')
        xyz = np.load(xyz_path)
        mag = np.load(mag_path)
        # Subsample 20 points for inversion (random subset)
        rng = np.random.default_rng(42)
        idx = rng.choice(len(xyz), size=min(20, len(xyz)), replace=False)
        sampled_xyz = xyz[idx]
        # IGRF: forward model outputs TMI anomalies, so subtract IGRF from raw TMI values
        igrf_intensity_nT = 47500.0
        observed_magnetics = mag[idx] - igrf_intensity_nT
        return sampled_xyz, observed_magnetics

    @staticmethod
    def _setup_magnetic_geomodel(magnetic_xyz):
        geo_model = _create_soricom_geomodel()

        xy_ravel = magnetic_xyz  # Already [X, Y, Z] in EPSG:32634

        gp.set_centered_grid(
            grid=geo_model.grid,
            centers=xy_ravel,
            resolution=np.array([10, 10, 15]),
            radius=np.array([2000, 5000, 2000]),
        )

        gradient_tensor_dict = calculate_magnetic_gradient_tensor(
            centered_grid=geo_model.grid.centered_grid,
            igrf_params={
                    "inclination": 57.0,
                    "declination": 4.0,
                    "intensity"  : 47500.0,
            },
            compute_tmi=True,
            units_nT=True,
        )

        geo_model.geophysics_input = gp.data.GeophysicsInput(
            magnetics_input=MagneticsInput(
                mag_kernel=gradient_tensor_dict['tmi_kernel'],
                susceptibilities=np.array([0.0001, 0.0001, 0.5, 0.0001]),
                igrf_params={
                        "inclination": gradient_tensor_dict['inclination'],
                        "declination": gradient_tensor_dict['declination'],
                        "intensity"  : gradient_tensor_dict['intensity'],
                },
            ),
        )

        geo_model.interpolation_options.mesh_extraction = False

        gp.set_active_grid(
            grid=geo_model.grid,
            grid_type=[geo_model.grid.GridTypes.CENTERED],
            reset=True,
        )
        gp.compute_model(geo_model)

        return geo_model, xy_ravel

    @staticmethod
    def _baseline_magnetics(geo_model):
        sol = gp.compute_model(geo_model)
        if sol.magnetics is None:
            raise RuntimeError(
                "Magnetics forward model returned None — check susceptibility "
                "array length matches number of surfaces + basement"
            )
        return sol.magnetics.detach().cpu().numpy() if hasattr(sol.magnetics, 'detach') else sol.magnetics

    @staticmethod
    def _set_magnetic_priors(samples: dict, geo_model: gp.data.GeoModel):
        from gempy.modules.data_manipulation import interpolation_input_from_structural_frame
        interp_input = interpolation_input_from_structural_frame(geo_model)

        if TestMagneticInversion.prior_key_surface_points_z in samples:
            shifts = samples[TestMagneticInversion.prior_key_surface_points_z]
            coords = interp_input.surface_points.sp_coords.clone()

            # Index 1:13 are the 12 host_rock points -> shifted by shifts[0] (layer wide)
            coords[1:13, 2] = coords[1:13, 2] + shifts[0]
            # Index 13:22 are the 9 chromite lense points -> shifted by shifts[1:] (independent)
            coords[13:22, 2] = coords[13:22, 2] + shifts[1:]

            interp_input.surface_points.sp_coords = coords

        if TestMagneticInversion.prior_key_susceptibility in samples:
            susceptibilities = samples[TestMagneticInversion.prior_key_susceptibility]
            if geo_model.geophysics_input and geo_model.geophysics_input.magnetics_input:
                geo_model.geophysics_input.magnetics_input.susceptibilities = susceptibilities

        return interp_input

    # --- Analysis tests ---

    def test_run_predictive_analysis(self, geophysical_dir):
        if not os.path.exists(nc_path):
            pytest.skip(f"Pre-computed inversion NetCDF file not found at: {nc_path}")
        data = az.from_netcdf(nc_path)
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        observed_norm = observed_magnetics_nt
        magnetic_response = r'$\mu_{magnetics}$'
        forward_norm = data.prior[magnetic_response].mean(axis=1)
        many_forward_norm = data.posterior_predictive[magnetic_response].values[
            0, -40:-20
        ]

        plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)

    def test_run_kde_sections(self, geophysical_dir):
        if not os.path.exists(nc_path):
            pytest.skip(f"Pre-computed inversion NetCDF file not found at: {nc_path}")
        data = az.from_netcdf(nc_path)
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)
        geo_model.interpolation_options.number_octree_levels = 5

        gempy_viz(geo_model, data, n_samples=49)

    def test_run_outlier_detection(self, geophysical_dir):
        if not os.path.exists(nc_path):
            pytest.skip(f"Pre-computed inversion NetCDF file not found at: {nc_path}")
        data = az.from_netcdf(nc_path)

        posterior_sigmas = data.posterior_predictive["sigma_stations"].values
        station_noise_mean = posterior_sigmas.mean(axis=(0, 1))
        sigma_global_mean = station_noise_mean.mean()
        problematic = np.where(station_noise_mean > 2 * sigma_global_mean)[0]
        print(f"Check stations: {problematic}")

        axes = az.plot_density(
            data=[data, data.prior],
            var_names=["sigma_stations"],
            filter_vars="like",
            hdi_prob=0.9999,
            shade=.2,
            data_labels=["Posterior", "Prior"],
            colors=[default_red, default_blue],
        )

        plt.show()

    def test_run_analysis(self, geophysical_dir):
        if not os.path.exists(nc_path):
            pytest.skip(f"Pre-computed inversion NetCDF file not found at: {nc_path}")
        data = az.from_netcdf(nc_path)

        az.plot_posterior(data, var_names=[self.prior_key_surface_points_z, "susceptibility"])
        plt.show()

        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        plt.rcParams['figure.dpi'] = 72

        axes = az.plot_density(
            data=[data, data.prior],
            var_names=[self.prior_key_surface_points_z, "susceptibility"],
            filter_vars="like",
            hdi_prob=0.9999,
            shade=.2,
            data_labels=["Posterior", "Prior"],
            colors=[default_red, default_blue],
        )

        plt.show()

        magnetic_response_key = r'$\mu_{magnetics}$'
        plot_geophysics_comparison(
            forward_norm=data.prior[magnetic_response_key].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,
            xy_ravel=xy_ravel,
        )

        plot_geophysics_comparison(
            forward_norm=data.posterior_predictive[magnetic_response_key].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,
            xy_ravel=xy_ravel,
        )

    def test_probability_plots(self, geophysical_dir):
        """Generate prior and posterior probability fields for the Soricom magnetic inversion."""
        if not os.path.exists(nc_path):
            pytest.skip(f"Pre-computed inversion NetCDF file not found at: {nc_path}")
        data = az.from_netcdf(nc_path)

        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)
        geo_model.interpolation_options.number_octree_levels = 5

        topography_path = paths.get_soricom_dem_path()

        print("\nComputing prior probability fields...")
        original_z = geo_model.surface_points_copy.df['Z'].to_numpy(copy=True)
        global _original_z_coords
        _original_z_coords = None
        probability_fields_for(
            geo_model=geo_model,
            inference_data=data.prior,
            topography_path=topography_path,
            var_name=self.prior_key_surface_points_z,
            update_model_fn=_update_model_for_plotting,
            ve=1,
        )

        if hasattr(data, 'posterior'):
            # Restore model Z so posterior OnlineProbability init uses baseline lithology
            gp.modify_surface_points(geo_model=geo_model, Z=original_z)
            _original_z_coords = None
            print("\nComputing posterior probability fields...")
            probability_fields_for(
                geo_model=geo_model,
                inference_data=data.posterior,
                topography_path=topography_path,
                var_name=self.prior_key_surface_points_z,
                update_model_fn=_update_model_for_plotting,
                ve=1,
            )
            
    
    def test_basic_models(self):
        import gempy_viewer as gpv
        geo_model = _create_soricom_geomodel()
        gp.compute_model(geo_model)
        gpv.plot_3d(geo_model)

