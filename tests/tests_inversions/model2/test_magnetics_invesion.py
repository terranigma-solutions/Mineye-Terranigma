"""Test magnetic inversion using the Soricom fault structural model.

Uses the Soricom fault model from ``examples/01_basic_examples/04_soricom_fault_model.py``
and magnetic point observations extracted from the B1B2 merged raster
(see ``test_extract_magnetic_points.py`` for extraction).
"""

import os

import arviz as az
import numpy as np
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist
from torch.xpu import device

import gempy as gp
import gempy_probability as gpp
from gempy_engine.core.data.geophysics_input import MagneticsInput
from gempy_engine.modules.geophysics.magnetic_gradient import calculate_magnetic_gradient_tensor
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.probabilistic_model import normalize, _modify_orientations
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import generate_multimagnetic_likelihood_hierarchical_per_station
from mineye.GeoModel.model_one.visualization import gempy_viz, plot_many_observed_vs_forward
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig

nc_path = "../../../examples/02_probabilistic_modeling/arviz_data_magnetic_soricom.nc"


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
    prior_key_dips = r'dips'
    prior_key_susceptibility = r'susceptibility'

    def test_run_predictive(self, geophysical_dir, n_samples=50):
        import time
        
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir=geophysical_dir,
        )

        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)

        compute_prior_predictive = True
        if compute_prior_predictive:
            start_time = time.time()
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=magnetic_observations_tensor,
                # n_samples=100,
                n_samples=40,
                plot_trace=False
            )
            elapsed_time = time.time() - start_time
            print(f"Prior predictive sampling completed in {elapsed_time:.2f} seconds ({elapsed_time / 60:.2f} minutes)")

    def test_magnetic_inversion(
            self, geophysical_dir, n_samples=50,
            arviz_data_filename="arviz_data_magnetic_soricom.nc",
    ):
        """Test magnetic inversion using the Soricom fault model and NUTS sampler."""
        print("Soricom magnetic inversion with NUTS...")
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir,
        )

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
        # import pykeops
        # pykeops.config.verbose = True
        # Compute baseline forward magnetics
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
        n_orientations = geo_model.orientations_copy.xyz.shape[0]
        model_priors = {
                self.prior_key_dips          : dist.Normal(
                    loc=torch.full((n_orientations,), 10.0, dtype=float_, device=device),
                    scale=torch.tensor(10.0, dtype=float_, device=device),
                    validate_args=True,
                ),
                self.prior_key_susceptibility: dist.Normal(
                    loc=torch.tensor([0.0, 0.05, 0.001, 0.001], dtype=float_, device=device),
                    scale=torch.tensor(0.03, dtype=float_, device=device),
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

        # 5) Likelihood
        likelihood_fn = generate_multimagnetic_likelihood_hierarchical_per_station(
            norm_params=norm_params,
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
                susceptibilities=np.array([0.0, 0.05, 0.001, 0.001]),
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
        interpolation_input = _modify_orientations(
            samples=samples,
            geo_model=geo_model,
            key=r'dips',
        )

        if TestMagneticInversion.prior_key_susceptibility in samples:
            susceptibilities = samples[TestMagneticInversion.prior_key_susceptibility]
            if geo_model.geophysics_input and geo_model.geophysics_input.magnetics_input:
                geo_model.geophysics_input.magnetics_input.susceptibilities = susceptibilities

        return interpolation_input

    # --- Analysis tests ---

    def test_run_predictive_analysis(self, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), nc_path))
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
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), nc_path))
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        gempy_viz(geo_model, data, n_samples=100)

    def test_run_outlier_detection(self, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), nc_path))

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
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), nc_path))
        data.posterior

        az.plot_posterior(data, var_names=["dips", "susceptibility"])
        plt.show()

        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        plt.rcParams['figure.dpi'] = 72

        axes = az.plot_density(
            data=[data, data.prior],
            var_names=["dips", "susceptibility"],
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

        if hasattr(data, 'prior') and magnetic_response_key in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[magnetic_response_key].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel,
            )

        if hasattr(data, 'posterior') and magnetic_response_key in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.posterior_predictive[magnetic_response_key].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel,
            )

        # gempy_viz(geo_model, data)
