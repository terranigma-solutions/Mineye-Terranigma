"""Test magnetic inversion using the Soricom fault structural model.

Uses the Soricom fault model from ``examples/01_basic_examples/04_soricom_fault_model.py``
and the cleaned magnetic point data from ``cleaned_magnetic_data.geojson``.
"""

import os

import arviz as az
import geopandas as gpd
import numpy as np
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
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality
from mineye.GeoModel.model_one.probabilistic_model import normalize, _modify_orientations
from mineye.GeoModel.model_one.visualization import gempy_viz, plot_many_observed_vs_forward
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison
from mineye.config import paths
from mineye.config.example_parameters import SoricomModelConfig


nc_path = "arviz_data_magnetic_soricom.nc"


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
        base_resolution=np.array([2, 2, 4]),
    )
    geo_model.interpolation_options.number_octree_levels_surface = 5

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
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir=geophysical_dir,
        )

        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)
        prior_inference_data: az.InferenceData = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=magnetic_observations_tensor,
            n_samples=100,
            plot_trace=True,
        )

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

        # Compute baseline forward magnetics
        baseline_fw_magnetics_np = self._baseline_magnetics(geo_model)

        norm_params = normalize(
            baseline_fw_gravity_np=baseline_fw_magnetics_np,
            observed_gravity=observed_magnetics_nt,
            method="align_to_reference",
            extrapolation_buffer=0.3,
        )

        # 3) Priors
        n_orientations = geo_model.orientations_copy.xyz.shape[0]
        model_priors = {
            self.prior_key_dips: dist.Normal(
                loc=torch.full((n_orientations,), 10.0, dtype=torch.float64),
                scale=torch.tensor(10.0, dtype=torch.float64),
                validate_args=True,
            ),
            self.prior_key_susceptibility: dist.Normal(
                loc=torch.tensor([0.05, 0.001], dtype=torch.float64),
                scale=torch.tensor(0.01, dtype=torch.float64),
            ).to_event(1),
        }

        # 4) Deterministics
        post_forward_dets = {
            "magnetic_response_raw": lambda samples, gm, sol: sol.magnetics,
            "magnetic_response": lambda samples, gm, sol: align_forward_to_observed(
                sol.magnetics, norm_params,
            ),
            "mean_magnetics": lambda samples, gm, sol: torch.mean(
                align_forward_to_observed(sol.magnetics, norm_params),
            ),
            "max_magnetics": lambda samples, gm, sol: torch.max(
                align_forward_to_observed(sol.magnetics, norm_params), 0,
            ),
        }

        # 5) Likelihood
        likelihood_fn = self._generate_multimagnetic_likelihood_hierarchical_per_station(
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
        magnetic_data = gpd.read_file(
            os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson'),
        )
        sampled_magnetic = magnetic_data.sample(n=20, random_state=42)
        observed_magnetics = sampled_magnetic['MAG'].values
        return sampled_magnetic, observed_magnetics

    @staticmethod
    def _setup_magnetic_geomodel(magnetic_data):
        geo_model = _create_soricom_geomodel()

        # Reproject to Soricom UTM if needed (the magnetic data is in EPSG:3857)
        if magnetic_data.crs is not None and magnetic_data.crs.to_string() != 'EPSG:3857':
            magnetic_data = magnetic_data.to_crs('EPSG:3857')

        xy_ravel = np.column_stack([
            np.array(magnetic_data.geometry.x.values),
            np.array(magnetic_data.geometry.y.values),
            np.full(len(magnetic_data), 500),
        ])

        gp.set_centered_grid(
            grid=geo_model.grid,
            centers=xy_ravel,
            resolution=np.array([10, 10, 15]),
            radius=np.array([2000, 5000, 2000]),
        )

        gradient_tensor_dict = calculate_magnetic_gradient_tensor(
            centered_grid=geo_model.grid.centered_grid,
            igrf_params={
                "inclination": 60.79,
                "declination": 1.26,
                "intensity": 47258.4,
            },
            compute_tmi=True,
            units_nT=True,
        )

        geo_model.geophysics_input = gp.data.GeophysicsInput(
            magnetics_input=MagneticsInput(
                mag_kernel=gradient_tensor_dict['tmi_kernel'],
                susceptibilities=np.array([0.05, 0.001]),
                igrf_params={
                    "inclination": gradient_tensor_dict['inclination'],
                    "declination": gradient_tensor_dict['declination'],
                    "intensity": gradient_tensor_dict['intensity'],
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

    @staticmethod
    def _generate_multimagnetic_likelihood_hierarchical_per_station(norm_params):
        def likelihood_fn(solutions: gp.data.Solutions) -> dist.Distribution:
            simulated_magnetics = align_forward_to_observed(
                solutions.magnetics, norm_params,
            )
            import pyro
            pyro.deterministic(r'$\mu_{magnetics}$', simulated_magnetics)
            n_stations = simulated_magnetics.shape[0]

            mu_log_sigma = pyro.sample(
                "mu_log_sigma",
                dist.Normal(
                    torch.tensor(np.log(50.0), dtype=torch.float64),
                    torch.tensor(0.5, dtype=torch.float64),
                ),
            )

            tau_log_sigma = pyro.sample(
                "tau_log_sigma",
                dist.HalfNormal(torch.tensor(0.5, dtype=torch.float64)),
            )

            log_sigma_stations = pyro.sample(
                "log_sigma_stations",
                dist.Normal(
                    mu_log_sigma.expand([n_stations]),
                    tau_log_sigma,
                ).to_event(1),
            )

            sigma_stations = torch.exp(log_sigma_stations)
            pyro.deterministic("sigma_stations", sigma_stations)

            return dist.Normal(simulated_magnetics, sigma_stations).to_event(1)

        return likelihood_fn

    # --- Analysis tests ---

    def test_run_predictive_analysis(self, geophysical_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), nc_path))
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data)

        observed_norm = observed_magnetics_nt
        forward_norm = data.prior[r'magnetic_response'].mean(axis=1)
        many_forward_norm = data.posterior_predictive[r'magnetic_response'].values[
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

        plot_geophysics_comparison(
            forward_norm=data.prior[r'magnetic_response'].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,
            xy_ravel=xy_ravel,
        )

        plot_geophysics_comparison(
            forward_norm=data.posterior_predictive[r'magnetic_response'].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,
            xy_ravel=xy_ravel,
        )

        if hasattr(data, 'prior') and r'magnetic_response' in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[r'magnetic_response'].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel,
            )

        if hasattr(data, 'posterior') and r'magnetic_response' in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.posterior_predictive[r'magnetic_response'].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel,
            )

        gempy_viz(geo_model, data)