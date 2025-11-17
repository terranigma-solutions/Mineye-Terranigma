import os

import arviz as az
import geopandas as gpd
import numpy as np
import torch
from matplotlib import pyplot as plt
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_one.inference_diagnostics import check_mcmc_quality
from mineye.GeoModel.model_one.probabilistic_model import normalize, _modify_orientations
from mineye.GeoModel.model_one.visualization import gempy_viz, plot_many_observed_vs_forward
from mineye.GeoModel.plotting.probabilistic_analysis import plot_geophysics_comparison


# noinspection PyUnusedImport


class TestMagneticInversion:
    prior_key_dips = r'dips'
    prior_key_susceptibility = r'susceptibility'

    def test_run_predictive(self, simple_geo_model, geophysical_dir, n_samples=50):
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir=geophysical_dir,
            simple_geo_model=simple_geo_model
        )

        # * 7) Run predictive
        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=magnetic_observations_tensor,
                n_samples=100,
                plot_trace=True
            )

    def test_magnetic_inversion(
            self, simple_geo_model, geophysical_dir, n_samples=50,
            arviz_data_filename="arviz_data_magnetic_Nov10_I_hierarchical.nc"
    ):
        """Test magnetic inversion using NUTS sampler."""
        print("Test magnetic inversion...")
        # Use actual magnetic measurement device locations
        geo_model, observed_magnetics_nt, prob_model = self._create_probabilistic_model(
            geophysical_dir, simple_geo_model
        )

        # * 7) Run predictive
        magnetic_observations_tensor = torch.tensor(observed_magnetics_nt)
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=magnetic_observations_tensor,
                n_samples=100,
                plot_trace=True
            )

        # * 8) Run inference
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
                num_chains=1
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        if compute_prior_predictive:
            data.extend(prior_inference_data)

        data.to_netcdf(os.path.join(os.path.dirname(__file__), arviz_data_filename))

    @staticmethod
    def _create_probabilistic_model(geophysical_dir, simple_geo_model):
        # * 1) Read magnetic data
        magnetic_data, observed_magnetics_nt = TestMagneticInversion._read_magnetics(geophysical_dir)

        # * 2) Setup initial Geomodel and normalize forward magnetics to observed magnetics
        geo_model, xy_ravel = TestMagneticInversion._setup_magnetic_geomodel(
            magnetic_data, simple_geo_model
        )
        geo_model.interpolation_options.sigmoid_slope = 100

        # Compute baseline forward model
        baseline_fw_magnetics_np = TestMagneticInversion._baseline_magnetics(geo_model)

        norm_params = normalize(
            baseline_fw_gravity_np=baseline_fw_magnetics_np,  # Reusing gravity normalization
            observed_gravity=observed_magnetics_nt,
            method="align_to_reference",
            extrapolation_buffer=0.3
        )

        # * 3) Setup Priors
        model_priors = {
                TestMagneticInversion.prior_key_dips          : dist.Normal(
                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # Dip 10 degrees
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                ),
                TestMagneticInversion.prior_key_susceptibility: dist.Normal(
                    loc=(torch.tensor([
                            0.05,  # plutonites - typical magnetic susceptibility (SI units)
                            0.001  # sedimentary host - low susceptibility
                    ])),
                    scale=torch.tensor(0.01),  # Reasonable uncertainty in susceptibility
                ).to_event(1)
        }

        # * 4) Set up Deterministics
        post_forward_dets = {
                "magnetic_response_raw": lambda samples, gm, sol: sol.magnetics,  # Store raw magnetics
                "magnetic_response"    : lambda samples, gm, sol: align_forward_to_observed(
                    sol.magnetics, norm_params
                ),  # Normalized!
                "mean_magnetics"       : lambda samples, gm, sol: torch.mean(
                    align_forward_to_observed(sol.magnetics, norm_params)
                ),
                "max_magnetics"        : lambda samples, gm, sol: torch.max(
                    align_forward_to_observed(sol.magnetics, norm_params), 0
                ),
        }

        # * 5) Set up likelihood functions
        likelihood_fn = TestMagneticInversion._generate_multimagnetic_likelihood_hierarchical_per_station(
            norm_params=norm_params
        )

        # * 6) Set up Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=TestMagneticInversion._set_magnetic_priors,
            likelihood_fn=likelihood_fn,
            pre_forward_deterministics={},
            post_forward_deterministics=post_forward_dets,
            obs_name="Magnetic Measurement"
        )
        return geo_model, observed_magnetics_nt, prob_model

    @staticmethod
    def _read_magnetics(geophysical_dir):
        """Read magnetic data from geophysical directory."""
        magnetic_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_magnetic_data.geojson'))
        # Sample 20 points to match gravity setup
        sampled_magnetic = magnetic_data.sample(n=20, random_state=42)
        observed_magnetics = sampled_magnetic['MAG'].values  # in nT
        return sampled_magnetic, observed_magnetics

    @staticmethod
    def _setup_magnetic_geomodel(magnetic_data, simple_geo_model):
        """Setup geomodel with magnetic measurement locations."""
        xy_ravel = np.column_stack([
                np.array(magnetic_data.geometry.x.values),
                np.array(magnetic_data.geometry.y.values),
                np.full(len(magnetic_data), 0)  # Set Z to surface elevation
        ])

        # Setup centered grid for magnetics
        gp.set_centered_grid(
            grid=simple_geo_model.grid,
            centers=xy_ravel,
            resolution=np.array([10, 10, 15]),
            radius=np.array([5000, 5000, 5000])
        )

        # Calculate magnetic gradient tensor
        from gempy_engine.modules.geophysics.magnetic_gradient import calculate_magnetic_gradient_tensor
        from gempy_engine.core.data.geophysics_input import MagneticsInput

        gradient_tensor_dict = calculate_magnetic_gradient_tensor(
            centered_grid=simple_geo_model.grid.centered_grid,
            igrf_params={
                    "inclination": 60.79,  # Huelva, Spain (2025)
                    "declination": 1.26,  # Huelva, Spain (2025)
                    "intensity"  : 47258.4  # Earth's field strength in nT
            },
            compute_tmi=True,
            units_nT=True,
        )

        # Configure geophysics input
        simple_geo_model.geophysics_input = gp.data.GeophysicsInput(
            magnetics_input=MagneticsInput(
                mag_kernel=gradient_tensor_dict['tmi_kernel'],
                susceptibilities=np.array([0.05, 0.001]),  # Initial values
                igrf_params={
                        "inclination": gradient_tensor_dict['inclination'],
                        "declination": gradient_tensor_dict['declination'],
                        "intensity"  : gradient_tensor_dict['intensity']
                }
            )
        )

        simple_geo_model.interpolation_options.mesh_extraction = False

        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CENTERED],
            reset=True
        )
        gp.compute_model(simple_geo_model)

        return simple_geo_model, xy_ravel

    @staticmethod
    def _baseline_magnetics(geo_model):
        """Compute baseline forward magnetics model."""
        sol = gp.compute_model(geo_model)
        return sol.magnetics.detach().cpu().numpy() if hasattr(sol.magnetics, 'detach') else sol.magnetics

    @staticmethod
    def _set_magnetic_priors(samples: dict, geo_model: gp.data.GeoModel):
        """Set priors for magnetic inversion - modifies susceptibilities and orientations."""
        # Modify orientations (dips)

        interpolation_input = _modify_orientations(
            samples=samples,
            geo_model=geo_model,
            key=r"dips"
        )

        # Modify susceptibilities
        if TestMagneticInversion.prior_key_susceptibility in samples:
            susceptibilities = samples[TestMagneticInversion.prior_key_susceptibility]
            if geo_model.geophysics_input and geo_model.geophysics_input.magnetics_input:
                geo_model.geophysics_input.magnetics_input.susceptibilities = susceptibilities

        return interpolation_input

    @staticmethod
    def _generate_multimagnetic_likelihood_hierarchical_per_station(norm_params):
        """
        Per-station noise with hierarchical structure for magnetics.
        Adapted from gravity likelihood.
        """

        def likelihood_fn(solutions: gp.data.Solutions) -> dist.Distribution:
            simulated_magnetics = align_forward_to_observed(solutions.magnetics, norm_params)
            import pyro
            pyro.deterministic(r'$\mu_{magnetics}$', simulated_magnetics)
            n_stations = simulated_magnetics.shape[0]

            # Global hyperprior on typical noise level (magnetics typically has higher noise)
            mu_log_sigma = pyro.sample(
                "mu_log_sigma",
                dist.Normal(
                    torch.tensor(np.log(50.0), dtype=torch.float64),  # ~50 nT typical noise
                    torch.tensor(0.5, dtype=torch.float64)
                )
            )

            # How much stations vary from each other
            tau_log_sigma = pyro.sample(
                "tau_log_sigma",
                dist.HalfNormal(torch.tensor(0.5, dtype=torch.float64))
            )

            # Per-station noise
            log_sigma_stations = pyro.sample(
                "log_sigma_stations",
                dist.Normal(
                    mu_log_sigma.expand([n_stations]),
                    tau_log_sigma
                ).to_event(1)
            )

            sigma_stations = torch.exp(log_sigma_stations)
            pyro.deterministic("sigma_stations", sigma_stations)

            return dist.Normal(simulated_magnetics, sigma_stations).to_event(1)

        return likelihood_fn

    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(
            os.path.dirname(__file__), "arviz_data_magnetic.nc"
        ))
        # Run comprehensive diagnostics with plots
        check_mcmc_quality(data, verbose=True, plot=True)

    def test_run_predictive_analysis(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(
            os.path.dirname(__file__), "arviz_data_magnetic_Nov10_I_hierarchical.nc"
        ))
        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data, simple_geo_model)

        # Prepare data
        observed_norm = observed_magnetics_nt
        forward_norm = data.prior[r'magnetic_response'].mean(axis=1)
        if PRIOR := True:
            many_forward_norm = data.prior[r'magnetic_response'].values[0, -20:]
        else:
            many_forward_norm = data.posterior_predictive[r'magnetic_response'].values[0, -40:-20]

        plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)

    def test_run_kde_sections(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(
            os.path.dirname(__file__), "arviz_data_magnetic_Nov10_I_hierarchical.nc"
        ))

        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data, simple_geo_model)

        gempy_viz(geo_model, data, n_samples=100)

    def test_run_outlier_detection(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(
            os.path.dirname(__file__), "arviz_data_magnetic_Nov10_I_hierarchical.nc"
        ))

        posterior_sigmas = data.posterior_predictive["sigma_stations"].values  # shape: (chains, samples, 20)

        # Stations with consistently high σ are either:
        # - Actually noisier (geology, instrument, location)
        # - Outliers / model misspecification
        # - Areas where your forward model is wrong

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

    def test_run_analysis(self, simple_geo_model, geophysical_dir):
        data = az.from_netcdf(os.path.join(
            os.path.dirname(__file__), "arviz_data_magnetic_Nov10_I_hierarchical.nc"
        ))
        data.posterior

        az.plot_posterior(data, var_names=["dips", "susceptibility"])
        plt.show()

        magnetic_data, observed_magnetics_nt = self._read_magnetics(geophysical_dir)
        geo_model, xy_ravel = self._setup_magnetic_geomodel(magnetic_data, simple_geo_model)

        plt.rcParams['figure.dpi'] = 72  # Lower DPI for faster rendering

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

        # Plot comparisons
        plot_geophysics_comparison(
            forward_norm=data.prior[r'magnetic_response'].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,  # Note: function name references gravity but works for magnetics
            xy_ravel=xy_ravel
        )

        plot_geophysics_comparison(
            forward_norm=data.posterior_predictive[r'magnetic_response'].mean(axis=1),
            normalization_method='align_to_reference',
            observed_ugal=observed_magnetics_nt,
            xy_ravel=xy_ravel
        )

        # * 9) Analysis inference
        if hasattr(data, 'prior') and r'magnetic_response' in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            # Reuse gravity uncertainty plotting for magnetics (just different units)
            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.prior[r'magnetic_response'].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel
            )

        if hasattr(data, 'posterior') and r'magnetic_response' in data.prior:
            from mineye.GeoModel.model_one.visualization import generate_gravity_uncertainty_plots

            magnetic_samples_norm, unit_label = generate_gravity_uncertainty_plots(
                gravity_samples_norm=data.posterior_predictive[r'magnetic_response'].values[0, :],
                observed_gravity_ugal=observed_magnetics_nt,
                xy_ravel=xy_ravel
            )

        # * 9) Analysis Gempy Model
        gempy_viz(geo_model, data)
