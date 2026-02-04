import os
import time
import pytest
import numpy as np
import torch
import arviz as az
import matplotlib.pyplot as plt
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from gempy_probability.modules.plot.plot_posterior import default_red, default_blue

# noinspection PyUnusedImports
from tests.tests_inversions.conftest import simple_geo_model, topography_dir, base_dir
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations
from mineye.GeoModel.model_one.visualization import (gempy_viz,
                                                     plot_many_observed_vs_forward,
                                                     compute_probability_density_fields)

from mineye.GeoModel.model_one.probabilistic_model_likelihoods import enmap_likelihood_fn


def set_priors_enmap(samples, geo_model):
    """A version of set_priors that only modifies orientations for EnMap."""
    return _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key="dips"
    )


class TestEnMapInversion:
    def test_run_predictive(self, simple_geo_model, base_dir, n_samples=50):
        geo_model, labels_enmap_tensor, prob_model = self._create_probabilistic_model(base_dir, simple_geo_model)

        # * 7) Run predictive
        print("Running prior predictive...")
        prior_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            # y_obs_list=labels_enmap_tensor,
            y_obs_list=None,
            n_samples=n_samples,
            plot_trace=True
        )
        # Extract predicted labels (shape: n_samples x n_points)
        # Check the exact key name in your prior_data dict (usually 'obs' or the name you gave it)
        pred_labels = prior_data.prior['EnMap Labels'][0]
        obs_labels = labels_enmap_tensor  # Your ground truth

        plt.figure(figsize=(10, 5))
        # Flatten pred_labels since it's (samples, labels) - we want all predictions in one histogram
        plt.hist([pred_labels.values.flatten(), obs_labels.numpy()],
                 label=['Prior Predicted', 'Observed Ground Truth'],
                 bins=np.arange(4) - 0.5, density=True, rwidth=0.8)
        plt.xticks([0, 1])
        plt.xlabel("Rock Type Label")
        plt.ylabel("Density")
        plt.title("Are my priors exploring all rock types?")
        plt.legend()
        plt.show()

        def plot_probability_heatmap(data, group='prior'):
            import seaborn as sns
            # Get the data array from ArviZ InferenceData
            # Standard Shape: (chain, draw, n_points, n_classes)
            # Example: (1, 50, 57, 3)
            probs_da = data[group]['probs_pred']

            # 1. Average over 'chain' and 'draw' dimensions to get mean per point
            # Result Shape: (n_points, n_classes)
            mean_probs = probs_da.mean(dim=["chain", "draw"]).values

            # 2. Transpose for the Heatmap
            # We want Y-axis = Classes, X-axis = Points
            # Result Shape: (n_classes, n_points)
            heatmap_data = mean_probs.T

            # 3. Plot
            plt.figure(figsize=(14, 4))
            sns.heatmap(heatmap_data, cmap="Blues", vmin=0, vmax=1,
                        annot=False, # Set True if you want numbers inside cells
                        cbar_kws={'label': 'Probability'})

            plt.title(f"Average Class Probabilities per Point ({group.capitalize()})")
            plt.xlabel("Point Index (0 to 56)")
            plt.ylabel("Rock Unit (Class)")

            # Fix Y-axis labels to be integers (0, 1, 2)
            plt.yticks(ticks=np.arange(heatmap_data.shape[0]) + 0.5,
                       labels=np.arange(heatmap_data.shape[0]),
                       rotation=0)

            # plt.tight_layout()
            plt.show()

        plot_probability_heatmap(prior_data, 'prior')


    def test_enmap_inversion(self, simple_geo_model, base_dir, n_samples=50,
                             arviz_data_filename="arviz_data_enmap_Feb04_2026.nc"):
        """Test EnMap inversion following the gravity inversion structure."""
        geo_model, labels_enmap_tensor, prob_model = self._create_probabilistic_model(base_dir, simple_geo_model)

        # 7. Run predictive
        print("Running prior predictive...")
        prior_inference_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=None,
            n_samples=n_samples,
            plot_trace=True
        )

        # 8. Run inference (NUTS)
        print("Running NUTS inference...")
        # Note: NUTS might struggle with discrete labels due to lack of gradients.
        # We decrease the number of samples for the test
        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=labels_enmap_tensor,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                # num_samples=20,
                # warmup_steps=20,
                num_samples=200,
                warmup_steps=200,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        data.extend(prior_inference_data)
        data.to_netcdf(os.path.join(os.path.dirname(__file__), arviz_data_filename))

    @staticmethod
    def _create_probabilistic_model(base_dir, simple_geo_model):
        # 1. Load EnMap extracted points and labels
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')

        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")

        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        labels_enmap[labels_enmap == 2] = 1 # * Normalize the labels

        print(f"\nLoaded {len(xyz_enmap)} points from EnMap extraction.")

        # 2. Set custom grid in GemPy model
        simple_geo_model.interpolation_options.mesh_extraction = False
        simple_geo_model.interpolation_options.evaluation_options.number_octree_levels = 1
        gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
        gp.set_active_grid(
            grid=simple_geo_model.grid,
            grid_type=[simple_geo_model.grid.GridTypes.CUSTOM],
            reset=True
        )

        # 3. Define Priors
        model_priors = {
                'dips': dist.Normal(
                    loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                )
        }
        # 5. Set up Pyro model
        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=set_priors_enmap,
            likelihood_fn=enmap_likelihood_fn,
            obs_name="EnMap Labels"
        )

        # 6. Prepare observed data
        labels_enmap_tensor = torch.tensor(labels_enmap, dtype=torch.float64)

        return simple_geo_model, labels_enmap_tensor, prob_model

    def test_run_predictive_analysis(self, simple_geo_model, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_enmap.nc"))
        az.plot_trace(data.prior)
        plt.show()

        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')
        labels_enmap = np.load(labels_path)

        # Prepare data
        observed_norm = labels_enmap
        forward_norm = data.prior["EnMap Labels"].mean(axis=1)

        # We take some samples from the prior
        many_forward_norm = data.prior["EnMap Labels"].values[0, -20:]

        plot_many_observed_vs_forward(forward_norm, many_forward_norm, observed_norm)

    def test_run_kde_sections(self, simple_geo_model, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_enmap.nc"))
        gempy_viz(simple_geo_model, data, n_samples=100, ve=3)

    def test_run_analysis(self, simple_geo_model, base_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_enmap.nc"))

        az.plot_posterior(data, var_names=["dips"])
        plt.show()

        plt.rcParams['figure.dpi'] = 72  # Lower DPI for faster rendering

        axes = az.plot_density(
            data=[data, data.prior],
            var_names=["dips"],
            filter_vars="like",
            hdi_prob=0.9999,
            shade=.2,
            data_labels=["Posterior", "Prior"],
            colors=[default_red, default_blue],
        )
        plt.show()

        # * Analysis Gempy Model
        gempy_viz(simple_geo_model, data)

    def test_probability_plots(self, simple_geo_model, base_dir, topography_dir):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data_enmap.nc"))

        self._probability_fields_for(simple_geo_model, data.prior, topography_dir)
        self._probability_fields_for(simple_geo_model, data.posterior)

    @staticmethod
    def _probability_fields_for(geo_model, inference_data, topography_dir=None):
        online_prob = compute_probability_density_fields(
            geo_model,
            inference_data,
            n_samples=20
        )
        import gempy_viewer as gpv

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.probability_field[0],
            show_data=True,
            ve=5,
            kwargs_lithology={
                    'cmap': 'viridis',
                    'norm': None
            }
        )

        gpv.plot_2d(
            geo_model,
            override_regular_grid=online_prob.entropy,
            show_data=True,
            ve=5,
            kwargs_lithology={
                    'cmap': 'magma',
                    'norm': None
            }
        )

        if topography_dir:
            import gempy as gp
            topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')
            gp.set_topography_from_file(
                grid=geo_model.grid,
                filepath=topography_path,
                crop_to_extent=[
                        geo_model.grid.extent[0],
                        geo_model.grid.extent[2],
                        geo_model.grid.extent[1],
                        geo_model.grid.extent[3]
                ]
            )

            gp.compute_model(geo_model)

            # * We inject the entropy field into the model
            geo_model.solutions.raw_arrays.scalar_field_matrix[0] = online_prob.entropy
            p3d = gpv.plot_3d(
                model=geo_model,
                active_scalar_field="sf_0",
                show_scalar=True,
                show_lith=False,
                show_topography=True,
                image=False,
                ve=4,
                threshold_kwargs={'value': [0.1, 0.9], 'invert': False},
                kwargs_pyvista_bounds={
                        'show_xlabels': False,
                        'show_ylabels': False,
                        'show_zlabels': False,
                }
            )
