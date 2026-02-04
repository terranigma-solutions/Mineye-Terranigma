import os
import time
import pytest
import numpy as np
import torch
import arviz as az
from pyro import distributions as dist

import gempy as gp
import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig

# noinspection PyUnusedImports
from tests.tests_inversions.conftest import simple_geo_model, topography_dir, base_dir
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations


def set_priors_enmap(samples, geo_model):
    """A version of set_priors that only modifies orientations for EnMap."""
    return _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key="dips"
    )


class TestEnMapInversion:
    def test_enmap_inversion(self, simple_geo_model, base_dir, n_samples=50):
        """Test EnMap inversion following the gravity inversion structure."""
        
        # 1. Load EnMap extracted points and labels
        xyz_path = os.path.join(base_dir, 'central_xyz.npy')
        labels_path = os.path.join(base_dir, 'central_labels.npy')
        
        if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
            pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")
            
        xyz_enmap = np.load(xyz_path)
        labels_enmap = np.load(labels_path)
        
        print(f"\nLoaded {len(xyz_enmap)} points from EnMap extraction.")
        
        # 2. Set custom grid in GemPy model
        gp.set_custom_grid(simple_geo_model.grid, xyz_enmap)
        
        # 3. Define Priors
        # Using similar priors as in gravity inversion
        model_priors = {
            'dips': dist.Normal(
                loc=(torch.ones(simple_geo_model.orientations_copy.xyz.shape[0]) * 10),
                scale=torch.tensor(10, dtype=torch.float64),
                validate_args=True
            )
        }
        
        # 4. Define EnMap Likelihood
        # The goal is to minimize the mismatch of labels.
        # Since GemPy labels are discrete, we use a Normal distribution with small sigma 
        # as a continuous proxy to follow the structure, noting that gradients will be zero.
        def enmap_likelihood_fn(solutions):
            labels_gempy = solutions.raw_arrays.custom
            
            # Ensure labels_gempy is a torch tensor and has the right type
            if not isinstance(labels_gempy, torch.Tensor):
                labels_gempy = torch.tensor(labels_gempy, dtype=torch.float64)
            else:
                labels_gempy = labels_gempy.to(torch.float64)
                
            return dist.Normal(loc=labels_gempy, scale=0.1).to_event(1)

        # 5. Set up Pyro model
        prob_model = gpp.make_gempy_pyro_model(
            priors=model_priors,
            set_interp_input_fn=set_priors_enmap,
            likelihood_fn=enmap_likelihood_fn,
            obs_name="EnMap Labels"
        )
        
        # 6. Prepare observed data
        labels_enmap_tensor = torch.tensor(labels_enmap, dtype=torch.float64)
        
        # 7. Run predictive
        print("Running prior predictive...")
        prior_data = gpp.run_predictive(
            prob_model=prob_model,
            geo_model=simple_geo_model,
            y_obs_list=labels_enmap_tensor,
            n_samples=n_samples,
            plot_trace=False
        )
        
        # 8. Run inference (NUTS)
        print("Running NUTS inference...")
        # Note: NUTS might struggle with discrete labels due to lack of gradients.
        # This follows the requested structure.
        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=simple_geo_model,
            y_obs_list=labels_enmap_tensor,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                num_samples=20,
                warmup_steps=20,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )
        
        # 9. Save or check results
        assert data is not None
        print("\nInference completed successfully.")
        
        if hasattr(data, 'posterior'):
            print(f"Posterior samples for dips: {data.posterior['dips'].shape}")

