import os
import torch
import pyro
from pyro import distributions as dist
import gempy as gp
from mineye.GeoModel.model_one.probabilistic_model import _modify_orientations, _modify_densities
from mineye.GeoModel.model_one.probabilistic_model_likelihoods import (
    generate_multigravity_likelihood_hierarchical_per_station,
    enmap_likelihood_fn
)


def joint_set_priors(samples, geo_model):
    """Combined version of set_priors that handles both orientations and densities."""
    # This is essentially what the original set_priors does, 
    # but we make sure it's compatible with both models.
    interpolation_input = _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key="dips"
    )

    if "density" in samples:
        _modify_densities(
            samples=samples,
            geo_model=geo_model,
            key="density"
        )

    return interpolation_input


def generate_joint_likelihood(norm_params, n_gravity_points: int):
    """
    Factory for a joint likelihood function that combines gravity and EnMap.
    """
    gravity_likelihood_gen = generate_multigravity_likelihood_hierarchical_per_station(
        norm_params=norm_params
    )

    def joint_likelihood_fn(solutions: gp.data.Solutions):
        # 1. Slice gravity solutions
        # solutions.gravity usually matches the custom grid if it's the only one.
        # But here solutions.gravity will have size (n_gravity_points + n_enmap_points).
        # We only want the first n_gravity_points for the gravity likelihood.

        # Backup original gravity
        full_gravity = solutions.gravity
        solutions.gravity = full_gravity[:n_gravity_points]

        # 2. Gravity Likelihood
        gravity_dist = gravity_likelihood_gen(solutions)

        # Restore full gravity just in case (though it might not be needed)
        solutions.gravity = full_gravity

        # 3. EnMap Likelihood
        # enmap_likelihood_fn uses solutions.octrees_output[0].last_output_center.exported_fields.scalar_field
        # with output_center.grid.custom_grid_slice.
        # This slice should correctly point to all custom points.
        # However, EnMap data starts AFTER gravity data in our combined custom grid.

        # We might need to modify enmap_likelihood_fn to take a slice or 
        # modify the solutions object temporarily.

        enmap_dist = enmap_likelihood_fn(solutions, slice_idx=slice(n_gravity_points, None))
        # Note: enmap_dist is Categorical(probs=probs).to_event(1)
        # We need to make sure 'probs' only covers the EnMap points.

        return [gravity_dist, enmap_dist]

    return joint_likelihood_fn
