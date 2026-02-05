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


def generate_joint_likelihood(norm_params):
    def joint_likelihood_fn(solutions: gp.data.Solutions):
        gravity_dist = generate_multigravity_likelihood_hierarchical_per_station(norm_params)
        enmap_dist = enmap_likelihood_fn(solutions)

        return [gravity_dist, enmap_dist]

    return joint_likelihood_fn
