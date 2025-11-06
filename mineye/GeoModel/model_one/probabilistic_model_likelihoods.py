from functools import partial

import numpy as np
import pyro
import torch
from pyro import distributions as dist

import gempy as gp
from mineye.GeoModel.geophysics import align_forward_to_observed


def generate_multigravity_likelihood_diagonal(norm_params):
    """
    Use independent Normal distributions instead of multivariate.
    Much more stable and faster.
    """

    def likelihood_fn(solutions: gp.data.Solutions) -> dist.Distribution:
        simulated_geophysics = align_forward_to_observed(-solutions.gravity, norm_params)
        pyro.deterministic(r'$\mu_{gravity}$', simulated_geophysics)
        n_stations = simulated_geophysics.shape[0]
        # ? Half normal is too unstable. Try with some prior that is somehow more stable
        
        # Sample noise standard deviation
        # sigma = pyro.sample(
        #     "sigma",
        #     dist.HalfNormal(torch.tensor(5_000.0, dtype=torch.float64)).expand([n_stations]).to_event(1)  # 100 µGal noise
        # )

        # sigma = torch.tensor(2000.0, dtype=torch.float64)
        
        # Independent Normal likelihood (much more stable!)
        sigma = pyro.sample(
            "sigma",
            dist.LogNormal(
                torch.tensor(np.log(1000.0), dtype=torch.float64),
                torch.tensor(0.5, dtype=torch.float64)
            )
        )
        
        return dist.Normal(simulated_geophysics, sigma).to_event(1)

    return likelihood_fn


def generate_multigravity_likelihood(covariance_matrix, norm_params):
    return partial(multigravity_likelihood, covariance_matrix=covariance_matrix, norm_params=norm_params)


def generate_multigravity_likelihood_hierarchical(xy_locations: torch.Tensor, norm_params):
    """
    Generate hierarchical likelihood with hyperparameters sampled INSIDE.
    
    This is the correct pattern for Pyro/NUTS.
    """

    def likelihood_fn(solutions: gp.data.Solutions) -> dist.Distribution:
        # Normalize the forward model output
        simulated_geophysics = align_forward_to_observed(-solutions.gravity, norm_params)
        pyro.deterministic(r'$\mu_{gravity}$', simulated_geophysics)

        # ✓ Sample hyperparameters HERE, inside the likelihood
        length_scale = pyro.sample(
            "length_scale",
            dist.LogNormal(
                loc=torch.tensor(np.log(2000.0), dtype=torch.float64),
                scale=torch.tensor(0.8, dtype=torch.float64)
            )
        )

        variance = pyro.sample(
            "variance",
            dist.InverseGamma(
                concentration=torch.tensor(3.0, dtype=torch.float64),
                rate=torch.tensor(75000.0, dtype=torch.float64)
            )
        )

        nu = pyro.sample(
            "nu",
            dist.Exponential(torch.tensor(0.2, dtype=torch.float64))
        ) + 2.0

        # Build covariance matrix with sampled hyperparameters
        covariance_matrix = gaussian_kernel(xy_locations, length_scale, variance)

        # Compute Cholesky
        try:
            scale_tril = torch.linalg.cholesky(covariance_matrix)
        except torch._C._LinAlgError as e:
            print(f"Cholesky failed with length_scale={length_scale.item():.2f}, "
                  f"variance={variance.item():.2f}")
            raise

        # Return Student-t likelihood
        return dist.MultivariateStudentT(
            df=nu,
            loc=simulated_geophysics,
            scale_tril=scale_tril
        )

    return likelihood_fn


def gaussian_kernel(locations, length_scale, variance, nugget=None):
    """
    Numerically stable Gaussian kernel with automatic jitter.
    """
    import torch

    # Type safety
    if not isinstance(locations, torch.Tensor):
        locations = torch.tensor(locations, dtype=torch.float64)
    else:
        locations = locations.to(dtype=torch.float64)

    length_scale = torch.as_tensor(length_scale, dtype=torch.float64)
    variance = torch.as_tensor(variance, dtype=torch.float64)

    # Default nugget: 0.1% of signal variance
    if nugget is None:
        nugget = 0.001 * variance
    else:
        nugget = torch.as_tensor(nugget, dtype=torch.float64)

    n_stations = locations.shape[0]

    # Compute distances
    distance_squared = torch.cdist(locations, locations, p=2).pow(2)

    # Stabilized exponential
    exponent = -0.5 * distance_squared / (length_scale.pow(2) + 1e-10)
    exponent = torch.clamp(exponent, min=-30.0)

    # Kernel
    K = variance * torch.exp(exponent)

    # Add nugget
    K = K + torch.eye(n_stations, dtype=torch.float64, device=K.device) * nugget

    return K
