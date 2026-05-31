from functools import partial
from typing import Union, Sequence, Tuple

import torch
from pyro.distributions import Distribution

import gempy as gp
from gempy_engine.core.backend_tensor import BackendTensor
from gempy_engine.core.data.interpolation_input import InterpolationInput


def normalize(baseline_fw_gravity_np, observed_gravity, method="quantile_align", extrapolation_buffer=0.3):
    # Convert observed gravity from mGal to μGal for comparison
    from mineye.GeoModel.geophysics import compute_alignment_params
    norm_params = compute_alignment_params(
        observed=observed_gravity,
        baseline_forward=baseline_fw_gravity_np,
        method=method,
        extrapolation_buffer=extrapolation_buffer
    )

    return norm_params


def set_priors(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
)-> InterpolationInput:
    interpolation_input = _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key=r"dips"
    )

    _modify_densities(
        samples=samples,
        geo_model=geo_model,
        key=r"density"
    )
    
    return interpolation_input


def set_magnetic_priors(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
) -> InterpolationInput:
    """Set priors for magnetic inversion - modifies susceptibilities and orientations."""
    interpolation_input = _modify_orientations(
        samples=samples,
        geo_model=geo_model,
        key=r"dips"
    )

    if "susceptibility" in samples:
        susceptibilities = samples["susceptibility"]
        if geo_model.geophysics_input and geo_model.geophysics_input.magnetics_input:
            geo_model.geophysics_input.magnetics_input.susceptibilities = susceptibilities

    return interpolation_input


def _modify_densities(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
        key: str
) -> None:
    """Modify densities for different keys."""
    geo_model.geophysics_input.densities = samples[key]


def create_orientation_modifier(key: str):
    """Factory function that creates orientation modifier functions for different keys."""
    return partial(_modify_orientations, key=key)


def _modify_orientations(
        samples: dict[str, Distribution],
        geo_model: gp.data.GeoModel,
        key: str
) -> InterpolationInput:
    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame

    interp_input: InterpolationInput = interpolation_input_from_structural_frame(geo_model)
    samples_value = samples[key]

    og_gradients = interp_input.orientations.dip_gradients
    azimuth, dip, polarity = compute_adp_from_gradients(
        G_x=og_gradients[:, 0],
        G_y=og_gradients[:, 1],
        G_z=og_gradients[:, 2],
    )

    gradients = convert_orientation_to_pole_vector(
        azimuth=azimuth,
        dip=samples_value,
        polarity=polarity
    )

    interp_input.orientations.dip_gradients = gradients
    return interp_input


def convert_orientation_to_pole_vector(
        azimuth: Union[torch.Tensor, Sequence[float]],
        dip: Union[torch.Tensor, Sequence[float]],
        polarity: Union[torch.Tensor, Sequence[float]]
) -> torch.Tensor:
    """
    Convert orientation parameters (azimuth, dip, polarity) to pole vectors (gradients).
    
    PyTorch version that preserves gradient flow for automatic differentiation.
    
    Parameters
    ----------
    azimuth : torch.Tensor or Sequence[float]
        Azimuth angles in degrees [0, 360]
    dip : torch.Tensor or Sequence[float]
        Dip angles in degrees [0, 180]
    polarity : torch.Tensor or Sequence[float]
        Polarity values (typically ±1)
    
    Returns
    -------
    gradients : torch.Tensor
        Pole vectors with shape (n, 3) where columns are [G_x, G_y, G_z]
    
    Notes
    -----
    - Converts inputs to tensors if they aren't already
    - All trigonometric operations preserve gradients
    - Output shape is (n_orientations, 3)
    
    Example
    -------
    >>> azimuth = torch.tensor([90.0, 180.0], dtype=torch.float64)
    >>> dip = torch.tensor([45.0, 30.0], dtype=torch.float64)
    >>> polarity = torch.tensor([1.0, 1.0], dtype=torch.float64)
    >>> gradients = convert_orientation_to_pole_vector(azimuth, dip, polarity)
    >>> gradients.shape
    torch.Size([2, 3])
    """
    # Convert to tensors if needed (preserves gradients if already tensors)
    float_ = BackendTensor.dtype_obj
    if not isinstance(azimuth, torch.Tensor):
        azimuth = torch.as_tensor(azimuth, dtype=float_)
    if not isinstance(dip, torch.Tensor):
        dip = torch.as_tensor(dip, dtype=float_)
    if not isinstance(polarity, torch.Tensor):
        polarity = torch.as_tensor(polarity, dtype=float_)

    # Ensure consistent dtype
    azimuth = azimuth.to(dtype=float_)
    dip = dip.to(dtype=float_)
    polarity = polarity.to(dtype=float_)

    # Convert degrees to radians
    azimuth_rad = torch.deg2rad(azimuth)
    dip_rad = torch.deg2rad(dip)

    # Calculate gradient components (all differentiable operations)
    G_x = torch.sin(dip_rad) * torch.sin(azimuth_rad) * polarity
    G_y = torch.sin(dip_rad) * torch.cos(azimuth_rad) * polarity
    G_z = torch.cos(dip_rad) * polarity

    # Stack into (n, 3) array
    # Use torch.stack instead of vstack for better control
    gradients = torch.stack([G_x, G_y, G_z], dim=-1)

    return gradients


def compute_adp_from_gradients(
        G_x: torch.Tensor,
        G_y: torch.Tensor,
        G_z: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """
    Compute azimuth, dip, and polarity from gradient components.
    
    PyTorch version that preserves gradient flow for automatic differentiation.
    
    Parameters
    ----------
    G_x : torch.Tensor
        Gradient in X direction
    G_y : torch.Tensor
        Gradient in Y direction
    G_z : torch.Tensor
        Gradient in Z direction
    
    Returns
    -------
    azimuth : torch.Tensor
        Azimuth in degrees [0, 360]
    dip : torch.Tensor
        Dip angle in degrees [0, 180]
    polarity : torch.Tensor
        Polarity values (all ones in this implementation)
    
    Notes
    -----
    - Uses torch.where() instead of boolean indexing to preserve gradients
    - Clamps values for numerical stability
    - All operations are differentiable
    """
    # Calculate polarity (assumed to be 1 for all)
    from gempy_engine.core.backend_tensor import BackendTensor
    polarity = torch.ones_like(G_x, dtype=BackendTensor.dtype_obj)

    # Calculate dip
    # Clamp G_z/polarity to [-1, 1] for numerical stability in arccos
    cos_dip = torch.clamp(G_z / polarity, min=-1.0, max=1.0)
    dip = torch.rad2deg(torch.acos(cos_dip))

    # Replace NaN with 0 (torch equivalent of np.nan_to_num)
    dip = torch.where(torch.isnan(dip), torch.zeros_like(dip, dtype=BackendTensor.dtype_obj), dip)

    # Calculate azimuth
    azimuth = torch.rad2deg(torch.atan2(G_x / polarity, G_y / polarity))

    # Replace NaN with 0
    azimuth = torch.where(torch.isnan(azimuth), torch.zeros_like(azimuth, dtype=BackendTensor.dtype_obj), azimuth)

    # Shift values from [-180, 0] to [180, 360]
    # Use torch.where to preserve gradients (instead of boolean indexing)
    azimuth = torch.where(azimuth < 0, azimuth + 360.0, azimuth)

    # Adjust azimuth where dip is nearly zero (azimuth is undefined)
    azimuth = torch.where(dip < 0.001, torch.zeros_like(azimuth, dtype=BackendTensor.dtype_obj), azimuth)

    return azimuth, dip, polarity
