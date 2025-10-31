import numpy as np
from typing import Literal, Tuple, Dict, Optional


def compute_normalization_params(
        reference_data: np.ndarray,
        method: Literal['zscore', 'minmax', 'mean_center', 'relative', 'robust_zscore', 'align_to_reference'] = 'zscore',
        baseline_forward_model: Optional[np.ndarray] = None,
        verbose: bool = True
) -> Dict[str, float]:
    """
    Compute normalization parameters from reference data (typically observed data).

    These parameters can then be applied consistently to multiple datasets (e.g., all Pyro samples).

    Args:
        reference_data: Reference gravity data array (e.g., observed data)
        method: Normalization method
        baseline_forward_model: Baseline forward model (e.g., with mean prior parameters).
                               Required for 'align_to_reference' method to preserve prior variability.
        verbose: Whether to print statistics

    Returns:
        Dictionary containing normalization parameters
    """
    if verbose:
        print(f"Computing {method} normalization parameters from reference data...")

    if method == 'zscore':
        params = {
            'method': method,
            'mean': float(np.mean(reference_data)),
            'std': float(np.std(reference_data))
        }

    elif method == 'robust_zscore':
        median_val = float(np.median(reference_data))
        mad_val = float(np.median(np.abs(reference_data - median_val)))
        params = {
            'method': method,
            'median': median_val,
            'mad': mad_val
        }

    elif method == 'minmax':
        params = {
            'method': method,
            'min': float(np.min(reference_data)),
            'max': float(np.max(reference_data))
        }

    elif method == 'mean_center':
        params = {
            'method': method,
            'mean': float(np.mean(reference_data))
        }

    elif method == 'relative':
        data_range = float(np.max(reference_data) - np.min(reference_data))
        params = {
            'method': method,
            'range': data_range
        }

    elif method == 'align_to_reference':
        # RECOMMENDED for forward models that are shifted/scaled differently from observations
        # CRITICAL: Stores BOTH reference (observed) AND baseline forward model statistics
        # This preserves prior variability while aligning scales

        if baseline_forward_model is None:
            raise ValueError(
                "align_to_reference method requires baseline_forward_model to preserve prior variability. "
                "Compute forward model with mean prior parameters before inference."
            )

        params = {
            'method': method,
            'reference_mean': float(np.mean(reference_data)),
            'reference_std': float(np.std(reference_data)),
            'reference_min': float(np.min(reference_data)),
            'reference_max': float(np.max(reference_data)),
            # CRITICAL: Fixed baseline statistics to preserve prior variability
            'baseline_forward_mean': float(np.mean(baseline_forward_model)),
            'baseline_forward_std': float(np.std(baseline_forward_model))
        }

    else:
        raise ValueError(f"Invalid normalization method: {method}")

    if verbose:
        print(f"  Normalization parameters: {params}")

    return params



def apply_normalization_torch(
        data,  # torch.Tensor or any array-like
        norm_params: Dict[str, float]
):
    """
    Apply pre-computed normalization parameters to data (PyTorch compatible).

    This version works with PyTorch tensors and can be used inside Pyro models.

    Args:
        data: Data tensor to normalize (PyTorch tensor or numpy array)
        norm_params: Normalization parameters from compute_normalization_params()

    Returns:
        Normalized data (same type as input)
    """
    method = norm_params['method']

    if method == 'zscore':
        normalized = (data - norm_params['mean']) / norm_params['std']

    elif method == 'robust_zscore':
        normalized = (data - norm_params['median']) / (1.4826 * norm_params['mad'])

    elif method == 'minmax':
        normalized = (data - norm_params['min']) / (norm_params['max'] - norm_params['min'])

    elif method == 'mean_center':
        normalized = data - norm_params['mean']

    elif method == 'relative':
        normalized = data / norm_params['range']

    elif method == 'align_to_reference':
        # Align data to match reference (observed) distribution
        # This is the RECOMMENDED method when forward and observed are shifted/scaled differently

        # CRITICAL: Use FIXED baseline statistics to preserve prior variability!
        # DO NOT compute mean/std from each sample - that erases the prior's effect!

        baseline_mean = norm_params['baseline_forward_mean']
        baseline_std = norm_params['baseline_forward_std']

        # Step 1: Standardize using FIXED baseline mean/std from initial model
        data_standardized = (data - baseline_mean) / baseline_std

        # Step 2: Scale to match reference std and shift to match reference mean
        # This transformation is LINEAR, so it preserves the variability introduced by priors
        normalized = data_standardized * norm_params['reference_std'] + norm_params['reference_mean']

    else:
        raise ValueError(f"Invalid normalization method: {method}")

    return normalized



def normalize_gravity_pair(
        observed: np.ndarray,
        forward_model: np.ndarray,
        method: Literal['zscore', 'minmax', 'mean_center', 'relative',
        'zscore_independent', 'robust_zscore', 'align_to_reference'] = 'zscore',
        verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, str, Dict[str, float]]:
    """
    Normalize a pair of observed and forward model gravity data using SHARED normalization parameters.
    
    This ensures both datasets are normalized using the same scale (computed from combined data),
    making visual comparison with colorbars more meaningful.
    
    Args:
        observed: Observed gravity data array
        forward_model: Forward model gravity data array
        method: Normalization method to apply
        verbose: Whether to print statistics
    
    Returns:
        Tuple of (observed_normalized, forward_model_normalized, unit_label, norm_params)
        where norm_params contains the shared normalization parameters used
    """
    if verbose:
        print(f"Applying {method} normalization with SHARED parameters...")

    # Compute normalization parameters from COMBINED data
    combined_data = np.concatenate([observed, forward_model])
    
    if method == 'zscore_independent':
        # RECOMMENDED: Normalize each dataset independently, then scale to common range
        # This preserves spatial patterns while making magnitudes comparable

        # Z-score normalize each independently
        obs_mean, obs_std = np.mean(observed), np.std(observed)
        fwd_mean, fwd_std = np.mean(forward_model), np.std(forward_model)

        observed_z = (observed - obs_mean) / obs_std
        forward_z = (forward_model - fwd_mean) / fwd_std

        # Scale both to common range for visualization
        combined_min = min(np.min(observed_z), np.min(forward_z))
        combined_max = max(np.max(observed_z), np.max(forward_z))

        observed_norm = observed_z
        forward_norm = forward_z

        unit_label = 'Z-score (independent)'
        norm_params = {
                'obs_mean': obs_mean, 'obs_std': obs_std,
                'fwd_mean': fwd_mean, 'fwd_std': fwd_std,
                'vmin'    : combined_min, 'vmax': combined_max
        }

    elif method == 'robust_zscore':
        # Robust normalization using median and MAD (Median Absolute Deviation)
        # Less sensitive to outliers than standard z-score

        obs_median = np.median(observed)
        fwd_median = np.median(forward_model)

        obs_mad = np.median(np.abs(observed - obs_median))
        fwd_mad = np.median(np.abs(forward_model - fwd_median))

        # Robust z-score: (x - median) / (1.4826 * MAD)
        # Factor 1.4826 makes MAD comparable to std for normal distribution
        observed_norm = (observed - obs_median) / (1.4826 * obs_mad)
        forward_norm = (forward_model - fwd_median) / (1.4826 * fwd_mad)

        unit_label = 'Robust Z-score'
        norm_params = {
                'obs_median': obs_median, 'obs_mad': obs_mad,
                'fwd_median': fwd_median, 'fwd_mad': fwd_mad
        }

    elif method == 'align_to_reference':
        # RECOMMENDED: Align forward model to match observed data distribution
        # This is perfect when forward and observed are shifted/scaled differently

        # Observed data stays as-is (this is the reference)
        observed_norm = observed

        # Forward model is transformed to match observed distribution
        obs_mean, obs_std = np.mean(observed), np.std(observed)
        fwd_mean, fwd_std = np.mean(forward_model), np.std(forward_model)

        # Standardize forward, then scale/shift to match observed
        forward_standardized = (forward_model - fwd_mean) / fwd_std
        forward_norm = forward_standardized * obs_std + obs_mean

        unit_label = 'Aligned (μGal)'
        norm_params = {
                'method': 'align_to_reference',
                'reference_mean': obs_mean,
                'reference_std': obs_std,
                'forward_original_mean': fwd_mean,
                'forward_original_std': fwd_std
        }

    elif method == 'zscore':
        # Z-score normalization (mean=0, std=1) - use combined mean/std
        data_mean = np.mean(combined_data)
        data_std = np.std(combined_data)

        observed_norm = (observed - data_mean) / data_std
        forward_norm = (forward_model - data_mean) / data_std

        unit_label = 'Z-score'
        norm_params = {'mean': data_mean, 'std': data_std}

    elif method == 'minmax':
        # Min-max normalization (0 to 1) - use combined min/max
        data_min = np.min(combined_data)
        data_max = np.max(combined_data)

        observed_norm = (observed - data_min) / (data_max - data_min)
        forward_norm = (forward_model - data_min) / (data_max - data_min)

        unit_label = 'Normalized [0-1]'
        norm_params = {'min': data_min, 'max': data_max}

    elif method == 'mean_center':
        # Mean centering (subtract mean) - use combined mean
        data_mean = np.mean(combined_data)

        observed_norm = observed - data_mean
        forward_norm = forward_model - data_mean

        unit_label = 'Mean-centered (μGal)'
        norm_params = {'mean': data_mean}

    elif method == 'relative':
        # Relative to range - use combined range
        data_range = np.max(combined_data) - np.min(combined_data)

        observed_norm = observed / data_range
        forward_norm = forward_model / data_range

        unit_label = 'Relative to range'
        norm_params = {'range': data_range}

    else:
        raise ValueError(f"Invalid normalization method: {method}")

    if verbose:
        print(f"  Shared normalization parameters: {norm_params}")
        print(f"  Observed stats (normalized): mean={np.mean(observed_norm):.3f}, std={np.std(observed_norm):.3f}, "
              f"min={np.min(observed_norm):.3f}, max={np.max(observed_norm):.3f}")
        print(f"  Forward stats (normalized):  mean={np.mean(forward_norm):.3f}, std={np.std(forward_norm):.3f}, "
              f"min={np.min(forward_norm):.3f}, max={np.max(forward_norm):.3f}")

    return observed_norm, forward_norm, unit_label, norm_params
