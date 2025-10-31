import numpy as np
from typing import Literal, Tuple, Dict


def normalize_gravity_data(
        data: np.ndarray,
        method: Literal['zscore', 'minmax', 'mean_center', 'relative'] = 'zscore'
) -> Tuple[np.ndarray, str]:
    """
    Normalize gravity data using various methods.
    
    Args:
        data: Input gravity data array
        method: Normalization method to apply:
            - 'zscore': Z-score normalization (mean=0, std=1)
            - 'minmax': Min-max normalization (0 to 1)
            - 'mean_center': Mean centering (subtract mean)
            - 'relative': Relative to range
    
    Returns:
        Tuple of (normalized_data, unit_label)
    
    Raises:
        ValueError: If invalid normalization method is provided
    """
    if method == 'zscore':
        # Z-score normalization (mean=0, std=1)
        data_mean, data_std = np.mean(data), np.std(data)
        normalized = (data - data_mean) / data_std
        unit_label = 'Z-score'

    elif method == 'minmax':
        # Min-max normalization (0 to 1)
        data_min, data_max = np.min(data), np.max(data)
        normalized = (data - data_min) / (data_max - data_min)
        unit_label = 'Normalized [0-1]'

    elif method == 'mean_center':
        # Mean centering (subtract mean)
        normalized = data - np.mean(data)
        unit_label = 'Mean-centered (μGal)'

    elif method == 'relative':
        # Relative to range
        data_range = np.max(data) - np.min(data)
        normalized = data / data_range
        unit_label = 'Relative to range'

    else:
        raise ValueError(f"Invalid normalization method: {method}")

    return normalized, unit_label


def normalize_gravity_pair(
        observed: np.ndarray,
        forward_model: np.ndarray,
        method: Literal['zscore', 'minmax', 'mean_center', 'relative',
        'zscore_independent', 'robust_zscore'] = 'zscore',
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
