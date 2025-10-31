import numpy as np
from typing import Literal, Tuple


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
        method: Literal['zscore', 'minmax', 'mean_center', 'relative'] = 'zscore',
        verbose: bool = True
) -> Tuple[np.ndarray, np.ndarray, str]:
    """
    Normalize a pair of observed and forward model gravity data.
    
    Args:
        observed: Observed gravity data array
        forward_model: Forward model gravity data array
        method: Normalization method to apply
        verbose: Whether to print statistics
    
    Returns:
        Tuple of (observed_normalized, forward_model_normalized, unit_label)
    """
    if verbose:
        print(f"Applying {method} normalization...")

    observed_norm, unit_label = normalize_gravity_data(observed, method)
    forward_norm, _ = normalize_gravity_data(forward_model, method)

    if verbose:
        print(f"  Observed stats (normalized): mean={np.mean(observed_norm):.3f}, std={np.std(observed_norm):.3f}")
        print(f"  Forward stats (normalized):  mean={np.mean(forward_norm):.3f}, std={np.std(forward_norm):.3f}")

    return observed_norm, forward_norm, unit_label