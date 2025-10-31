import numpy as np
from typing import Literal, Tuple, Dict, Optional

import numpy as np
from typing import Dict, Optional

import torch


def compute_alignment_params(
        observed: np.ndarray,
        baseline_forward: Optional[np.ndarray] = None,
        verbose: bool = True
) -> Dict[str, float]:
    """
    Compute parameters to linearly align any forward gravity field to the observed distribution.
    If baseline_forward is provided, its fixed stats are stored to preserve prior variability.
    """
    if verbose:
        print("Computing alignment parameters from observed data" +
              (" and baseline forward" if baseline_forward is not None else "") + "...")

    params = {
            "method"        : "align_to_reference",
            "reference_mean": float(np.mean(observed)),
            "reference_std" : float(np.std(observed)),
    }

    if baseline_forward is not None:
        params["baseline_forward_mean"] = float(np.mean(baseline_forward))
        params["baseline_forward_std"] = float(np.std(baseline_forward))

    if verbose:
        print(f"  Alignment params: {params}")

    return params


# ... existing code ...

def align_forward_to_observed(
        forward: torch.Tensor,
        params: Dict[str, float]
) -> np.ndarray:
    """
    Linearly map a forward gravity field to the observed distribution using precomputed params.
    If baseline stats exist, they are used to preserve prior variability; otherwise,
    the forward’s own stats are used (still a linear alignment).
    """
    ref_mu = params["reference_mean"]
    ref_sigma = params["reference_std"]

    # Prefer fixed baseline stats if available
    base_mu = params.get("baseline_forward_mean", float(torch.mean(forward)))
    base_sigma = params.get("baseline_forward_std", float(torch.std(forward)))

    # Standardize with baseline (or forward) stats, then scale/shift to observed
    standardized = (forward - base_mu) / base_sigma
    aligned = standardized * ref_sigma + ref_mu
    return aligned


# ... existing code ...
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
                'mean'  : float(np.mean(reference_data)),
                'std'   : float(np.std(reference_data))
        }

    elif method == 'robust_zscore':
        median_val = float(np.median(reference_data))
        mad_val = float(np.median(np.abs(reference_data - median_val)))
        params = {
                'method': method,
                'median': median_val,
                'mad'   : mad_val
        }

    elif method == 'minmax':
        params = {
                'method': method,
                'min'   : float(np.min(reference_data)),
                'max'   : float(np.max(reference_data))
        }

    elif method == 'mean_center':
        params = {
                'method': method,
                'mean'  : float(np.mean(reference_data))
        }

    elif method == 'relative':
        data_range = float(np.max(reference_data) - np.min(reference_data))
        params = {
                'method': method,
                'range' : data_range
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
                'method'               : method,
                'reference_mean'       : float(np.mean(reference_data)),
                'reference_std'        : float(np.std(reference_data)),
                'reference_min'        : float(np.min(reference_data)),
                'reference_max'        : float(np.max(reference_data)),
                # CRITICAL: Fixed baseline statistics to preserve prior variability
                'baseline_forward_mean': float(np.mean(baseline_forward_model)),
                'baseline_forward_std' : float(np.std(baseline_forward_model))
        }

    else:
        raise ValueError(f"Invalid normalization method: {method}")

    if verbose:
        print(f"  Normalization parameters: {params}")

    return params


