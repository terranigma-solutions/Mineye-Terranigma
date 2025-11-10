from typing import Dict, Optional, Literal

import numpy as np
import torch


def compute_alignment_params(
        observed: np.ndarray,
        baseline_forward: Optional[np.ndarray] = None,
        verbose: bool = True,
        method: Literal["align_to_reference", "quantile_align"] = "quantile_align",
        n_quantiles: int = 11,
        extrapolation_buffer: float = 0.2
) -> Dict[str, np.ndarray | float | str]:
    """
    Compute parameters to align forward gravity fields to observed distribution.
    method:
      - "align_to_reference": linear mean/std alignment (assumes near-normal)
      - "quantile_align": robust monotonic alignment via piecewise linear quantile mapping

    extrapolation_buffer: fraction to extend the quantile range beyond baseline (e.g., 0.2 = 20%)
    """
    if verbose:
        print(f"Computing {method} alignment parameters...")

    if method == "align_to_reference":
        params = {
            "method": method,
            "reference_mean": float(np.mean(observed)),
            "reference_std": float(np.std(observed)),
        }
        if baseline_forward is not None:
            params["baseline_forward_mean"] = float(np.mean(baseline_forward))
            params["baseline_forward_std"] = float(np.std(baseline_forward))
        if verbose:
            print(f"  Alignment params: {params}")
        return params

    # Robust quantile alignment
    if baseline_forward is None:
        raise ValueError("quantile_align requires baseline_forward.")

    q = np.linspace(0, 1, n_quantiles)
    obs_qv = np.quantile(observed, q)
    base_qv = np.quantile(baseline_forward, q)

    # Ensure strict monotonicity to avoid flat segments
    eps = 1e-9
    for arr in (obs_qv, base_qv):
        # Enforce increasing
        for i in range(1, len(arr)):
            if arr[i] <= arr[i-1]:
                arr[i] = arr[i-1] + eps

    # Add extrapolation buffer to allow exploration beyond baseline range
    if extrapolation_buffer > 0:
        base_range = base_qv[-1] - base_qv[0]
        obs_range = obs_qv[-1] - obs_qv[0]

        # Extend baseline quantiles
        base_qv_extended = np.concatenate([
            [base_qv[0] - extrapolation_buffer * base_range],
            base_qv,
            [base_qv[-1] + extrapolation_buffer * base_range]
        ])

        # Extend observed quantiles (linear extrapolation)
        obs_qv_extended = np.concatenate([
            [obs_qv[0] - extrapolation_buffer * obs_range],
            obs_qv,
            [obs_qv[-1] + extrapolation_buffer * obs_range]
        ])

        base_qv = base_qv_extended
        obs_qv = obs_qv_extended

    params = {
        "method": "quantile_align",
        "quantiles": q,
        "baseline_quantile_values": base_qv.astype(float),
        "observed_quantile_values": obs_qv.astype(float),
    }
    if verbose:
        print("  Quantile alignment params ready "
              f"(n_knots={n_quantiles}, buffer={extrapolation_buffer:.1%}): ranges "
              f"[{base_qv[0]:.2f},{base_qv[-1]:.2f}] -> [{obs_qv[0]:.2f},{obs_qv[-1]:.2f}]")
    return params


def align_forward_to_observed(
        forward: torch.Tensor,
        params: Dict[str, np.ndarray | float | str],
        eps: float = 1e-6,
        negate_forward: bool = False
) -> torch.Tensor:
    """
    Align forward gravity to observed using params from compute_alignment_params.
    Supports:
      - method="align_to_reference": linear mean/std
      - method="quantile_align": robust piecewise-linear quantile mapping
    """
    x = -forward if negate_forward else forward
    method = params.get("method", "align_to_reference")

    if method == "align_to_reference":
        ref_mu = float(params["reference_mean"])
        ref_sigma = float(params["reference_std"])
        base_mu = float(params.get("baseline_forward_mean"))
        base_sigma = float(params.get("baseline_forward_std"))
        base_sigma = max(base_sigma, eps)
        standardized = (x - base_mu) / base_sigma
        return standardized * ref_sigma + ref_mu

    # Quantile alignment (robust)
    base_qv = torch.tensor(params["baseline_quantile_values"], dtype=x.dtype, device=x.device)
    obs_qv = torch.tensor(params["observed_quantile_values"], dtype=x.dtype, device=x.device)

    # Piecewise-linear interpolation/extrapolation: map baseline value -> observed value
    # Extended quantiles allow exploration beyond baseline range
    idx = torch.searchsorted(base_qv, x, right=True).clamp(min=1, max=base_qv.numel()-1)
    x0 = base_qv[idx - 1]
    x1 = base_qv[idx]
    y0 = obs_qv[idx - 1]
    y1 = obs_qv[idx]
    w = (x - x0) / torch.clamp(x1 - x0, min=eps)
    y = y0 + w * (y1 - y0)
    return y