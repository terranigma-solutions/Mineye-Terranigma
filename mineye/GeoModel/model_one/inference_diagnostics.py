import arviz as az
from matplotlib import pyplot as plt


def check_mcmc_quality(data: az.InferenceData, min_ess: float = 400, max_rhat: float = 1.01,
                       verbose: bool = True, plot: bool = False) -> bool:
    """
    Automated MCMC quality checks following 2025 best practices.

    Args:
        data: ArviZ InferenceData object
        min_ess: Minimum acceptable ESS (bulk)
        max_rhat: Maximum acceptable R-hat value
        verbose: If True, print detailed diagnostics
        plot: If True, show diagnostic plots

    Returns:
        True if all checks pass, False otherwise
    """
    import numpy as np

    issues = []
    warnings = []

    if verbose:
        print("\n" + "=" * 60)
        print("MCMC DIAGNOSTICS")
        print("=" * 60)

    # Check divergences
    n_divergences = int(data.sample_stats.diverging.sum().values)
    if verbose:
        print(f"Divergences: {n_divergences}")
    if n_divergences > 0:
        issues.append(f"❌ {n_divergences} divergent transitions")
        if verbose:
            print(f"  ⚠️  WARNING: {n_divergences} divergent transitions detected!")

    # Check ESS - properly handle xarray Dataset with multi-dimensional variables
    if verbose:
        print(f"\nESS (bulk):")
    ess_bulk = az.ess(data, method="bulk")
    for var in ess_bulk.data_vars:
        values = np.array(ess_bulk[var].values)  # Ensure it's a numpy array
        # Use minimum ESS if multi-dimensional (worst case)
        ess_val = float(values.min())

        if verbose:
            if values.size == 1:
                print(f"  {var}: {ess_val:.1f}")
            else:
                print(f"  {var}: min={ess_val:.1f}, mean={float(values.mean()):.1f}, max={float(values.max()):.1f}")

        if ess_val < min_ess:
            if ess_val < min_ess / 2:
                issues.append(f"❌ {var}: ESS = {ess_val:.1f} (< {min_ess})")
            else:
                warnings.append(f"⚠️  {var}: ESS = {ess_val:.1f} (< {min_ess})")

    # ESS tail
    if verbose:
        ess_tail = az.ess(data, method="tail")
        print(f"\nESS (tail):")
        for var in ess_tail.data_vars:
            values = ess_tail[var].values
            if values.size == 1:
                print(f"  {var}: {float(values):.1f}")
            else:
                print(f"  {var}: min={float(values.min()):.1f}, mean={float(values.mean()):.1f}, max={float(values.max()):.1f}")

    # Check R-hat - properly handle xarray Dataset with multi-dimensional variables
    if verbose:
        print(f"\nR-hat:")
    rhat = az.rhat(data)
    for var in rhat.data_vars:
        values = np.array(rhat[var].values)  # Ensure it's a numpy array
        # Use maximum R-hat if multi-dimensional (worst case)
        rhat_val = float(values.max())

        if verbose:
            if values.size == 1:
                warning_sym = " ⚠️" if rhat_val > 1.01 else ""
                print(f"  {var}: {rhat_val:.4f}{warning_sym}")
            else:
                warning_sym = " ⚠️" if rhat_val > 1.01 else ""
                print(f"  {var}: max={rhat_val:.4f}{warning_sym}")

        if rhat_val > max_rhat:
            if rhat_val > 1.1:
                issues.append(f"❌ {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")
            else:
                warnings.append(f"⚠️  {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")

    # MCSE (Monte Carlo Standard Error)
    if verbose:
        mcse = az.mcse(data)
        print(f"\nMCSE:")
        for var in mcse.data_vars:
            values = mcse[var].values
            if values.size == 1:
                print(f"  {var}: {float(values):.6f}")
            else:
                print(f"  {var}: mean={float(values.mean()):.6f}")

    # Comprehensive summary table
    if verbose:
        summary = az.summary(
            data,
            var_names=None,  # All variables
            hdi_prob=0.94,  # 94% HDI is standard
            kind="stats"  # Can also use "diagnostics"
        )
        print(f"\nSummary Statistics:\n{summary}")

    # Report results
    if issues:
        print("\n❌ CRITICAL MCMC QUALITY ISSUES:")
        for issue in issues:
            print(f"  {issue}")

    if warnings:
        print("\n⚠️  MCMC QUALITY WARNINGS:")
        for warning in warnings:
            print(f"  {warning}")

    if not issues and not warnings:
        print("\n✅ All MCMC quality checks passed!")
        all_passed = True
    elif not issues:
        print("\n⚠️  Some warnings, but no critical issues")
        all_passed = True
    else:
        print("\n❌ Critical issues detected - consider re-running with adjusted parameters")
        all_passed = False

    # Visual diagnostics
    if plot:
        az.plot_trace(data, compact=True)
        plt.tight_layout()
        plt.show()

        # Rank plots (better than trace plots for convergence)
        try:
            az.plot_rank(data, var_names=["dips"])
            plt.show()
        except Exception as e:
            print(f"Could not create rank plot: {e}")

        # Energy plot (only if energy data is available from HMC/NUTS)
        try:
            if hasattr(data.sample_stats, 'energy'):
                az.plot_energy(data)
                plt.show()
            else:
                print("\nNote: Energy diagnostics not available (requires Pyro with log_prob tracking)")
        except Exception as e:
            print(f"Could not create energy plot: {e}")

    if verbose:
        print("=" * 60 + "\n")

    return all_passed


import pyro.poutine as poutine


def check_likelihood_balance(prob_model, geo_model, y_obs_list):
    """
    Runs the model once to inspect the log_prob magnitudes of specific observation sites.
    """
    print("\n--- Likelihood Balance Check ---")

    # 1. Capture the trace (Run the model and record execution)
    # We wrap the model in a trace and execute it.
    # Note: We don't need 'guide' here because we just want to check the model's 
    # probabilities at the initial state (or prior mean).
    trace = poutine.trace(prob_model).get_trace(
        geo_model=geo_model,
        obs_data=y_obs_list,
        # Assuming you added gravity argument to your model
    )

    # 2. Extract Log Probs for Observation Nodes
    # You need to know the exact "name" you gave the sample site in your model.
    # e.g., dist.Normal(...).to_event(1), name="Gravity Data"

    # -- EnMap Likelihood --
    try:
        enmap_node = trace.nodes["probs_pred"]  # Replace with your exact string name
        enmap_lp = enmap_node["log_prob_sum"]
        print(f"EnMap Log-Likelihood:   {enmap_lp:.2f}")
    except KeyError:
        print("Could not find node 'EnMap Labels'. Check your obs_name string.")

    # -- Gravity Likelihood --
    try:
        gravity_node = trace.nodes[r'$\mu_{gravity}$']  # Replace with your exact string name
        gravity_lp = gravity_node["log_prob_sum"]
        print(f"Gravity Log-Likelihood: {gravity_lp:.2f}")
    except KeyError:
        print("Could not find node 'Gravity Data'.")

    # 3. Compare
    if 'enmap_lp' in locals() and 'gravity_lp' in locals():
        ratio = abs(gravity_lp / enmap_lp)
        print(f"\nRatio (Gravity / EnMap): {ratio:.1f}x")

        if ratio > 100:
            print("WARNING: Gravity is dominating. EnMap will be ignored.")
        elif ratio < 0.01:
            print("WARNING: EnMap is dominating. Gravity will be ignored.")
        else:
            print("STATUS: Balanced. Both datasets should contribute.")