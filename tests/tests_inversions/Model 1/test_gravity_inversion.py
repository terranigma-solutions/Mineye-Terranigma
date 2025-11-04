import os
from functools import partial

import arviz
import arviz as az
import geopandas as gpd
import numpy as np
import pyro
import torch
from matplotlib import pyplot as plt

import gempy_probability as gpp
from gempy_probability.core.samplers_data import NUTSConfig
from mineye.GeoModel.geophysics import align_forward_to_observed
from mineye.GeoModel.model_1_aux import setup_geomodel, baseline, normalize, plot, gempy_viz, create_orientation_modifier
# noinspection PyUnusedImports
from tests import conftest


class TestProbabilisticInversion:
    prior_key = r'dips'

    def test_gravity_inversion(self, simple_geo_model, geophysical_dir, n_samples=50):
        """Test reading and computing a geological model."""
        print("Test gravity inversion...")
        # Use actual gravity measurement device locations
        # * 1) Read gravity data
        gravity_data, observed_gravity_ugal = self._read_gravity(geophysical_dir)

        # * 2) Setup initial Geomodel and normalize forward gravity to the observed gravity
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)
        geo_model.interpolation_options.sigmoid_slope = 100
        baseline_fw_gravity_np = baseline(geo_model)
        norm_params = normalize(baseline_fw_gravity_np, observed_gravity_ugal)

        # * 3) Setup Priors
        model_priors = {
                r'dips': dist.Normal(
                    loc=(torch.ones(geo_model.orientations_copy.xyz.shape[0]) * 10),  # This is just dip 10 degrees
                    scale=torch.tensor(10, dtype=torch.float64),
                    validate_args=True
                )
        }

        # TODO: Here we could add density and range

        # * 4) Set up Deterministics
        pre_forward_dets = {
                "dips_degrees": lambda samples, gm: samples["dips"],  # Just pass through
        }

        post_forward_dets = {
                "gravity_response_raw": lambda samples, gm, sol: sol.gravity,  # Store raw gravity
                "gravity_response"    : lambda samples, gm, sol: align_forward_to_observed(-sol.gravity, norm_params),  # Normalized!
                "mean_gravity"        : lambda samples, gm, sol: torch.mean(align_forward_to_observed(-sol.gravity, norm_params)),
                "max_gravity"         : lambda samples, gm, sol: torch.max(align_forward_to_observed(-sol.gravity, norm_params), 0),
        }

        # * 5) Set up likelihood functions
        likelihood_fn = generate_multigravity_likelihood_diagonal(
            norm_params=norm_params
        )

        # * 6) Set up Pyro model
        prob_model: gpp.GemPyPyroModel = gpp.make_gempy_pyro_model_extended(
            priors=model_priors,
            set_interp_input_fn=create_orientation_modifier(key=TestProbabilisticInversion.prior_key),
            likelihood_fn=likelihood_fn,
            pre_forward_deterministics=pre_forward_dets,
            post_forward_deterministics=post_forward_dets,
            obs_name="Gravity Measurement"
        )

        # * 7) Run predictive
        gravity_observations_tensor = torch.tensor(observed_gravity_ugal)
        trace = trace_pyro_model(prob_model, geo_model, torch.tensor(observed_gravity_ugal, dtype=torch.float64))
        compute_prior_predictive = True
        if compute_prior_predictive:
            prior_inference_data: az.InferenceData = gpp.run_predictive(
                prob_model=prob_model,
                geo_model=geo_model,
                y_obs_list=gravity_observations_tensor,
                n_samples=10,
                plot_trace=True
            )

        # * 8) Run inference

        data = gpp.run_nuts_inference(
            prob_model=prob_model,
            geo_model=geo_model,
            y_obs_list=gravity_observations_tensor,
            config=NUTSConfig(
                step_size=0.0001,
                adapt_step_size=True,
                target_accept_prob=0.65,
                max_tree_depth=5,
                init_strategy='median',
                num_samples=20,
                warmup_steps=5,
            ),
            plot_trace=True,
            run_posterior_predictive=True
        )

        if compute_prior_predictive:
            data.extend(prior_inference_data)

        data.to_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data.nc"))

    def _read_gravity(self, geophysical_dir):
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

        # Take a spatially distributed subset of measurements
        n_points = 20  # Adjust this number to control how many points you want
        xy_coords = gravity_data.geometry.apply(lambda p: (p.x, p.y)).to_list()
        xy_array = np.array(xy_coords)

        # Use K-means clustering to get well-distributed points
        from sklearn.cluster import KMeans
        kmeans = KMeans(n_clusters=n_points, random_state=42)
        kmeans.fit(xy_array)

        # Find the closest points to cluster centers
        from scipy.spatial.distance import cdist
        centers = kmeans.cluster_centers_
        distances = cdist(centers, xy_array)
        indices = [np.argmin(dist) for dist in distances]

        # Filter the gravity data
        observed_gravity = observed_gravity[indices]
        gravity_data = gravity_data.iloc[indices]

        observed_gravity_ugal = observed_gravity * 1000
        return gravity_data, observed_gravity_ugal

    def test_run_diagnostics(self):
        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data.nc"))
        check_mcmc_quality(data)

        # After MCMC - Modern ArviZ diagnostics (2025)
        print("\n" + "=" * 60)
        print("MCMC DIAGNOSTICS")
        print("=" * 60)

        # 1. Divergences
        n_divergences = int(data.sample_stats.diverging.sum().values)
        print(f"Divergences: {n_divergences}")
        if n_divergences > 0:
            print(f"  ⚠️  WARNING: {n_divergences} divergent transitions detected!")

        # 2. ESS (Effective Sample Size) - use bulk and tail
        ess_bulk = az.ess(data, method="bulk")
        ess_tail = az.ess(data, method="tail")
        print(f"\nESS (bulk):")
        for var in ess_bulk.data_vars:
            values = ess_bulk[var].values
            if values.size == 1:
                print(f"  {var}: {float(values):.1f}")
            else:
                # For multi-dimensional variables, show summary
                print(f"  {var}: min={float(values.min()):.1f}, mean={float(values.mean()):.1f}, max={float(values.max()):.1f}")

        print(f"\nESS (tail):")
        for var in ess_tail.data_vars:
            values = ess_tail[var].values
            if values.size == 1:
                print(f"  {var}: {float(values):.1f}")
            else:
                print(f"  {var}: min={float(values.min()):.1f}, mean={float(values.mean()):.1f}, max={float(values.max()):.1f}")

        # 3. R-hat (should be < 1.01, ideally < 1.05)
        rhat = az.rhat(data)
        print(f"\nR-hat:")
        for var in rhat.data_vars:
            values = rhat[var].values
            if values.size == 1:
                rhat_val = float(values)
                warning = " ⚠️" if rhat_val > 1.01 else ""
                print(f"  {var}: {rhat_val:.4f}{warning}")
            else:
                max_rhat = float(values.max())
                warning = " ⚠️" if max_rhat > 1.01 else ""
                print(f"  {var}: max={max_rhat:.4f}{warning}")

        # 4. MCSE (Monte Carlo Standard Error)
        mcse = az.mcse(data)
        print(f"\nMCSE:")
        for var in mcse.data_vars:
            values = mcse[var].values
            if values.size == 1:
                print(f"  {var}: {float(values):.6f}")
            else:
                print(f"  {var}: mean={float(values.mean()):.6f}")

        # 5. Comprehensive summary table
        summary = az.summary(
            data,
            var_names=None,  # All variables
            hdi_prob=0.94,  # 94% HDI is standard
            kind="stats"  # Can also use "diagnostics"
        )
        print(f"\nSummary Statistics:\n{summary}")

        # 6. Visual diagnostics (recommended)
        if True:  # Set to False to skip plots
            az.plot_trace(data, compact=True)
            plt.tight_layout()
            plt.show()

            # 7. Rank plots (better than trace plots for convergence)
            try:
                az.plot_rank(data, var_names=["dips"])
                plt.show()
            except Exception as e:
                print(f"Could not create rank plot: {e}")

            # 8. Energy plot (only if energy data is available from HMC/NUTS)
            try:
                if hasattr(data.sample_stats, 'energy'):
                    az.plot_energy(data)
                    plt.show()
                else:
                    print("\nNote: Energy diagnostics not available (requires Pyro with log_prob tracking)")
            except Exception as e:
                print(f"Could not create energy plot: {e}")

        print("=" * 60 + "\n")

    def test_run_analysis(self, simple_geo_model, geophysical_dir):

        data = az.from_netcdf(os.path.join(os.path.dirname(__file__), "arviz_data.nc"))

        gravity_data, observed_gravity_ugal = self._read_gravity(geophysical_dir)
        geo_model, xy_ravel = setup_geomodel(gravity_data, simple_geo_model)

        # # Posterior predictive checks
        az.plot_ppc(data, num_pp_samples=20)
        # * 9) Analysis inference
        gravity_samples_norm, unit_label = plot(
            gravity_samples_norm=data.prior[r'gravity_response'].values[0, :],  # (n_samples, n_devices)
            observed_gravity_ugal=observed_gravity_ugal,
            xy_ravel=xy_ravel
        )

        # * 9) Analysis Gempy Model
        gempy_viz(geo_model, data)

    def test_gravity_duplicates(self, geophysical_dir):
        """Test reading and computing a geological model."""

        # Use actual gravity measurement device locations
        # * 1) Read gravity data
        gravity_data = gpd.read_file(os.path.join(geophysical_dir, 'cleaned_gravity_data.geojson'))
        observed_gravity = gravity_data['VALU_BOU267'].values  # in mGal

        # * 1b) CHECK FOR DUPLICATES
        xy_coords = gravity_data[['geometry']].apply(lambda row: (row.geometry.x, row.geometry.y), axis=1)
        xy_array = np.array(xy_coords.tolist())

        # Find unique locations (within 1m tolerance)
        from scipy.spatial.distance import pdist, squareform
        distances = squareform(pdist(xy_array))
        np.fill_diagonal(distances, np.inf)  # Ignore self-distances

        duplicates = np.any(distances < 1.0, axis=1)  # Locations within 1m
        if np.any(duplicates):
            print(f"WARNING: Found {duplicates.sum()} stations with duplicates/near-duplicates")
            print("Removing duplicates...")

            # Keep first occurrence of each duplicate set
            keep_mask = np.ones(len(xy_array), dtype=bool)
            for i in np.where(duplicates)[0]:
                if keep_mask[i]:
                    # Find all duplicates of this point
                    dups = distances[i] < 1.0
                    # Keep first, remove rest
                    dup_indices = np.where(dups)[0]
                    if len(dup_indices) > 0:
                        keep_mask[dup_indices[1:]] = False

            # Filter data
            gravity_data = gravity_data[keep_mask]
            observed_gravity = gravity_data['VALU_BOU267'].values
            observed_gravity_ugal = observed_gravity * 1000
            print(f"Kept {keep_mask.sum()} unique stations out of {len(keep_mask)}")


import gempy as gp
import pyro.distributions as dist


def check_mcmc_quality(data: az.InferenceData, min_ess: float = 400, max_rhat: float = 1.01) -> bool:
    """
    Automated MCMC quality checks following 2025 best practices.

    Args:
        data: ArviZ InferenceData object
        min_ess: Minimum acceptable ESS (bulk)
        max_rhat: Maximum acceptable R-hat value

    Returns:
        True if all checks pass, False otherwise
    """
    import matplotlib.pyplot as plt
    import numpy as np

    issues = []
    warnings = []

    # Check divergences
    n_divergences = int(data.sample_stats.diverging.sum().values)
    if n_divergences > 0:
        issues.append(f"❌ {n_divergences} divergent transitions")

    # Check ESS - properly handle xarray Dataset with multi-dimensional variables
    ess_bulk = az.ess(data, method="bulk")
    for var in ess_bulk.data_vars:
        values = np.array(ess_bulk[var].values)  # Ensure it's a numpy array
        # Use minimum ESS if multi-dimensional (worst case)
        ess_val = float(values.min())

        if ess_val < min_ess:
            if ess_val < min_ess / 2:
                issues.append(f"❌ {var}: ESS = {ess_val:.1f} (< {min_ess})")
            else:
                warnings.append(f"⚠️  {var}: ESS = {ess_val:.1f} (< {min_ess})")

    # Check R-hat - properly handle xarray Dataset with multi-dimensional variables
    rhat = az.rhat(data)
    for var in rhat.data_vars:
        values = np.array(rhat[var].values)  # Ensure it's a numpy array
        # Use maximum R-hat if multi-dimensional (worst case)
        rhat_val = float(values.max())

        if rhat_val > max_rhat:
            if rhat_val > 1.1:
                issues.append(f"❌ {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")
            else:
                warnings.append(f"⚠️  {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")

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
        return True
    elif not issues:
        print("\n⚠️  Some warnings, but no critical issues")
        return True
    else:
        print("\n❌ Critical issues detected - consider re-running with adjusted parameters")
        return False


def check_mcmc_quality_(data: az.InferenceData, min_ess: float = 400, max_rhat: float = 1.01) -> bool:
    """
    Automated MCMC quality checks following 2025 best practices.

    Args:
        data: ArviZ InferenceData object
        min_ess: Minimum acceptable ESS (bulk)
        max_rhat: Maximum acceptable R-hat value

    Returns:
        True if all checks pass, False otherwise
    """
    import matplotlib.pyplot as plt

    issues = []
    warnings = []

    # Check divergences
    n_divergences = int(data.sample_stats.diverging.sum().values)
    if n_divergences > 0:
        issues.append(f"❌ {n_divergences} divergent transitions")

    # Check ESS - properly handle xarray Dataset
    ess_bulk = az.ess(data, method="bulk")
    for var in ess_bulk.data_vars:
        ess_val = float(ess_bulk[var].values)
        if ess_val < min_ess:
            if ess_val < min_ess / 2:
                issues.append(f"❌ {var}: ESS = {ess_val:.1f} (< {min_ess})")
            else:
                warnings.append(f"⚠️  {var}: ESS = {ess_val:.1f} (< {min_ess})")

    # Check R-hat - properly handle xarray Dataset
    rhat = az.rhat(data)
    for var in rhat.data_vars:
        rhat_val = float(rhat[var].values)
        if rhat_val > max_rhat:
            if rhat_val > 1.1:
                issues.append(f"❌ {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")
            else:
                warnings.append(f"⚠️  {var}: R-hat = {rhat_val:.4f} (> {max_rhat})")

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
        return True
    elif not issues:
        print("\n⚠️  Some warnings, but no critical issues")
        return True
    else:
        print("\n❌ Critical issues detected - consider re-running with adjusted parameters")
        return False


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


# Keep your gaussian_kernel function as-is
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


def debug_gradient_flow(likelihood_fn, solutions, obs_data):
    """
    Systematically check where gradient flow breaks.
    
    Call this BEFORE running NUTS to diagnose issues.
    """
    import torch

    print("\n" + "=" * 80)
    print("GRADIENT FLOW DIAGNOSTIC")
    print("=" * 80)

    # Get the distribution
    lik_dist = likelihood_fn(solutions)

    # Compute log probability
    try:
        log_prob = lik_dist.log_prob(obs_data)
        print(f"✓ Log probability computed: {log_prob.item():.4f}")
    except Exception as e:
        print(f"✗ Failed to compute log_prob: {e}")
        return

    # Check if we can compute gradients
    print(f"\nLog prob requires_grad: {log_prob.requires_grad}")
    print(f"Log prob grad_fn: {log_prob.grad_fn}")

    if not log_prob.requires_grad:
        print("\n⚠️  WARNING: log_prob does not require gradients!")
        print("This means NUTS cannot compute gradients for sampling.")

        # Trace back through the computation graph
        print("\nTracing back through computation:")
        print(f"  Distribution type: {type(lik_dist).__name__}")

        # Check distribution parameters
        if hasattr(lik_dist, 'loc'):
            print(f"  loc requires_grad: {lik_dist.loc.requires_grad}")
            print(f"  loc grad_fn: {lik_dist.loc.grad_fn}")

        if hasattr(lik_dist, 'scale_tril'):
            print(f"  scale_tril requires_grad: {lik_dist.scale_tril.requires_grad}")
            print(f"  scale_tril grad_fn: {lik_dist.scale_tril.grad_fn}")

        if hasattr(lik_dist, 'covariance_matrix'):
            print(f"  covariance_matrix requires_grad: {lik_dist.covariance_matrix.requires_grad}")
            print(f"  covariance_matrix grad_fn: {lik_dist.covariance_matrix.grad_fn}")

        if hasattr(lik_dist, 'df'):
            print(f"  df requires_grad: {lik_dist.df.requires_grad}")
            print(f"  df grad_fn: {lik_dist.df.grad_fn}")
    else:
        print("\n✓ Gradient flow is intact!")

        # Try computing gradients
        try:
            # Create dummy parameters to test backward pass
            test_grad = torch.autograd.grad(
                outputs=log_prob,
                inputs=[p for p in lik_dist.parameters() if p.requires_grad],
                allow_unused=True
            )
            print(f"✓ Backward pass successful, got {len(test_grad)} gradients")
        except Exception as e:
            print(f"✗ Backward pass failed: {e}")

    print("=" * 80 + "\n")


# Usage in your test:
def test_gravity_inversion(self, simple_geo_model, geophysical_dir, n_samples=50):
    # ... setup code ...

    # Before running NUTS, test gradient flow:
    from gempy_probability.modules.forwards import run_gempy_forward
    from gempy.modules.data_manipulation import interpolation_input_from_structural_frame

    # Create a test forward pass
    test_interp_input = interpolation_input_from_structural_frame(geo_model)
    test_solutions = run_gempy_forward(test_interp_input, geo_model)

    # Debug the likelihood
    debug_gradient_flow(likelihood_fn, test_solutions, torch.tensor(observed_gravity_ugal, dtype=torch.float64))

    # ... continue with NUTS ...


def trace_pyro_model(prob_model, geo_model, obs_data, print_full=True):
    """
    Use Pyro's trace functionality to see all sample sites and deterministics.

    This shows you EXACTLY what Pyro sees.
    """
    from pyro import poutine
    import torch

    print("\n" + "=" * 80)
    print("PYRO MODEL TRACE")
    print("=" * 80 + "\n")

    # Trace the model execution
    trace = poutine.trace(prob_model).get_trace(geo_model, obs_data)

    print(f"{'Site Name':<30} | {'Type':<15} | {'Shape':<15} | {'Grad?':<10} | {'Grad Fn'}")
    print("-" * 100)

    for name, node in trace.nodes.items():
        if node['type'] == 'sample':
            value = node['value']

            if isinstance(value, torch.Tensor):
                shape_str = str(tuple(value.shape))
                grad_str = "✓" if value.requires_grad else "✗"
                grad_fn_str = str(value.grad_fn)[:30] if value.grad_fn else "None"

                print(f"{name:<30} | {node['type']:<15} | {shape_str:<15} | {grad_str:<10} | {grad_fn_str}")

                # Flag problematic ones
                if not value.requires_grad and name not in ['obs', 'Gravity Measurement']:
                    print(f"  ⚠️  WARNING: Sample site '{name}' has no gradient flow!")
            else:
                print(f"{name:<30} | {node['type']:<15} | {str(type(value)):<15} | {'N/A':<10} | N/A")

    if print_full:
        print("\n" + "=" * 80)
        print("FULL TRACE DETAILS")
        print("=" * 80 + "\n")
        print(trace.format_shapes())

    print("\n" + "=" * 80 + "\n")

    return trace


def generate_multigravity_likelihood_diagonal(norm_params):
    """
    Use independent Normal distributions instead of multivariate.
    Much more stable and faster.
    """

    def likelihood_fn(solutions: gp.data.Solutions) -> dist.Distribution:
        simulated_geophysics = align_forward_to_observed(-solutions.gravity, norm_params)
        pyro.deterministic(r'$\mu_{gravity}$', simulated_geophysics)

        # Sample noise standard deviation
        sigma = pyro.sample(
            "sigma",
            dist.HalfNormal(torch.tensor(500.0, dtype=torch.float64))  # 100 µGal noise
        )

        # Independent Normal likelihood (much more stable!)
        return dist.Normal(simulated_geophysics, sigma).to_event(1)

    return likelihood_fn
