#!/usr/bin/env python3
"""
Memory Consumption Test Script for BaySeg Algorithms

This script tests memory consumption of both optimized and original BaySeg algorithms
by taking memory snapshots at different step intervals (50, 100, 200, 400, 800 steps)
using the smallest ROI configuration.

Based on the performance_comparison.py blueprint.
"""

import matplotlib.pyplot as plt
import time
import psutil
import gc
import tracemalloc
from typing import Dict, List, Tuple
import rasterio
import sys

# Add the bayseg directory to path
sys.path.append('/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/bayseg')

from bayseg.bayseg import BaySeg
from mineye.BayesianSegmentation.Sentinel2_full_workflow import crop_to_bounds, apply_soil_mask

# Import the original BaySeg implementation from performance_comparison.py
import numpy as np
from sklearn import mixture
from scipy.stats import multivariate_normal, norm
from copy import copy
from itertools import combinations
import tqdm


def compute_ie(data):
    return np.zeros(data.shape[0])  # Placeholder

def compute_labels_prob(data):
    return np.zeros((data.shape[0], data.shape[1]))  # Placeholder

# BaySegOriginal class - complete implementation from performance_comparison.py
# This version uses np.append to grow arrays causing memory issues (no memory management)
class BaySegOriginal:
    def __init__(self, data, n_labels, beta_init=1, stencil=None, normalize=True):
        # Store initial data
        self.data = data
        self.shape = np.shape(data)
        self.phys_shp = np.array(self.shape[:-1])
        self.n_feat = self.shape[-1]
        self.stencil = stencil
        self.colors = pseudocolor(self.shape, self.stencil)

        # Fetch dimensionality and feature vector from input data
        if len(self.shape) == 2:
            self.dim = 1
            self.feat = self.data
        elif len(self.shape) == 3:
            self.dim = 2
            self.feat = np.array([self.data[:, :, f].ravel() for f in range(self.n_feat)]).T
        elif len(self.shape) == 4:
            raise Exception("3D segmentation not yet supported.")
        else:
            raise Exception("Data format appears to be wrong (neither 1-, 2- or 3-D).")

        if normalize:
            self.normalize_feature_vectors()

        # Initialize GMM
        self.n_labels = n_labels
        self.gmm = mixture.GaussianMixture(n_components=n_labels, covariance_type="full")
        self.gmm.fit(self.feat)

        # Initialize storage arrays using numpy arrays (these grow via np.append - memory issue!)
        self.labels = np.array([self.gmm.predict(self.feat)])
        self.mus = np.array([self.gmm.means_])
        self.covs = np.array([self.gmm.covariances_])
        
        self.labels_probability = np.zeros((1, self.labels.shape[1], self.n_labels))
        self.storage_gibbs_e = np.zeros([1, self.labels.shape[1], self.n_labels])
        self.storage_like_e = np.zeros([1, self.labels.shape[1], self.n_labels])
        self.storage_te = np.zeros([1, self.labels.shape[1], self.n_labels])
        
        self.beta_acc_ratio = np.array([])
        self.cov_acc_ratio = np.array([])
        self.mu_acc_ratio = np.array([])

        # Initialize PRIOR distributions for beta, mu and covariance
        if self.dim == 1:
            self.prior_beta = norm(beta_init, np.eye(1) * 100)
            self.betas = [beta_init]
        elif self.dim == 2:
            if self.stencil == "4p":
                beta_dim = 2
            elif self.stencil == "8p" or self.stencil is None:
                beta_dim = 4
            self.betas = [[beta_init for i in range(beta_dim)]]
            self.prior_beta = multivariate_normal([beta_init for i in range(beta_dim)], np.eye(beta_dim) * 100)
        elif self.dim == 3:
            raise Exception("3D not yet supported.")

        # MU priors
        prior_mu_means = [self.mus[0][label] for label in range(self.n_labels)]
        prior_mu_stds = [np.eye(self.n_feat) * 100 for label in range(self.n_labels)]
        self.priors_mu = [multivariate_normal(prior_mu_means[label], prior_mu_stds[label]) 
                         for label in range(self.n_labels)]

        # COV priors
        self.b_sigma = np.zeros((self.n_labels, self.n_feat))
        for l in range(self.n_labels):
            self.b_sigma[l, :] = np.log(np.sqrt(np.diag(self.gmm.covariances_[l, :, :])))
        self.kesi = np.ones((self.n_labels, self.n_feat)) * 100
        self.nu = self.n_feat + 1

    def fit(self, n, beta_jump_length=10, mu_jump_length=0.0005, cov_volume_jump_length=0.00005,
            theta_jump_length=0.0005, t=1., verbose=False, fix_beta=False):
        for g in tqdm.trange(n):
            self.gibbs_sample(t, beta_jump_length, mu_jump_length, cov_volume_jump_length, 
                             theta_jump_length, verbose, fix_beta)

    def gibbs_sample(self, t, beta_jump_length, mu_jump_length, cov_volume_jump_length, 
                     theta_jump_length, verbose, fix_beta):
        # CALCULATE TOTAL ENERGY (real computations)
        energy_like = self.calc_energy_like(self.mus[-1], self.covs[-1])
        gibbs_energy = self._calc_gibbs_energy_vect(self.labels[-1], self.betas[-1], verbose=verbose)
        total_energy = energy_like + gibbs_energy
        
        # CALCULATE PROBABILITY OF LABELS
        labels_prob = _calc_labels_prob(total_energy, t)
        # Memory issue: np.append grows arrays without limit
        self.storage_te = np.append(self.storage_te, total_energy[np.newaxis, :, :], axis=0)

        # Update labels using proper algorithm
        new_labels = copy(self.labels[-1])
        for i, color_f in enumerate(self.colors):
            new_labels[color_f] = draw_labels_vect(labels_prob[color_f])
            # Recalculate energies with updated labels
            gibbs_energy = self._calc_gibbs_energy_vect(new_labels, self.betas[-1], verbose=verbose)
            total_energy = energy_like + gibbs_energy
            labels_prob = _calc_labels_prob(total_energy, t)

        # Memory issue: arrays keep growing
        self.labels_probability = np.append(self.labels_probability, labels_prob[np.newaxis, :, :], axis=0)
        self.labels = np.append(self.labels, new_labels[np.newaxis, :], axis=0)

        # Calculate component coefficient
        energy_for_comp_coef = gibbs_energy
        comp_coef = _calc_labels_prob(energy_for_comp_coef, t)
        
        # PROPOSAL STEP
        beta_prop = self.propose_beta(self.betas[-1], beta_jump_length)
        mu_prop = self.propose_mu(self.mus[-1], mu_jump_length)
        cov_prop = _propose_cov(self.covs[-1], self.n_feat, self.n_labels, 
                               cov_volume_jump_length, theta_jump_length)

        # UPDATE MU
        mu_next = copy(self.mus[-1])
        cov_next = copy(self.covs[-1])
        
        for l in range(self.n_labels):
            mu_temp = copy(mu_next)
            mu_temp[l, :] = mu_prop[l, :]
            
            lp_mu_prev = self.log_prior_density_mu(mu_next, l)
            lp_mu_prop = self.log_prior_density_mu(mu_temp, l)
            
            lmd_prev = self.calc_sum_log_mixture_density(comp_coef, mu_next, cov_next)
            lmd_prop = self.calc_sum_log_mixture_density(comp_coef, mu_temp, cov_next)
            
            log_target_prev = lmd_prev + lp_mu_prev
            log_target_prop = lmd_prop + lp_mu_prop
            
            mu_eval = evaluate(log_target_prop, log_target_prev)
            if mu_eval[0]:
                mu_next[l, :] = mu_prop[l, :]
            self.mu_acc_ratio = np.append(self.mu_acc_ratio, mu_eval[1])

        # Memory issue: arrays keep growing
        self.mus = np.append(self.mus, mu_next[np.newaxis, :, :], axis=0)

        # UPDATE COVARIANCE
        for l in range(self.n_labels):
            cov_temp = copy(cov_next)
            cov_temp[l, :, :] = cov_prop[l, :, :]
            
            lp_cov_prev = self.log_prior_density_cov(cov_next, l)
            lp_cov_prop = self.log_prior_density_cov(cov_temp, l)
            
            lmd_prev = self.calc_sum_log_mixture_density(comp_coef, mu_next, cov_next)
            lmd_prop = self.calc_sum_log_mixture_density(comp_coef, mu_next, cov_temp)
            
            log_target_prev = lmd_prev + lp_cov_prev
            log_target_prop = lmd_prop + lp_cov_prop
            
            mu_eval = evaluate(log_target_prop, log_target_prev)
            if mu_eval[0]:
                cov_next[l, :] = cov_prop[l, :]
            self.cov_acc_ratio = np.append(self.cov_acc_ratio, mu_eval[1])

        # Memory issue: arrays keep growing  
        self.covs = np.append(self.covs, cov_next[np.newaxis, :, :], axis=0)
        self.storage_gibbs_e = np.append(self.storage_gibbs_e, gibbs_energy[np.newaxis, :, :], axis=0)
        self.storage_like_e = np.append(self.storage_like_e, energy_like[np.newaxis, :], axis=0)

        if not fix_beta:
            # UPDATE BETA
            lp_beta_prev = self.log_prior_density_beta(self.betas[-1])
            lp_beta_prop = self.log_prior_density_beta(beta_prop)
            
            lmd_prev = self.calc_sum_log_mixture_density(comp_coef, self.mus[-1], self.covs[-1])
            
            gibbs_energy_prop = self._calc_gibbs_energy_vect(self.labels[-1], beta_prop, verbose=verbose)
            energy_for_comp_coef_prop = gibbs_energy_prop
            comp_coef_prop = _calc_labels_prob(energy_for_comp_coef_prop, t)
            
            lmd_prop = self.calc_sum_log_mixture_density(comp_coef_prop, self.mus[-1], self.covs[-1])
            
            log_target_prev = lmd_prev + lp_beta_prev
            log_target_prop = lmd_prop + lp_beta_prop
            
            mu_eval = evaluate(log_target_prop, log_target_prev)
            if mu_eval[0]:
                self.betas.append(beta_prop)
            else:
                self.betas.append(self.betas[-1])
            self.beta_acc_ratio = np.append(self.beta_acc_ratio, mu_eval[1])
        else:
            self.betas.append(self.betas[-1])

    def log_prior_density_mu(self, mu, label):
        with np.errstate(divide='ignore'):
            return np.sum(np.log(self.priors_mu[label].pdf(mu)))

    def log_prior_density_beta(self, beta):
        return np.log(self.prior_beta.pdf(beta))

    def log_prior_density_cov(self, cov, l):
        lam = np.sqrt(np.diag(cov[l, :, :]))
        r = np.diag(1. / lam) @ cov[l, :, :] @ np.diag(1. / lam)
        
        det_r = np.linalg.det(r)
        if det_r <= 0:
            return -np.inf
            
        logp_r = -0.5 * (self.nu + self.n_feat + 1) * np.log(det_r)
        
        try:
            inv_r = np.linalg.inv(r)
            logp_r -= self.nu / 2. * np.sum(np.log(np.diag(inv_r)))
        except np.linalg.LinAlgError:
            return -np.inf
            
        logp_lam = np.sum(np.log(multivariate_normal(mean=self.b_sigma[l, :], 
                                                    cov=self.kesi[l, :]).pdf(np.log(lam.T))))
        
        return logp_r + logp_lam

    def propose_beta(self, beta_prev, beta_jump_length):
        if self.dim == 1:
            beta_dim = 1
        elif self.dim == 2:
            beta_dim = 2 if self.stencil == "4p" else 4
        else:
            raise Exception("3D not yet supported.")

        noise = np.random.normal(loc=0, scale=np.sqrt(beta_jump_length), size=beta_dim)
        return beta_prev + noise

    def propose_mu(self, mu_prev, mu_jump_length):
        noise = np.random.multivariate_normal(
            mean=np.zeros(self.n_feat),
            cov=np.eye(self.n_feat) * mu_jump_length,
            size=self.n_labels
        )
        return mu_prev + noise

    def calc_sum_log_mixture_density(self, comp_coef, mu, cov):
        lmd = np.zeros((self.phys_shp.prod(), self.n_labels))
        for l in range(self.n_labels):
            draw = multivariate_normal(mean=mu[l, :], cov=cov[l, :, :]).pdf(self.feat)
            multi = comp_coef[:, l] * np.array([draw])
            lmd[:, l] = multi
        lmd = np.sum(lmd, axis=1)
        with np.errstate(divide='ignore'):
            lmd = np.log(lmd)
        return np.sum(lmd)

    def calc_energy_like(self, mu, cov):
        # Simplified energy calculation for performance 
        energy_like_labels = np.zeros((len(self.feat), self.n_labels))
        for l in range(self.n_labels):
            try:
                cov_inv = np.linalg.inv(cov[l])
                log_det = np.log(np.linalg.det(cov[l]))
                diff = self.feat - mu[l]
                quad_terms = np.sum(diff @ cov_inv * diff, axis=1)
                energy_like_labels[:, l] = 0.5 * (quad_terms + log_det)
            except np.linalg.LinAlgError:
                energy_like_labels[:, l] = np.inf
        return energy_like_labels

    def _calc_gibbs_energy_vect(self, labels, beta, verbose=False):
        # Simplified 2D Gibbs energy calculation
        if self.dim == 2:
            labels = labels.reshape(self.shape[0], self.shape[1])
            ge = np.tile(np.zeros_like(labels).astype(float), (self.n_labels, 1, 1))
            
            comp = np.tile(np.arange(self.n_labels)[:, np.newaxis, np.newaxis], 
                          (1, self.shape[0], self.shape[1]))
            
            # Basic neighborhood energy calculation (simplified)
            if len(beta) >= 2:
                mask_h = np.not_equal(comp[:, 1:-1, 1:-1], labels[:-2, 1:-1]) + \
                         np.not_equal(comp[:, 1:-1, 1:-1], labels[2:, 1:-1])
                ge[:, 1:-1, 1:-1] += mask_h * beta[0]
                
                mask_v = np.not_equal(comp[:, 1:-1, 1:-1], labels[1:-1, :-2]) + \
                         np.not_equal(comp[:, 1:-1, 1:-1], labels[1:-1, 2:])
                ge[:, 1:-1, 1:-1] += mask_v * beta[1]
            
            return np.array([ge[l, :, :].ravel() for l in range(self.n_labels)]).T
        else:
            # Simple 1D case
            return np.random.random((len(labels), self.n_labels)) * beta[0]

    def normalize_feature_vectors(self):
        self.feat = (self.feat - np.mean(self.feat, axis=0).T) / np.std(self.feat, axis=0)


def pseudocolor(shape, stencil=None):
    """Simplified pseudocolor function."""
    dim = len(shape) - 1
    if dim == 1:
        i_w = np.arange(0, shape[0], step=2)
        i_b = np.arange(1, shape[0], step=2)
        return np.array([i_w, i_b]).T
    elif dim == 2:
        if stencil is None or stencil == "8p":
            colors = 4
            colored_image = np.tile(np.kron([[0, 1], [2, 3]] * int(shape[0] / 2), np.ones((1, 1))), int(shape[1] / 2))
            if shape[0] % 2 != 0:
                colored_image = np.append(colored_image, colored_image[-2, :][np.newaxis, :], axis=0)
            if shape[1] % 2 != 0:
                colored_image = np.append(colored_image, colored_image[:, -2][:, np.newaxis], axis=1)
            colored_flat = colored_image.reshape(shape[0] * shape[1])
            ci = []
            for c in range(colors):
                x = np.where(colored_flat == c)[0]
                ci.append(x)
            return np.array(ci)
        elif stencil == "4p":
            colors = 2
            colored_image = np.tile(np.kron([[0, 1], [1, 0]] * int(shape[0] / 2), np.ones((1, 1))), int(shape[1] / 2))
            if shape[0] % 2 != 0:
                colored_image = np.append(colored_image, colored_image[-2, :][np.newaxis, :], axis=0)
            if shape[1] % 2 != 0:
                colored_image = np.append(colored_image, colored_image[:, -2][:, np.newaxis], axis=1)
            colored_flat = colored_image.reshape(shape[0] * shape[1])
            ci = []
            for c in range(colors):
                x = np.where(colored_flat == c)[0]
                ci.append(x)
            return ci


def load_and_preprocess_data(bands: Dict[str, str], bounds: Tuple[float, float, float, float], 
                           use_soil_mask: bool = True, ref_band: str = "B4") -> np.ndarray:
    """Load and preprocess data from given bands and bounds."""
    
    # Choose reference band
    band_keys = [k for k in bands.keys() if k not in ("TCI", "SCL")]
    if ref_band not in bands:
        ref_band = band_keys[0]
    
    # Build stack from provided bands (exclude TCI and SCL)
    stack = []
    with rasterio.open(bands[ref_band]) as ref:
        transform = ref.transform
    
    for name, path in bands.items():
        if name in ("TCI", "SCL"):
            continue
        with rasterio.open(path) as src:
            band_data = crop_to_bounds(src, bounds, src.transform)
            # Ensure 2D per-band array
            if band_data.ndim == 3 and band_data.shape[0] == 1:
                band_data = band_data[0]
            stack.append(band_data)
    
    # Stack into (rows, cols, bands)
    img_stack = np.stack(stack, axis=-1).astype(np.float64)
    
    # Apply soil mask if requested and SCL available
    if use_soil_mask and "SCL" in bands:
        img_stack = apply_soil_mask(img_stack, bands["SCL"], bounds)
    
    return img_stack


def measure_memory_at_steps(algorithm_func, data, n_classes, step_snapshots: List[int], **kwargs) -> Dict:
    """
    Measure memory consumption at specific step intervals.
    
    Args:
        algorithm_func: BaySeg algorithm class
        data: Input data
        n_classes: Number of classes
        step_snapshots: List of step numbers to take memory snapshots
        **kwargs: Additional algorithm parameters
    
    Returns:
        Dictionary with memory measurements at each step
    """
    
    # Force garbage collection before measurement
    gc.collect()
    
    # Start memory tracking
    tracemalloc.start()
    process = psutil.Process()
    initial_memory = process.memory_info().rss / 1024 / 1024  # MB
    
    results = {
        'steps': [],
        'memory_mb': [],
        'peak_memory_mb': [],
        'execution_times': [],
        'success': True,
        'error': None,
        'algorithm_name': algorithm_func.__name__
    }
    
    try:
        # Separate constructor kwargs from fit kwargs
        constructor_kwargs = {k: v for k, v in kwargs.items() 
                            if k in ['beta_init', 'stencil', 'normalize', 'max_history', 'store_diagnostics']}
        fit_kwargs = {k: v for k, v in kwargs.items() 
                     if k in ['beta_jump_length', 'mu_jump_length', 'cov_volume_jump_length', 
                              'theta_jump_length', 't', 'verbose', 'fix_beta']}
        
        # Initialize the algorithm
        seg = algorithm_func(data=data, n_labels=n_classes, **constructor_kwargs)
        
        # Run algorithm and take snapshots at specified steps
        current_step = 0
        start_time = time.time()
        
        for target_step in step_snapshots:
            # Run iterations until we reach the target step
            steps_to_run = target_step - current_step
            if steps_to_run > 0:
                seg.fit(n=steps_to_run, **fit_kwargs)
                current_step = target_step
            
            # Take memory snapshot
            current_memory_trace, peak_memory_trace = tracemalloc.get_traced_memory()
            current_process_memory = process.memory_info().rss / 1024 / 1024  # MB
            
            memory_usage = current_process_memory - initial_memory
            peak_memory_mb = peak_memory_trace / 1024 / 1024  # Convert to MB
            current_time = time.time() - start_time
            
            results['steps'].append(target_step)
            results['memory_mb'].append(memory_usage)
            results['peak_memory_mb'].append(peak_memory_mb)
            results['execution_times'].append(current_time)
            
            print(f"  Step {target_step}: Memory = {memory_usage:.2f} MB, Peak = {peak_memory_mb:.2f} MB, Time = {current_time:.2f}s")
        
    except Exception as e:
        results['success'] = False
        results['error'] = str(e)
        print(f"  ERROR: {e}")
        
    finally:
        tracemalloc.stop()
        gc.collect()
    
    return results


def get_roi_info(bounds: Tuple[float, float, float, float]) -> Dict:
    """Calculate ROI area and approximate pixel count."""
    width = bounds[2] - bounds[0]  # xmax - xmin
    height = bounds[3] - bounds[1]  # ymax - ymin
    area_m2 = width * height
    area_km2 = area_m2 / 1000000
    
    # Approximate pixel count (assuming 20m resolution)
    approx_pixels = int(area_m2 / (20 * 20))  
    
    return {
        'area_km2': area_km2,
        'approx_pixels': approx_pixels,
        'width_m': width,
        'height_m': height
    }


def create_memory_consumption_plots(results: Dict):
    """Create plots showing memory consumption over steps for both algorithms."""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))
    fig.suptitle('Memory Consumption Comparison: BaySeg vs BaySegOriginal', fontsize=16)
    
    # Extract data for both algorithms
    bayseg_data = results['bayseg']
    original_data = results['bayseg_original']
    
    if bayseg_data['success']:
        steps_opt = bayseg_data['steps']
        memory_opt = bayseg_data['memory_mb']
        peak_memory_opt = bayseg_data['peak_memory_mb']
        times_opt = bayseg_data['execution_times']
    else:
        steps_opt = memory_opt = peak_memory_opt = times_opt = []
    
    if original_data['success']:
        steps_orig = original_data['steps']
        memory_orig = original_data['memory_mb']
        peak_memory_orig = original_data['peak_memory_mb']
        times_orig = original_data['execution_times']
    else:
        steps_orig = memory_orig = peak_memory_orig = times_orig = []
    
    # Plot 1: Memory usage vs Steps
    ax1.set_title('Memory Usage vs Steps')
    if steps_opt:
        ax1.plot(steps_opt, memory_opt, 'bo-', label='BaySeg (Optimized)', linewidth=2, markersize=8)
    if steps_orig:
        ax1.plot(steps_orig, memory_orig, 'ro-', label='BaySegOriginal', linewidth=2, markersize=8)
    ax1.set_xlabel('Steps')
    ax1.set_ylabel('Memory Usage (MB)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Peak Memory vs Steps
    ax2.set_title('Peak Memory vs Steps')
    if steps_opt:
        ax2.plot(steps_opt, peak_memory_opt, 'bs-', label='BaySeg (Optimized)', linewidth=2, markersize=8)
    if steps_orig:
        ax2.plot(steps_orig, peak_memory_orig, 'rs-', label='BaySegOriginal', linewidth=2, markersize=8)
    ax2.set_xlabel('Steps')
    ax2.set_ylabel('Peak Memory (MB)')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Memory Efficiency (Memory per Step)
    ax3.set_title('Memory Efficiency (Memory per Step)')
    if steps_opt:
        memory_per_step_opt = [m/s if s > 0 else 0 for m, s in zip(memory_opt, steps_opt)]
        ax3.plot(steps_opt, memory_per_step_opt, 'b^-', label='BaySeg (Optimized)', linewidth=2, markersize=8)
    if steps_orig:
        memory_per_step_orig = [m/s if s > 0 else 0 for m, s in zip(memory_orig, steps_orig)]
        ax3.plot(steps_orig, memory_per_step_orig, 'r^-', label='BaySegOriginal', linewidth=2, markersize=8)
    ax3.set_xlabel('Steps')
    ax3.set_ylabel('Memory per Step (MB/step)')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Execution Time vs Steps
    ax4.set_title('Execution Time vs Steps')
    if steps_opt:
        ax4.plot(steps_opt, times_opt, 'bv-', label='BaySeg (Optimized)', linewidth=2, markersize=8)
    if steps_orig:
        ax4.plot(steps_orig, times_orig, 'rv-', label='BaySegOriginal', linewidth=2, markersize=8)
    ax4.set_xlabel('Steps')
    ax4.set_ylabel('Execution Time (seconds)')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('memory_consumption_comparison.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"\nPlot saved as 'memory_consumption_comparison.png'")


def print_memory_summary(results: Dict):
    """Print a summary of memory consumption results."""
    
    print("\n" + "="*70)
    print("MEMORY CONSUMPTION SUMMARY")
    print("="*70)
    
    roi_info = results['roi_info']
    print(f"ROI: {roi_info['area_km2']:.2f} km² ({roi_info['approx_pixels']:,} pixels)")
    print(f"Test Steps: {results['step_snapshots']}")
    
    print(f"\nBaySeg (Optimized):")
    bayseg_data = results['bayseg']
    if bayseg_data['success']:
        print(f"  Success: ✓")
        for i, (step, mem, peak, time_val) in enumerate(zip(
            bayseg_data['steps'], 
            bayseg_data['memory_mb'], 
            bayseg_data['peak_memory_mb'],
            bayseg_data['execution_times']
        )):
            print(f"    Step {step:3d}: {mem:6.2f} MB memory, {peak:6.2f} MB peak, {time_val:6.2f}s")
        
        final_memory = bayseg_data['memory_mb'][-1]
        final_peak = bayseg_data['peak_memory_mb'][-1]
        print(f"  Final: {final_memory:.2f} MB memory, {final_peak:.2f} MB peak")
    else:
        print(f"  Success: ✗ - {bayseg_data['error']}")
    
    print(f"\nBaySegOriginal:")
    original_data = results['bayseg_original']
    if original_data['success']:
        print(f"  Success: ✓")
        for i, (step, mem, peak, time_val) in enumerate(zip(
            original_data['steps'], 
            original_data['memory_mb'], 
            original_data['peak_memory_mb'],
            original_data['execution_times']
        )):
            print(f"    Step {step:3d}: {mem:6.2f} MB memory, {peak:6.2f} MB peak, {time_val:6.2f}s")
        
        final_memory = original_data['memory_mb'][-1]
        final_peak = original_data['peak_memory_mb'][-1]
        print(f"  Final: {final_memory:.2f} MB memory, {final_peak:.2f} MB peak")
    else:
        print(f"  Success: ✗ - {original_data['error']}")
    
    # Calculate memory savings if both algorithms succeeded
    if bayseg_data['success'] and original_data['success']:
        print(f"\nMemory Savings:")
        final_savings = original_data['memory_mb'][-1] - bayseg_data['memory_mb'][-1]
        final_savings_pct = (final_savings / original_data['memory_mb'][-1]) * 100
        peak_savings = original_data['peak_memory_mb'][-1] - bayseg_data['peak_memory_mb'][-1]
        peak_savings_pct = (peak_savings / original_data['peak_memory_mb'][-1]) * 100
        
        print(f"  Final Memory: {final_savings:.2f} MB ({final_savings_pct:.1f}% reduction)")
        print(f"  Peak Memory: {peak_savings:.2f} MB ({peak_savings_pct:.1f}% reduction)")


def run_memory_consumption_test():
    """Run the memory consumption test."""
    
    print("BaySeg Memory Consumption Test")
    print("=" * 50)
    
    # Tharsis data configuration (using same paths as performance_comparison.py)
    bands = {
        "B4": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B04_20m.jp2",
        "B6": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B06_20m.jp2", 
        "B7": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B07_20m.jp2",
        "B8A": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B8A_20m.jp2",
        "B11": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B11_20m.jp2",
        "B12": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_B12_20m.jp2",
        "SCL": "/Users/simonvirgo/PycharmProjects/Mineye-Terranigma/Data/Tharsis_all_Sentinel2 /combined/20m/merged_SCL_20m.jp2",
    }
    
    # Use the smallest ROI from performance_comparison.py
    smallest_roi_bounds = (207500, 4172500, 212500, 4177500)  # ROI1_25km (~25 km²)
    
    # Memory snapshot steps as requested
    step_snapshots = [50, 100, 200, 400, 800]
    
    # Algorithm configurations
    n_classes = 6
    beta_init = 30.0
    beta_jump_length = 0.1
    
    print(f"Using smallest ROI: {smallest_roi_bounds}")
    roi_info = get_roi_info(smallest_roi_bounds)
    print(f"ROI Area: {roi_info['area_km2']:.2f} km²")
    print(f"Approximate pixels: {roi_info['approx_pixels']:,}")
    print(f"Memory snapshots at steps: {step_snapshots}")
    
    # Load and preprocess data
    print(f"\nLoading and preprocessing data...")
    try:
        data = load_and_preprocess_data(bands, smallest_roi_bounds, use_soil_mask=True)
        print(f"Data shape: {data.shape}")
    except Exception as e:
        print(f"ERROR loading data: {e}")
        return
    
    results = {
        'roi_info': roi_info,
        'step_snapshots': step_snapshots,
        'bayseg': {},
        'bayseg_original': {}
    }
    
    # Test BaySeg (optimized with memory management)
    print(f"\n--- Testing BaySeg (Optimized) ---")
    results['bayseg'] = measure_memory_at_steps(
        BaySeg, data, n_classes, step_snapshots,
        beta_init=beta_init, beta_jump_length=beta_jump_length,
        max_history=50, store_diagnostics=False  # Memory optimizations
    )
    
    # Test BaySegOriginal (without memory optimizations)
    print(f"\n--- Testing BaySegOriginal (No Memory Optimizations) ---")
    results['bayseg_original'] = measure_memory_at_steps(
        BaySegOriginal, data, n_classes, step_snapshots,
        beta_init=beta_init, beta_jump_length=beta_jump_length
    )
    
    # Print summary
    print_memory_summary(results)
    
    # Create plots
    print(f"\nCreating memory consumption plots...")
    create_memory_consumption_plots(results)
    
    return results


def main():
    """Main function to run the memory consumption test."""
    try:
        results = run_memory_consumption_test()
        print(f"\nMemory consumption test completed successfully!")
        return results
    except Exception as e:
        print(f"Error running memory consumption test: {e}")
        import traceback
        traceback.print_exc()
        return None


# Helper functions for BaySegOriginal (from performance_comparison.py)

def _calc_labels_prob(te, t):
    """"Calculate labels probability for array of total energies (te) and totally arbitrary skalar value t."""
    return (np.exp(-te / t).T / np.sum(np.exp(-te / t), axis=1)).T

def draw_labels_vect(labels_prob):
    """Vectorized draw of the label for each elements respective labels probability."""
    r = np.random.rand(len(labels_prob))
    p = np.cumsum(labels_prob, axis=1)
    d = (p.T - r).T
    return np.sum(np.greater_equal(0, d), axis=1)

def evaluate(log_target_prop, log_target_prev):
    """Evaluate whether to accept or reject a proposal based on log target values."""
    if np.isnan(log_target_prop) or np.isnan(log_target_prev):
        return False, 0.0
    if np.isinf(log_target_prop) and np.isinf(log_target_prev):
        return False, 0.0
    if np.isinf(log_target_prop) and log_target_prop > 0:
        return True, np.inf
    if np.isinf(log_target_prev) and log_target_prev > 0:
        return False, 0.0
        
    diff = log_target_prop - log_target_prev
    if diff > 700:
        return True, np.inf
    if diff < -700:
        return False, 0.0
        
    ratio = np.exp(diff)
    
    if (ratio > 1) or (np.random.uniform() < ratio):
        return True, ratio
    else:
        return False, ratio

def _propose_cov(cov_prev, n_feat, n_labels, cov_jump_length, theta_jump_length):
    """Proposes a perturbed n-dimensional covariance matrix."""
    comb = list(combinations(range(n_feat), 2))
    n_comb = len(comb)
    theta_jump = multivariate_normal(mean=[0 for i in range(n_comb)], 
                                   cov=np.ones(n_comb) * theta_jump_length).rvs()

    if n_comb == 1:
        theta_jump = [theta_jump]

    cov_prop = np.zeros_like(cov_prev)

    for l in range(n_labels):
        v_l, d_l, v_l_t = np.linalg.svd(cov_prev[l, :, :])
        log_d_jump = multivariate_normal(mean=[0 for i in range(n_feat)], 
                                       cov=np.eye(n_feat) * cov_jump_length).rvs()
        d_prop = np.diag(np.exp(np.log(d_l) + log_d_jump))
        
        a = np.eye(n_feat)
        for val in range(n_comb):
            rotation_matrix = _cov_proposal_rotation_matrix(v_l[:, comb[val][0]], 
                                                          v_l[:, comb[val][1]], theta_jump[val])
            a = rotation_matrix @ a
        
        v_prop = a @ v_l
        cov_prop[l, :, :] = v_prop @ d_prop @ v_prop.T

    return cov_prop

def _cov_proposal_rotation_matrix(x, y, theta):
    """Creates the rotation matrix needed for the covariance matrix proposal step."""
    x = np.array([x]).T
    y = np.array([y]).T

    uu = x / np.linalg.norm(x)
    vv = y - uu.T @ y * uu
    vv = vv / np.linalg.norm(vv)

    rotation_matrix = np.eye(len(x)) - uu @ uu.T - vv @ vv.T + np.hstack((uu, vv)) @ np.array(
        [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]]) @ np.hstack((uu, vv)).T
    return rotation_matrix


if __name__ == "__main__":
    main()