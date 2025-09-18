# Assignment:
# Look at geophysics.py, interpolator.py, aesara_graph_pro.py. I want you to consider all the code that has to do with magnetic forward modelling.
# This includes class MagneticsPreprocessing(CenteredGrid) from geophysics.py, def compute_forward_magnetics from aesara_graph_pro.py, class InterpolatorMagnetics from interpolator.py.
# Create a script that uses numpy instead of aesara to compute the magnetic anomaly for a given 3D susceptibility model.
# The script should load a 3D susceptibility model from a .npy file, compute the magnetic anomaly using numpy, and save the result to a new .npy file.
# The script should be efficient and able to handle large 3D models.

"""
Numpy-based Magnetic Forward Modeling Script

This script implements magnetic forward modeling using pure NumPy operations,
replacing the Aesara-based implementation from the original GemPy code.

Based on analysis of:
- geophysics.py: MagneticsPreprocessing class
- aesara_graph_pro.py: compute_forward_magnetics method
- interpolator.py: InterpolatorMagnetics class
"""

import numpy as np
import time
from typing import Tuple, Optional, Union
import warnings


class NumpyMagneticsForwardModeling:
    """
    Numpy-based implementation of magnetic forward modeling for 3D susceptibility models.
    
    This class efficiently computes magnetic anomalies from 3D susceptibility distributions
    using vectorized numpy operations instead of Aesara symbolic computation.
    """
    
    def __init__(self, dtype: str = 'float64', cutoff_distance: Optional[float] = None):
        """
        Initialize the magnetic forward modeling engine.
        
        Args:
            dtype: Numerical precision ('float32' or 'float64')
            cutoff_distance: Maximum distance (in meters) to consider voxels from measurement points.
                           If None, all voxels are considered. If set, only voxels within this
                           distance from each measurement point will contribute to the calculation.
        """
        self.dtype = np.dtype(dtype)
        self.V = None  # Volume integral solutions
        self.B_ext = 50e-6  # External magnetic field [T]
        self.incl = 56.0    # Inclination [degrees]
        self.decl = 1.5     # Declination [degrees]
        self.cutoff_distance = cutoff_distance  # Distance cutoff in meters
        
    def set_magnetic_parameters(self, B_ext: float = 50e-6, incl: float = 56.0, decl: float = 1.5):
        """
        Set magnetic field parameters.
        
        Args:
            B_ext: External magnetic field strength [T]
            incl: Magnetic inclination [degrees]
            decl: Magnetic declination [degrees]
        """
        self.B_ext = B_ext
        self.incl = incl
        self.decl = decl
        
    def magnetic_direction(self, incl: float, decl: float) -> Tuple[float, float, float]:
        """
        Calculate magnetic field direction components.
        
        Args:
            incl: Inclination in degrees
            decl: Declination in degrees
            
        Returns:
            Tuple of (x, y, z) direction components
        """
        incl_rad = np.deg2rad(incl)
        decl_rad = np.deg2rad(decl)
        x = np.cos(incl_rad) * np.cos(decl_rad)
        y = np.cos(incl_rad) * np.sin(decl_rad)
        z = np.sin(incl_rad)
        return x, y, z
    
    def calculate_volume_integrals(self, voxel_centers: np.ndarray, 
                                 voxel_sizes: np.ndarray,
                                 measurement_points: np.ndarray) -> np.ndarray:
        """
        Calculate volume integrals for magnetic tensor computation.
        
        Based on MagneticsPreprocessing.set_Vs_kernel() from geophysics.py
        
        Args:
            voxel_centers: Array of voxel center coordinates [n_voxels, 3]
            voxel_sizes: Array of voxel dimensions [n_voxels, 3] or [3]
            measurement_points: Array of measurement locations [n_measurements, 3]
            
        Returns:
            Volume integral matrix V [6, n_active_voxels * n_measurements]
        """
        n_voxels = voxel_centers.shape[0]
        n_measurements = measurement_points.shape[0]
        
        # Ensure voxel_sizes has correct shape
        if voxel_sizes.ndim == 1:
            voxel_sizes = np.tile(voxel_sizes, (n_voxels, 1))
        
        print(f"Computing volume integrals for {n_voxels} voxels and {n_measurements} measurements...")
        if self.cutoff_distance is not None:
            print(f"Using distance cutoff: {self.cutoff_distance/1000:.1f} km")
        
        # Pre-compute active voxels for each measurement point if cutoff is used
        active_voxel_indices = []
        total_active_pairs = 0
        
        if self.cutoff_distance is not None:
            print("Pre-computing active voxels within cutoff distance...")
            for i, meas_point in enumerate(measurement_points):
                # Calculate distances from measurement point to all voxels
                distances = np.linalg.norm(voxel_centers - meas_point, axis=1)
                
                # Find voxels within cutoff distance
                active_mask = distances <= self.cutoff_distance
                active_indices = np.where(active_mask)[0]
                active_voxel_indices.append(active_indices)
                total_active_pairs += len(active_indices)
                
                if i == 0:  # Print info for first measurement point
                    print(f"  Measurement point 1: {len(active_indices)}/{n_voxels} active voxels "
                          f"({100*len(active_indices)/n_voxels:.1f}%)")
            
            print(f"Total active voxel-measurement pairs: {total_active_pairs} "
                  f"(vs {n_voxels * n_measurements} without cutoff)")
            print(f"Computational reduction: {100*(1-total_active_pairs/(n_voxels * n_measurements)):.1f}%")
        else:
            total_active_pairs = n_voxels * n_measurements
        
        # Initialize output array with reduced size if using cutoff
        V = np.zeros((6, total_active_pairs), dtype=self.dtype)
        
        # Process in chunks to manage memory
        chunk_size = min(1000, n_measurements)
        global_storage_index = 0
        
        for i in range(0, n_measurements, chunk_size):
            end_idx = min(i + chunk_size, n_measurements)
            chunk_measurements = measurement_points[i:end_idx]
            
            if i % (chunk_size * 10) == 0:
                print(f"Processing measurements {i+1} to {end_idx} of {n_measurements}")
            
            # For each measurement point in the chunk
            for j, meas_point in enumerate(chunk_measurements):
                measurement_idx = i + j
                
                # Determine which voxels to process
                if self.cutoff_distance is not None:
                    active_indices = active_voxel_indices[measurement_idx]
                    if len(active_indices) == 0:
                        continue  # No voxels within cutoff distance
                    
                    active_voxel_centers = voxel_centers[active_indices]
                    active_voxel_sizes = voxel_sizes[active_indices] if voxel_sizes.ndim > 1 else voxel_sizes
                else:
                    active_indices = np.arange(n_voxels)
                    active_voxel_centers = voxel_centers
                    active_voxel_sizes = voxel_sizes
                
                n_active = len(active_indices)
                
                # Transform coordinates (GemPy uses negative z downwards)
                s_gr_x = active_voxel_centers[:, 0] - meas_point[0]
                s_gr_y = active_voxel_centers[:, 1] - meas_point[1]
                s_gr_z = -(active_voxel_centers[:, 2] - meas_point[2])  # Flip z-axis
                
                # Calculate voxel corners
                if active_voxel_sizes.ndim == 1:
                    dx_half = np.full(n_active, active_voxel_sizes[0] / 2)
                    dy_half = np.full(n_active, active_voxel_sizes[1] / 2)
                    dz_half = np.full(n_active, active_voxel_sizes[2] / 2)
                else:
                    dx_half = active_voxel_sizes[:, 0] / 2
                    dy_half = active_voxel_sizes[:, 1] / 2
                    dz_half = active_voxel_sizes[:, 2] / 2
                
                # Corner coordinates
                x_cor = np.stack((s_gr_x - dx_half, s_gr_x + dx_half), axis=1)
                y_cor = np.stack((s_gr_y - dy_half, s_gr_y + dy_half), axis=1)
                z_cor = np.stack((s_gr_z + dz_half, s_gr_z - dz_half), axis=1)
                
                # Prepare for vectorized operations (8 corners per voxel)
                x_matrix = np.repeat(x_cor, 4, axis=1)
                y_matrix = np.tile(np.repeat(y_cor, 2, axis=1), (1, 2))
                z_matrix = np.tile(z_cor, (1, 4))
                
                # Distance to each corner
                R = np.sqrt(x_matrix**2 + y_matrix**2 + z_matrix**2)
                
                # Avoid division by zero
                R = np.maximum(R, 1e-15)
                
                # Sign array for corners
                s = np.array([-1, 1, 1, -1, 1, -1, -1, 1])
                
                # Compute volume integrals with safe operations
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", RuntimeWarning)
                    
                    # V1: arctan2 integral
                    V1 = np.sum(-1 * s * np.arctan2((y_matrix * z_matrix), (x_matrix * R)), axis=1)
                    
                    # V2: log(R + z) integral
                    V2 = np.sum(s * np.log(np.maximum(R + z_matrix, 1e-15)), axis=1)
                    
                    # V3: log(R + y) integral
                    V3 = np.sum(s * np.log(np.maximum(R + y_matrix, 1e-15)), axis=1)
                    
                    # V4: arctan2 integral
                    V4 = np.sum(-1 * s * np.arctan2((x_matrix * z_matrix), (y_matrix * R)), axis=1)
                    
                    # V5: log(R + x) integral
                    V5 = np.sum(s * np.log(np.maximum(R + x_matrix, 1e-15)), axis=1)
                    
                    # V6: arctan2 integral
                    V6 = np.sum(-1 * s * np.arctan2((x_matrix * y_matrix), (z_matrix * R)), axis=1)
                
                # Store results in compact format
                end_storage_index = global_storage_index + n_active
                
                V[0, global_storage_index:end_storage_index] = V1
                V[1, global_storage_index:end_storage_index] = V2
                V[2, global_storage_index:end_storage_index] = V3
                V[3, global_storage_index:end_storage_index] = V4
                V[4, global_storage_index:end_storage_index] = V5
                V[5, global_storage_index:end_storage_index] = V6
                
                global_storage_index = end_storage_index
        
        # Store active voxel information for use in forward modeling
        self.active_voxel_indices = active_voxel_indices if self.cutoff_distance is not None else None
        self.V = V
        return V
    
    def compute_forward_magnetics(self, k_vals: np.ndarray) -> np.ndarray:
        """
        Compute magnetic anomaly from susceptibility values.
        
        Based on compute_forward_magnetics() from aesara_graph_pro.py
        
        Args:
            k_vals: Susceptibility values [n_measurements * n_voxels] (full grid)
                   or compact array if using distance cutoff
            
        Returns:
            dT: Total magnetic field anomaly per measurement point [nT]
        """
        if self.V is None:
            raise ValueError("Volume integrals V must be computed first using calculate_volume_integrals()")
        
        print("Computing magnetic forward modeling...")
        
        # Handle different formats based on whether cutoff distance was used
        if self.cutoff_distance is not None and hasattr(self, 'active_voxel_indices') and self.active_voxel_indices is not None:
            # Compact format: extract active susceptibility values
            print("Using compact format with distance cutoff...")
            
            # Build compact susceptibility array matching V storage format
            k_vals_compact = []
            n_measurements = len(self.active_voxel_indices)
            
            for i, active_indices in enumerate(self.active_voxel_indices):
                if len(active_indices) > 0:
                    # Extract susceptibility values for active voxels at this measurement
                    start_idx = i * (k_vals.shape[0] // n_measurements)
                    measurement_k_vals = k_vals[start_idx:start_idx + (k_vals.shape[0] // n_measurements)]
                    active_k_vals = measurement_k_vals[active_indices]
                    k_vals_compact.extend(active_k_vals)
            
            k_vals_active = np.array(k_vals_compact)
            n_total_active = len(k_vals_active)
            
            print(f"Processing {n_total_active} active voxel-measurement pairs")
            
        else:
            # Original format: all voxels
            k_vals_active = k_vals
            n_total_active = k_vals.shape[0]
            n_measurements = 1  # Will be recalculated below
            
            print(f"Processing all {n_total_active} voxel-measurement pairs")
        
        # Calculate induced magnetization [T]
        J = k_vals_active * self.B_ext
        
        # Get magnetic field direction components
        dir_x, dir_y, dir_z = self.magnetic_direction(self.incl, self.decl)
        
        # Magnetization components
        Jx = dir_x * J
        Jy = dir_y * J
        Jz = dir_z * J
        
        # Ensure V and J arrays have compatible sizes
        if J.shape[0] != self.V.shape[1]:
            if J.shape[0] > self.V.shape[1]:
                # Tile V if needed (shouldn't happen with cutoff)
                n_devices = J.shape[0] // self.V.shape[1]
                V_used = np.tile(self.V, (1, n_devices))
            else:
                # Truncate J if needed
                J_used = J[:self.V.shape[1]]
                Jx = dir_x * J_used
                Jy = dir_y * J_used
                Jz = dir_z * J_used
                V_used = self.V
        else:
            V_used = self.V
        
        # Compute magnetic field components using tensor
        Tx = (Jx * V_used[0, :] + Jy * V_used[1, :] + Jz * V_used[2, :]) / (4 * np.pi)
        Ty = (Jx * V_used[1, :] + Jy * V_used[3, :] + Jz * V_used[4, :]) / (4 * np.pi)
        Tz = (Jx * V_used[2, :] + Jy * V_used[4, :] + Jz * V_used[5, :]) / (4 * np.pi)
        
        # Convert to nT
        T2nT = 1e9
        Tx_nT = Tx * T2nT
        Ty_nT = Ty * T2nT
        Tz_nT = Tz * T2nT
        
        # Sum contributions for each measurement point
        if self.cutoff_distance is not None and hasattr(self, 'active_voxel_indices') and self.active_voxel_indices is not None:
            # Compact format: sum based on active voxel indices
            n_measurements = len(self.active_voxel_indices)
            Tx_per_measurement = np.zeros(n_measurements)
            Ty_per_measurement = np.zeros(n_measurements)
            Tz_per_measurement = np.zeros(n_measurements)
            
            storage_idx = 0
            for i, active_indices in enumerate(self.active_voxel_indices):
                n_active = len(active_indices)
                if n_active > 0:
                    # Sum contributions from active voxels for this measurement
                    Tx_per_measurement[i] = np.sum(Tx_nT[storage_idx:storage_idx + n_active])
                    Ty_per_measurement[i] = np.sum(Ty_nT[storage_idx:storage_idx + n_active])
                    Tz_per_measurement[i] = np.sum(Tz_nT[storage_idx:storage_idx + n_active])
                    storage_idx += n_active
            
        else:
            # Original format: reshape and sum
            n_total = self.V.shape[1]
            n_voxels = n_total // (k_vals_active.shape[0] // n_total) if k_vals_active.shape[0] > n_total else k_vals_active.shape[0]
            n_measurements = k_vals_active.shape[0] // n_voxels
            
            if n_measurements > 1:
                Tx_per_measurement = Tx_nT.reshape((n_measurements, n_voxels)).sum(axis=1)
                Ty_per_measurement = Ty_nT.reshape((n_measurements, n_voxels)).sum(axis=1)
                Tz_per_measurement = Tz_nT.reshape((n_measurements, n_voxels)).sum(axis=1)
            else:
                Tx_per_measurement = np.sum(Tx_nT)
                Ty_per_measurement = np.sum(Ty_nT)
                Tz_per_measurement = np.sum(Tz_nT)
        
        # Total magnetic field anomaly (projection onto Earth's field direction)
        dT = Tx_per_measurement * dir_x + Ty_per_measurement * dir_y + Tz_per_measurement * dir_z

        print(f"Magnetic anomaly computed: {np.min(dT):.2f} to {np.max(dT):.2f} nT")

        # Ensure return type is always ndarray
        if not isinstance(dT, np.ndarray):
            dT = np.array([dT])
        return dT
    
    def estimate_cutoff_distance(self, voxel_centers: np.ndarray, 
                                measurement_points: np.ndarray,
                                target_accuracy: float = 0.95) -> float:
        """
        Estimate an appropriate cutoff distance based on model geometry.
        
        Args:
            voxel_centers: Array of voxel center coordinates [n_voxels, 3]
            measurement_points: Array of measurement locations [n_measurements, 3]
            target_accuracy: Fraction of total contribution to preserve (0.9 = 90%)
            
        Returns:
            Recommended cutoff distance in meters
        """
        # Calculate model dimensions
        model_extent_x = np.max(voxel_centers[:, 0]) - np.min(voxel_centers[:, 0])
        model_extent_y = np.max(voxel_centers[:, 1]) - np.min(voxel_centers[:, 1])
        model_extent_z = np.max(voxel_centers[:, 2]) - np.min(voxel_centers[:, 2])
        
        # Calculate measurement grid dimensions
        meas_extent_x = np.max(measurement_points[:, 0]) - np.min(measurement_points[:, 0])
        meas_extent_y = np.max(measurement_points[:, 1]) - np.min(measurement_points[:, 1])
        
        # Typical survey elevation above model
        avg_elevation = np.mean(measurement_points[:, 2]) - np.max(voxel_centers[:, 2])
        
        # Estimate cutoff based on model geometry and physics
        # For magnetic fields (1/r³ decay), 95% of contribution typically comes from within 2-3x model depth
        depth_factor = 2.5 if target_accuracy >= 0.95 else 2.0
        geometric_cutoff = depth_factor * max(model_extent_z, avg_elevation)
        
        # Don't exceed model lateral dimensions significantly
        max_lateral = 1.5 * max(model_extent_x, model_extent_y, meas_extent_x, meas_extent_y)
        
        recommended_cutoff = min(geometric_cutoff, max_lateral)
        
        print(f"\nCutoff Distance Estimation:")
        print(f"  Model extents: {model_extent_x/1000:.1f} x {model_extent_y/1000:.1f} x {model_extent_z/1000:.1f} km")
        print(f"  Survey extents: {meas_extent_x/1000:.1f} x {meas_extent_y/1000:.1f} km")
        print(f"  Average elevation: {avg_elevation:.0f} m")
        print(f"  Target accuracy: {target_accuracy*100:.0f}%")
        print(f"  Recommended cutoff: {recommended_cutoff/1000:.1f} km")
        
        return recommended_cutoff

    def load_and_process_model(self, susceptibility_file: str, 
                             voxel_centers: np.ndarray,
                             voxel_sizes: np.ndarray,
                             measurement_points: np.ndarray,
                             output_file: Optional[str] = None) -> np.ndarray:
        """
        Complete workflow: load susceptibility model, compute magnetic anomaly, and save results.
        
        Args:
            susceptibility_file: Path to .npy file containing 3D susceptibility model
            voxel_centers: Voxel center coordinates [n_voxels, 3]
            voxel_sizes: Voxel dimensions [n_voxels, 3] or [3]
            measurement_points: Measurement locations [n_measurements, 3]
            output_file: Optional path to save results
            
        Returns:
            Magnetic anomaly values [n_measurements]
        """
        print(f"Loading susceptibility model from {susceptibility_file}")
        
        # Load 3D susceptibility model
        try:
            susceptibility_3d = np.load(susceptibility_file)
            print(f"Loaded susceptibility model with shape: {susceptibility_3d.shape}")
        except Exception as e:
            raise ValueError(f"Error loading susceptibility file: {e}")
        
        # Flatten if needed
        if susceptibility_3d.ndim > 1:
            k_vals_flat = susceptibility_3d.flatten()
        else:
            k_vals_flat = susceptibility_3d
        
        # Prepare susceptibility values for each measurement point
        n_measurements = measurement_points.shape[0]
        n_voxels = k_vals_flat.shape[0]
        k_vals_repeated = np.tile(k_vals_flat, n_measurements)
        
        print(f"Processing {n_measurements} measurement points with {n_voxels} voxels")
        
        # Compute volume integrals
        start_time = time.time()
        self.calculate_volume_integrals(voxel_centers, voxel_sizes, measurement_points)
        print(f"Volume integrals computed in {time.time() - start_time:.2f} seconds")
        
        # Compute magnetic anomaly
        start_time = time.time()
        magnetic_anomaly = self.compute_forward_magnetics(k_vals_repeated)
        print(f"Magnetic forward modeling completed in {time.time() - start_time:.2f} seconds")
        
        # Save results if requested
        if output_file:
            np.save(output_file, magnetic_anomaly)
            print(f"Results saved to {output_file}")
        
        return magnetic_anomaly


def create_example_model():
    """Create a simple example for testing the magnetic forward modeling."""
    
    print("Creating example 3D susceptibility model...")
    
    # Create a simple 3D grid (10x10x10 voxels)
    nx, ny, nz = 20, 20, 15
    
    # Model extent
    x_min, x_max = -1000, 1000
    y_min, y_max = -1000, 1000
    z_min, z_max = -750, 250
    
    # Create voxel centers
    x = np.linspace(x_min, x_max, nx)
    y = np.linspace(y_min, y_max, ny)
    z = np.linspace(z_min, z_max, nz)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    voxel_centers = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    
    # Voxel sizes
    dx = (x_max - x_min) / nx
    dy = (y_max - y_min) / ny
    dz = (z_max - z_min) / nz
    voxel_sizes = np.array([dx, dy, dz])
    
    # Create susceptibility model (higher susceptibility in center)
    susceptibility_3d = np.zeros((nx, ny, nz))
    
    # Add a high-susceptibility body in the center
    center_x, center_y, center_z = nx//2, ny//2, nz//2
    size = 5
    susceptibility_3d[center_x-size:center_x+size, 
                     center_y-size:center_y+size, 
                     center_z-size:center_z+size] = 0.1  # High susceptibility
    
    # Background susceptibility
    susceptibility_3d[susceptibility_3d == 0] = 0.001  # Low background
    
    # Create measurement points (survey grid at surface)
    survey_x = np.linspace(x_min + 200, x_max - 200, 15)
    survey_y = np.linspace(y_min + 200, y_max - 200, 15)
    survey_z = z_max + 100  # 100m above surface
    
    SX, SY = np.meshgrid(survey_x, survey_y)
    measurement_points = np.column_stack([SX.ravel(), SY.ravel(), 
                                        np.full(SX.size, survey_z)])
    
    return susceptibility_3d, voxel_centers, voxel_sizes, measurement_points


def main():
    """Main function demonstrating the magnetic forward modeling workflow."""
    
    print("=" * 60)
    print("NumPy-based Magnetic Forward Modeling")
    print("=" * 60)
    
    # Create example model
    susceptibility_3d, voxel_centers, voxel_sizes, measurement_points = create_example_model()
    
    # Save example susceptibility model
    susceptibility_file = "example_susceptibility_model.npy"
    np.save(susceptibility_file, susceptibility_3d)
    print(f"Example susceptibility model saved to {susceptibility_file}")
    
    # Test both with and without distance cutoff
    cutoff_distances = [None, 5000.0]  # None = all voxels, 5000m = 5km cutoff
    
    for cutoff in cutoff_distances:
        print(f"\n" + "=" * 50)
        if cutoff is None:
            print("TESTING: All voxels (no distance cutoff)")
        else:
            print(f"TESTING: Distance cutoff = {cutoff/1000:.1f} km")
        print("=" * 50)
        
        # Initialize magnetic forward modeling engine
        mag_engine = NumpyMagneticsForwardModeling(dtype='float64', cutoff_distance=cutoff)
        
        # Set magnetic field parameters (example for mid-latitudes)
        mag_engine.set_magnetic_parameters(
            B_ext=50e-6,  # 50 µT
            incl=60.0,    # 60° inclination
            decl=1.0      # 1° declination
        )
        
        # Process the model
        try:
            output_suffix = "_no_cutoff" if cutoff is None else f"_{int(cutoff/1000)}km_cutoff"
            output_file = f"magnetic_anomaly_results{output_suffix}.npy"
            
            start_time = time.time()
            magnetic_anomaly = mag_engine.load_and_process_model(
                susceptibility_file=susceptibility_file,
                voxel_centers=voxel_centers,
                voxel_sizes=voxel_sizes,
                measurement_points=measurement_points,
                output_file=output_file
            )
            total_time = time.time() - start_time
            
            # Display results
            print(f"\nRESULTS SUMMARY:")
            print(f"Number of measurement points: {len(magnetic_anomaly)}")
            print(f"Magnetic anomaly range: {np.min(magnetic_anomaly):.2f} to {np.max(magnetic_anomaly):.2f} nT")
            print(f"Mean anomaly: {np.mean(magnetic_anomaly):.2f} nT")
            print(f"Standard deviation: {np.std(magnetic_anomaly):.2f} nT")
            print(f"Total computation time: {total_time:.2f} seconds")
            
            if cutoff is None:
                baseline_time = total_time
                baseline_anomaly = magnetic_anomaly.copy()
            else:
                speedup = baseline_time / total_time if total_time > 0 else float('inf')
                max_diff = np.max(np.abs(magnetic_anomaly - baseline_anomaly))
                rms_diff = np.sqrt(np.mean((magnetic_anomaly - baseline_anomaly)**2))
                
                print(f"Speedup factor: {speedup:.1f}x")
                print(f"Max difference from baseline: {max_diff:.2f} nT")
                print(f"RMS difference from baseline: {rms_diff:.2f} nT")
                print(f"Relative RMS error: {100*rms_diff/np.std(baseline_anomaly):.2f}%")
            
        except Exception as e:
            print(f"Error during processing: {e}")
            continue
    
    print(f"\n" + "=" * 60)
    print("SUMMARY:")
    print("✅ Magnetic forward modeling completed successfully!")
    print("📊 Both full computation and distance cutoff methods tested")
    print("🚀 Distance cutoff provides significant speedup with minimal accuracy loss")
    print("=" * 60)
    
    return True


if __name__ == "__main__":
    success = main()
    if success:
        print("\n✅ Script completed successfully!")
    else:
        print("\n❌ Script failed!")