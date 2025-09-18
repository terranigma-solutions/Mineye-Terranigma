"""
This file is part of gempy (NumPy version).

Function that performs the equivalent of Aesara sparse operations using NumPy/SciPy.
"""

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve as sparse_solve

# Check if skcuda is installed (optional, not used in NumPy version)
try:
    import skcuda
    SKCUDA_IMPORT = True
except ImportError:
    SKCUDA_IMPORT = False


def as_sparse_variable(x, name=None):
    """
    Convert a matrix or array into a SciPy sparse matrix if it is not already.
    
    Parameters
    ----------
    x : np.ndarray or scipy.sparse.spmatrix
        Input dense or sparse matrix.
    
    Returns
    -------
    scipy.sparse.spmatrix or np.ndarray
        Sparse matrix if input is sparse-like, else dense array.
    """
    if sparse.issparse(x):
        return x
    try:
        return sparse.csr_matrix(x)
    except Exception as e:
        raise TypeError(f"Cannot convert {type(x)} to a sparse matrix.") from e


class SolveSparse:
    """
    Solve sparse linear system C * x = b using SciPy.
    """

    def __call__(self, C, b):
        """
        Solve C * x = b.

        Parameters
        ----------
        C : scipy.sparse.spmatrix or np.ndarray
            Sparse matrix of coefficients.
        b : np.ndarray
            Right-hand side vector.

        Returns
        -------
        np.ndarray
            Solution vector x.
        """
        C = as_sparse_variable(C)
        b = np.asarray(b).flatten()  # ensure b is 1D
        weights = sparse_solve(C, b)
        return weights



import numpy as np

class NumpyGraphPro(object):
    def __init__(self, optimizer='fast_compile', verbose=None, dtype=None,
                 output=None, **kwargs):
        """
        NumPy equivalent of aesaraGraphPro.__init__
        Most aesara shared/tensor placeholders have been replaced by numpy arrays.
        Some helper methods used in the original init are stubbed below; replace them
        with your real implementations as needed.
        """
        # basic scalars
        self.lenght_of_faults = np.int32(0)
        self.pi = np.pi

        # OPTIONS
        if output is None:
            output = ['geology']
        self.output = output

        if verbose is None:
            self.verbose = [None]
        else:
            self.verbose = verbose

        self.compute_type = output

        # choose dtype
        if dtype is None:
            # original used aesara.config.device to detect cuda -> we set CPU default
            dtype = 'float64'
        self.dtype = dtype

        # Trade speed for memory this will consume more memory
        self.max_speed = kwargs.get('max_speed', 1)
        self.sparse_version = kwargs.get('sparse_version', False)

        self.gradient = kwargs.get('gradient', False)
        # no aesara device concept here
        self.device = 'cpu'

        # CONSTANT PARAMETERS FOR ALL SERIES
        # KRIGING
        self.a_T = np.ones(3, dtype=self.dtype)
        self.a_T_scalar = self.a_T
        self.c_o_T = np.ones(3, dtype=self.dtype)
        self.c_o_T_scalar = self.c_o_T

        self.n_universal_eq_T = np.zeros(5, dtype='int32')
        self.n_universal_eq_T_op = 3

        # They weight the contribution of the surface_points against the orientations.
        self.i_reescale = np.array(4., dtype=self.dtype)
        self.gi_reescale = np.array(2., dtype=self.dtype)

        # Number of dimensions. Now it is not too variable anymore
        self.n_dimensions = 3

        # This is not accumulative
        self.number_of_points_per_surface_T = np.zeros(3, dtype='int32')
        # placeholder "op" version (originally symbolic)
        self.number_of_points_per_surface_T_op = self.number_of_points_per_surface_T.copy()

        # npf: cumulative start index for each surface (like original)
        self.npf = np.cumsum(np.concatenate(([0], self.number_of_points_per_surface_T[:-1])))
        self.npf_op = self.npf
        # name metadata removed (pure numpy)

        self.nugget_effect_grad_T = np.ones(4, dtype=self.dtype)
        self.nugget_effect_scalar_T = np.ones(4, dtype=self.dtype)
        self.nugget_effect_grad_T_op = self.nugget_effect_grad_T

        # COMPUTE WEIGHTS - VARIABLES
        # placeholders for input data (original used symbolic tensors)
        # Initialize as empty arrays with sensible shapes; user will set actual data later.
        self.dips_position_all = np.empty((0, 3), dtype=self.dtype)   # Nx3
        self.dip_angles_all = np.empty((0,), dtype=self.dtype)        # N
        self.azimuth_all = np.empty((0,), dtype=self.dtype)           # N
        self.polarity_all = np.empty((0,), dtype=self.dtype)          # N

        self.surface_points_all = np.empty((0, 3), dtype=self.dtype)  # Mx3

        # len_points uses shapes from above
        self.len_points = max(0, self.surface_points_all.shape[0] - self.number_of_points_per_surface_T.shape[0])

        # Tiling dips to the 3 spatial coordinations
        self.dips_position = self.dips_position_all
        # if dips_position is empty this will be shape (3*0,3) -> still empty
        self.dips_position_tiled = np.tile(self.dips_position, (self.n_dimensions, 1))

        # These are subsets of the data for each series. Start as whole arrays; user can override.
        self.dip_angles = self.dip_angles_all
        self.azimuth = self.azimuth_all
        self.polarity = self.polarity_all

        # call stubbed helper to create rest/ref matrices (replace this stub with real implementation)
        rest_ref_aux = self.set_rest_ref_matrix(self.number_of_points_per_surface_T)
        # rest_ref_aux expected to be an iterable with at least 4 entries (the original code indexes [0..3])
        # We set default empty arrays in stub so the following assignments work.
        self.ref_layer_points_all = rest_ref_aux[0]
        self.rest_layer_points_all = rest_ref_aux[1]

        self.nugget_effect_scalar_T_ref_rest = self.set_nugget_surface_points(
            rest_ref_aux[2], rest_ref_aux[3],
            self.number_of_points_per_surface_T)

        self.nugget_effect_scalar_T_op = self.nugget_effect_scalar_T_ref_rest

        self.ref_layer_points = self.ref_layer_points_all
        self.rest_layer_points = self.rest_layer_points_all

        # fault matrix placeholder
        self.fault_matrix = np.empty((0, 0), dtype=self.dtype)

        self.input_parameters_kriging = [self.dips_position_all, self.dip_angles_all,
                                         self.azimuth_all, self.polarity_all,
                                         self.surface_points_all, self.fault_matrix]

        # COMPUTE SCALAR FIELDS - VARIABLES
        self.grid_val_T = np.empty((0, 3), dtype=self.dtype)
        self.input_parameters_export = [self.dips_position_all,
                                        self.surface_points_all,
                                        self.fault_matrix,
                                        self.grid_val_T]

        self.input_parameters_kriging_export = [self.dips_position_all,
                                                self.dip_angles_all,
                                                self.azimuth_all,
                                                self.polarity_all,
                                                self.surface_points_all,
                                                self.fault_matrix, self.grid_val_T]

        # interface_loc translation (kept as 0 as in original)
        interface_loc = 0
        # slicing with possibly empty arrays - results are empty arrays too
        self.fault_drift_at_surface_points_rest = self.fault_matrix[
            :,
            interface_loc: interface_loc + self.len_points] if self.fault_matrix.size else np.empty((0, 0), dtype=self.dtype)

        self.fault_drift_at_surface_points_ref = self.fault_matrix[
            :,
            interface_loc + self.len_points:] if self.fault_matrix.size else np.empty((0, 0), dtype=self.dtype)

        # COMPUTE BLOCKS - VARIABLES
        if self.gradient:
            self.sig_slope = np.array(50, dtype=self.dtype)
        else:
            self.sig_slope = np.array(50000, dtype=self.dtype)
            self.not_l = np.array(50., dtype=self.dtype)
            self.ellipse_factor_exponent = np.array(2., dtype=self.dtype)

        # values_properties_op placeholder
        self.values_properties_op = np.empty((0, 0), dtype=self.dtype)

        self.n_surface = np.arange(1, 5000, dtype='int32')

        self.input_parameters_block = [self.dips_position_all, self.dip_angles_all,
                                       self.azimuth_all, self.polarity_all,
                                       self.surface_points_all, self.fault_matrix,
                                       self.grid_val_T, self.values_properties_op]

        # SHARED-like arrays (now plain numpy)
        self.is_finite_ctrl = np.zeros(3, dtype='int32')
        self.inf_factor = 0
        self.is_fault = np.zeros(5000, dtype=bool)

        # COMPUTE LOOP - SHARED
        self.fault_relation = np.array([[0, 0], [0, 0]], dtype='int32')

        # Structure (original used shared aranges)
        self.n_surfaces_per_series = np.arange(2, dtype='int32')
        self.len_series_i = np.arange(2, dtype='int32')
        self.len_series_o = np.arange(2, dtype='int32')
        self.len_series_w = np.arange(2, dtype='int32')

        # Control flow placeholders (originally symbolic vectors)
        self.compute_weights_ctrl = np.zeros(0, dtype=bool)
        self.compute_scalar_ctrl = np.zeros(0, dtype=bool)
        self.compute_block_ctrl = np.zeros(0, dtype=bool)

        self.is_finite_ctrl = np.zeros(3, dtype='int32')
        self.onlap_erode_ctrl = np.zeros(3, dtype='int32')

        self.input_parameters_loop = [self.dips_position_all, self.dip_angles_all,
                                      self.azimuth_all, self.polarity_all,
                                      self.surface_points_all, self.fault_matrix,
                                      self.grid_val_T, self.values_properties_op,
                                      self.compute_weights_ctrl, self.compute_scalar_ctrl,
                                      self.compute_block_ctrl]

        self.is_erosion = np.array([1, 0], dtype='int32')
        self.is_onlap = np.array([0, 1], dtype='int32')

        self.offset = np.array(10., dtype=self.dtype)
        self.shift = 0

        # gravity branch
        if 'gravity' in self.compute_type:
            self.lg0 = np.array(0, dtype='int64')
            self.lg1 = np.array(1, dtype='int64')
            self.tz = np.empty(0, dtype=self.dtype)
            self.pos_density = np.array(1, dtype='int64')

        # magnetics branch
        if 'magnetics' in self.compute_type:
            self.lg0 = np.array(0, dtype='int64')
            self.lg1 = np.array(1, dtype='int64')
            self.B_ext = np.array(50e-6, dtype=self.dtype)  # External magnetic field in [T]
            self.incl = np.array(1., dtype=self.dtype)
            self.decl = np.array(1., dtype=self.dtype)
            self.pos_magnetics = np.array(2, dtype='int64')
            self.V = np.ones((6, 10), dtype=self.dtype)

        # topology branch
        if 'topology' in self.compute_type:
            self.max_lith = np.array(10, dtype=int)
            self.pos_topology_id = np.array(-1, dtype=int)
            self.regular_grid_res = np.ones(3, dtype=int)
            self.dxdydz = np.ones(3, dtype=self.dtype)

        # Results matrix (numpy equivalents)
        self.weights_vector = np.zeros(10000, dtype=self.dtype)
        self.scalar_fields_matrix = np.zeros((3, 10000), dtype=self.dtype)
        self.block_matrix = np.zeros((3, 3, 10000), dtype=self.dtype)
        self.mask_matrix = np.zeros((3, 10000), dtype=bool)

        # sfai: shape = (is_erosion.shape[0], n_surfaces_per_series[-1])
        # ensure n_surfaces_per_series has at least one element
        if self.n_surfaces_per_series.size > 0:
            n_last = int(self.n_surfaces_per_series[-1])
        else:
            n_last = 0
        self.sfai = np.zeros((self.is_erosion.shape[0], n_last), dtype=self.dtype)

        self.new_block = self.block_matrix.copy()
        self.new_weights = self.weights_vector.copy()
        self.new_scalar = self.scalar_fields_matrix.copy()
        self.new_mask = self.mask_matrix.copy()
        self.new_sfai = self.sfai.copy()

    def compute_weights(self):
        return self.solve_kriging()

    def compute_scalar_field(self, weights, grid, fault_matrix):
        grid_val = self.x_to_interpolate(grid)
        return self.scalar_field_at_all(weights, grid_val, fault_matrix)

    def compute_formation_block(self, Z_x, scalar_field_at_surface_points, values):
        return self.export_formation_block(Z_x, scalar_field_at_surface_points, values)

    def compute_fault_block(self, Z_x, scalar_field_at_surface_points, values,
                            n_series, grid):
        grid_val = self.x_to_interpolate(grid)
        finite_faults_ellipse = self.select_finite_faults(grid_val)
        return self.export_fault_block(Z_x, scalar_field_at_surface_points, values,
                                       finite_faults_ellipse)

    def compute_final_block(self, mask, block):
        """
        NumPy equivalent of T.sum(T.stack([mask], axis=1) * block, axis=0)
        mask: (series, points)
        block: (series, props, points)
        """
        mask_expanded = np.expand_dims(mask, axis=1)  # (series,1,points)
        final_model = np.sum(mask_expanded * block, axis=0)  # (props, points)
        return final_model

    def compute_series(self, grid=None, shift=None):
        if grid is None:
            grid = self.grid_val_T
        if shift is None:
            shift = self.shift

        # reset placeholders
        self.mask_matrix_f = np.zeros_like(self.mask_matrix)
        self.fault_matrix = np.zeros_like(self.block_matrix)

        # --- LOOPING ---
        # In aesara this was scan; here we need a python loop.
        # Placeholder: call compute_a_series for each iteration.
        series_outputs = []
        for i in range(len(self.n_surfaces_per_series) - 1):
            out = self.compute_a_series(
                self.block_matrix,
                self.weights_vector,
                self.scalar_fields_matrix,
                self.sfai,
                self.mask_matrix,
                self.mask_matrix_f,
                self.fault_matrix,
                0,  # was int64 counter
                grid,
                shift
            )
            series_outputs.append(out)

        # take the last outputs (mimic aesara.scan final state)
        if series_outputs:
            (self.block_op, self.weights_op, self.scalar_op, self.sfai_op,
             mask, fault_mask, self.fault_matrix, _) = series_outputs[-1]
        else:
            mask = np.zeros_like(self.mask_matrix)
            fault_mask = np.zeros_like(self.mask_matrix)

        mask_rev_cumprod = np.vstack([
            mask[[-1]],
            np.cumprod(~mask[:-1], axis=0)
        ])
        self.mask_op2 = mask_rev_cumprod
        block_mask = mask * mask_rev_cumprod

        fault_block = self.compute_final_block(fault_mask, self.block_op)
        final_model = self.compute_final_block(block_mask, self.block_op)

        return [final_model, self.block_op, fault_block, self.weights_op,
                self.scalar_op, self.sfai_op, block_mask, fault_mask]

    def create_oct_voxels(self, xyz, level=1):
        x_ = np.repeat(
            np.stack((xyz[:, 0] - self.dxdydz[0] / level / 4,
                      xyz[:, 0] + self.dxdydz[0] / level / 4), axis=1), 4, axis=1)
        y_ = np.tile(
            np.repeat(np.stack((xyz[:, 1] - self.dxdydz[1] / level / 4,
                                xyz[:, 1] + self.dxdydz[1] / level / 4), axis=1),
                      2, axis=1), (1, 2))
        z_ = np.tile(
            np.stack((xyz[:, 2] - self.dxdydz[2] / level / 4,
                      xyz[:, 2] + self.dxdydz[2] / level / 4), axis=1),
            (1, 4))
        return np.stack((x_.ravel(), y_.ravel(), z_.ravel()), axis=1)

    def create_oct_level_dense(self, unique_val, grid):
        uv_3d = np.round(unique_val[0, :np.prod(self.regular_grid_res)]
                         .reshape(self.regular_grid_res)).astype('int32')

        new_shape = tuple(self.regular_grid_res) + (3,)
        xyz = grid[:np.prod(self.regular_grid_res)].reshape(new_shape)

        shift_x = uv_3d[1:, :, :] - uv_3d[:-1, :, :]
        x_edg = (xyz[:-1, :, :][shift_x != 0] + xyz[1:, :, :][shift_x != 0]) / 2

        shift_y = uv_3d[:, 1:, :] - uv_3d[:, :-1, :]
        y_edg = (xyz[:, :-1, :][shift_y != 0] + xyz[:, 1:, :][shift_y != 0]) / 2

        shift_z = uv_3d[:, :, 1:] - uv_3d[:, :, :-1]
        z_edg = (xyz[:, :, :-1][shift_z != 0] + xyz[:, :, 1:][shift_z != 0]) / 2

        new_xyz_edg = np.vstack((x_edg, y_edg, z_edg))
        return self.create_oct_voxels(new_xyz_edg)

    def create_oct_level_sparse(self, unique_val, grid):
        xyz_8 = grid.reshape((-1, 8, 3))
        uv_8 = np.round(unique_val[0, :].reshape((-1, 8)))

        shift_x = uv_8[:, :4] - uv_8[:, 4:]
        x_edg = (xyz_8[:, :4, :][shift_x != 0] +
                 xyz_8[:, 4:, :][shift_x != 0]) / 2

        shift_y = uv_8[:, [0, 1, 4, 5]] - uv_8[:, [2, 3, 6, 7]]
        y_edg = (xyz_8[:, [0, 1, 4, 5], :][shift_y != 0] +
                 xyz_8[:, [2, 3, 6, 7], :][shift_y != 0]) / 2

        shift_z = uv_8[:, ::2] - uv_8[:, 1::2]
        z_edg = (xyz_8[:, ::2, :][shift_z != 0] +
                 xyz_8[:, 1::2, :][shift_z != 0]) / 2

        new_xyz_edg = np.vstack((x_edg, y_edg, z_edg))
        return self.create_oct_voxels(new_xyz_edg, level=2)

    def compute_topology(self, unique_val):
        uv_3d = np.round(unique_val[:np.prod(self.regular_grid_res)]
                         .reshape(self.regular_grid_res)).astype('int32')

        uv_l = np.hstack([
            uv_3d[1:, :, :].reshape((1, -1)),
            uv_3d[:, 1:, :].reshape((1, -1)),
            uv_3d[:, :, 1:].reshape((1, -1))
        ])

        uv_r = np.hstack([
            uv_3d[:-1, :, :].reshape((1, -1)),
            uv_3d[:, :-1, :].reshape((1, -1)),
            uv_3d[:, :, :-1].reshape((1, -1))
        ])

        shift = uv_l - uv_r
        select_edges = (shift.reshape((1, -1)) != 0)
        select_edges_dir = select_edges.reshape((3, -1))

        select_voxels = np.zeros_like(uv_3d, dtype=int)

        # Update voxels where edges differ
        select_voxels[1:, :, :] += select_edges_dir[0].reshape(self.regular_grid_res - np.array([1, 0, 0]))
        select_voxels[:-1, :, :] += select_edges_dir[0].reshape(self.regular_grid_res - np.array([1, 0, 0]))

        select_voxels[:, 1:, :] += select_edges_dir[1].reshape(self.regular_grid_res - np.array([0, 1, 0]))
        select_voxels[:, :-1, :] += select_edges_dir[1].reshape(self.regular_grid_res - np.array([0, 1, 0]))

        select_voxels[:, :, 1:] += select_edges_dir[2].reshape(self.regular_grid_res - np.array([0, 0, 1]))
        select_voxels[:, :, :-1] += select_edges_dir[2].reshape(self.regular_grid_res - np.array([0, 0, 1]))

        uv_lr = np.vstack([uv_l.reshape((1, -1)), uv_r.reshape((1, -1))])
        uv_lr_boundaries = uv_lr[np.tile(select_edges.reshape((1, -1)), (2, 1))].reshape((2, -1)).T

        # unique edges and counts
        edges_id, count_edges = np.unique(uv_lr_boundaries, axis=0, return_counts=True)
        return select_voxels, edges_id, count_edges

    def get_boundary_voxels(self, unique_val):
        uv_3d = np.round(unique_val[0, :np.prod(self.regular_grid_res)]
                         .reshape(self.regular_grid_res)).astype('int32')

        uv_l = np.hstack([
            uv_3d[1:, :, :].reshape((1, -1)),
            uv_3d[:, 1:, :].reshape((1, -1)),
            uv_3d[:, :, 1:].reshape((1, -1))
        ])

        uv_r = np.hstack([
            uv_3d[:-1, :, :].reshape((1, -1)),
            uv_3d[:, :-1, :].reshape((1, -1)),
            uv_3d[:, :, :-1].reshape((1, -1))
        ])

        shift = uv_l - uv_r
        select_edges = (shift.reshape((1, -1)) != 0)
        select_edges_dir = select_edges.reshape((3, -1))

        select_voxels = np.zeros_like(uv_3d, dtype=int)

        select_voxels[1:, :, :] += select_edges_dir[0].reshape(self.regular_grid_res - np.array([1, 0, 0]))
        select_voxels[:-1, :, :] += select_edges_dir[0].reshape(self.regular_grid_res - np.array([1, 0, 0]))

        select_voxels[:, 1:, :] += select_edges_dir[1].reshape(self.regular_grid_res - np.array([0, 1, 0]))
        select_voxels[:, :-1, :] += select_edges_dir[1].reshape(self.regular_grid_res - np.array([0, 1, 0]))

        select_voxels[:, :, 1:] += select_edges_dir[2].reshape(self.regular_grid_res - np.array([0, 0, 1]))
        select_voxels[:, :, :-1] += select_edges_dir[2].reshape(self.regular_grid_res - np.array([0, 0, 1]))

        return select_voxels

    # region Geometry
    def repeat_list(self, val, r_0, r_1, repeated_array, n_col):
        """
        Repeat val into repeated_array[r_0:r_1, :]
        """
        repeated_array[r_0:r_1] = np.full((r_1 - r_0, n_col), val)
        return repeated_array

    def set_rest_ref_matrix(self, number_of_points_per_surface):
        ref_positions = np.cumsum(np.concatenate(([0], number_of_points_per_surface[:-1] + 1)))
        cum_rep = np.cumsum(np.concatenate(([0], number_of_points_per_surface)))

        ref_points = np.zeros((cum_rep[-1], 3))
        for i in range(len(ref_positions)):
            r0, r1 = cum_rep[i], cum_rep[i+1]
            ref_points = self.repeat_list(self.surface_points_all[ref_positions[i]], r0, r1, ref_points, 3)

        rest_mask = np.ones(self.surface_points_all.shape[0], dtype='int16')
        rest_mask[ref_positions] = 0
        rest_mask = np.nonzero(rest_mask)[0]
        rest_points = self.surface_points_all[rest_mask]
        return [ref_points, rest_points, ref_positions, rest_mask]

    def set_nugget_surface_points(self, ref_positions, rest_mask, number_of_points_per_surface):
        cum_rep = np.cumsum(np.concatenate(([0], number_of_points_per_surface)))
        ref_nugget = np.zeros((cum_rep[-1], 1))
        for i in range(len(ref_positions)):
            r0, r1 = cum_rep[i], cum_rep[i+1]
            ref_nugget = self.repeat_list(self.nugget_effect_scalar_T[ref_positions[i]], r0, r1, ref_nugget, 1)

        rest_nugget = self.nugget_effect_scalar_T[rest_mask]
        nugget_rest_ref = ref_nugget.reshape((1, -1))[0] + rest_nugget
        return nugget_rest_ref
    
    @staticmethod
    def squared_euclidean_distances(x_1, x_2):
        """
        Compute the Euclidean distances in 3D between all the points in x_1 and x_2.
        """
        sqd = np.sqrt(np.maximum(
            np.sum(x_1 ** 2, axis=1).reshape(-1, 1) +
            np.sum(x_2 ** 2, axis=1).reshape(1, -1) -
            2 * np.dot(x_1, x_2.T), 1e-12
        ))
        return sqd

    def matrices_shapes(self):
        """
        Get all the lengths of the matrices that form the covariance matrix.
        """
        length_of_CG = self.dips_position_tiled.shape[0]
        length_of_CGI = self.rest_layer_points.shape[0]
        length_of_U_I = self.n_universal_eq_T_op
        length_of_faults = self.lenght_of_faults
        length_of_C = length_of_CG + length_of_CGI + length_of_U_I + length_of_faults

        return length_of_CG, length_of_CGI, length_of_U_I, length_of_faults, length_of_C

    # region Kriging
    def cov_surface_points(self):
        """Create covariance function for the surface_points"""
        sed_rest_rest = self.squared_euclidean_distances(self.rest_layer_points,
                                                         self.rest_layer_points)
        sed_ref_rest = self.squared_euclidean_distances(self.ref_layer_points,
                                                        self.rest_layer_points)
        sed_rest_ref = self.squared_euclidean_distances(self.rest_layer_points,
                                                        self.ref_layer_points)
        sed_ref_ref = self.squared_euclidean_distances(self.ref_layer_points,
                                                       self.ref_layer_points)

        def cov_kernel(d):
            return ((d < self.a_T_scalar) *
                    (1 - 7 * (d / self.a_T_scalar) ** 2 +
                     35 / 4 * (d / self.a_T_scalar) ** 3 -
                     7 / 2 * (d / self.a_T_scalar) ** 5 +
                     3 / 4 * (d / self.a_T_scalar) ** 7))

        C_I = (self.c_o_T_scalar * self.i_reescale *
               (cov_kernel(sed_rest_rest) -
                cov_kernel(sed_ref_rest) -
                cov_kernel(sed_rest_ref) +
                cov_kernel(sed_ref_ref)))

        C_I += np.eye(C_I.shape[0]) * self.nugget_effect_scalar_T_op
        return C_I

    def cov_gradients(self, verbose=0):
        """Covariance function for the gradients"""
        sed_dips_dips = self.squared_euclidean_distances(self.dips_position_tiled,
                                                         self.dips_position_tiled)

        # Cartesian differences
        h_u = np.vstack([
            np.tile(self.dips_position[:, 0] - self.dips_position[:, 0].reshape(-1, 1), (self.n_dimensions, 1)),
            np.tile(self.dips_position[:, 1] - self.dips_position[:, 1].reshape(-1, 1), (self.n_dimensions, 1)),
            np.tile(self.dips_position[:, 2] - self.dips_position[:, 2].reshape(-1, 1), (self.n_dimensions, 1))
        ])
        h_v = h_u.T

        # Block diagonal mask
        perpendicularity_matrix = np.zeros_like(sed_dips_dips)
        n = self.dips_position.shape[0]
        perpendicularity_matrix[:n, :n] = 1
        perpendicularity_matrix[n:2*n, n:2*n] = 1
        perpendicularity_matrix[2*n:3*n, 2*n:3*n] = 1

        with np.errstate(divide='ignore', invalid='ignore'):
            C_G = np.where(
                sed_dips_dips == 0,
                0,
                ((h_u * h_v / sed_dips_dips ** 2) *
                 ((sed_dips_dips < self.a_T_scalar) *
                  (-self.c_o_T_scalar * (
                          -14 / self.a_T_scalar ** 2 +
                          105 / 4 * sed_dips_dips / self.a_T_scalar ** 3 -
                          35 / 2 * sed_dips_dips ** 3 / self.a_T_scalar ** 5 +
                          21 / 4 * sed_dips_dips ** 5 / self.a_T_scalar ** 7))
                  + (sed_dips_dips < self.a_T_scalar) *
                  self.c_o_T_scalar * 7 *
                  (9 * sed_dips_dips ** 5 -
                   20 * self.a_T_scalar ** 2 * sed_dips_dips ** 3 +
                   15 * self.a_T_scalar ** 4 * sed_dips_dips -
                   4 * self.a_T_scalar ** 5) /
                  (2 * self.a_T_scalar ** 7))
                 - (perpendicularity_matrix *
                    (sed_dips_dips < self.a_T_scalar) *
                    self.c_o_T_scalar *
                    (-14 / self.a_T_scalar ** 2 +
                     105 / 4 * sed_dips_dips / self.a_T_scalar ** 3 -
                     35 / 2 * sed_dips_dips ** 3 / self.a_T_scalar ** 5 +
                     21 / 4 * sed_dips_dips ** 5 / self.a_T_scalar ** 7)))
            )

        C_G += np.eye(C_G.shape[0]) * self.nugget_effect_grad_T_op
        return C_G

    def cov_interface_gradients(self):
        """Cross covariance gradients–surface points"""
        sed_dips_rest = self.squared_euclidean_distances(self.dips_position_tiled,
                                                         self.rest_layer_points)
        sed_dips_ref = self.squared_euclidean_distances(self.dips_position_tiled,
                                                        self.ref_layer_points)

        hu_rest = np.vstack([
            (self.dips_position[:, 0] - self.rest_layer_points[:, 0].reshape(-1, 1)).T,
            (self.dips_position[:, 1] - self.rest_layer_points[:, 1].reshape(-1, 1)).T,
            (self.dips_position[:, 2] - self.rest_layer_points[:, 2].reshape(-1, 1)).T
        ])

        hu_ref = np.vstack([
            (self.dips_position[:, 0] - self.ref_layer_points[:, 0].reshape(-1, 1)).T,
            (self.dips_position[:, 1] - self.ref_layer_points[:, 1].reshape(-1, 1)).T,
            (self.dips_position[:, 2] - self.ref_layer_points[:, 2].reshape(-1, 1)).T
        ])

        C_GI = self.gi_reescale * (
            hu_rest * (sed_dips_rest < self.a_T_scalar) *
            (-self.c_o_T_scalar * (
                    -14 / self.a_T_scalar ** 2 +
                    105 / 4 * sed_dips_rest / self.a_T_scalar ** 3 -
                    35 / 2 * sed_dips_rest ** 3 / self.a_T_scalar ** 5 +
                    21 / 4 * sed_dips_rest ** 5 / self.a_T_scalar ** 7))
            - hu_ref * (sed_dips_ref < self.a_T_scalar) *
            (-self.c_o_T_scalar * (
                    -14 / self.a_T_scalar ** 2 +
                    105 / 4 * sed_dips_ref / self.a_T_scalar ** 3 -
                    35 / 2 * sed_dips_ref ** 3 / self.a_T_scalar ** 5 +
                    21 / 4 * sed_dips_ref ** 5 / self.a_T_scalar ** 7))
        ).T
        return C_GI

    def universal_matrix(self):
        """Drift matrices for potential field and gradients"""
        n = self.dips_position.shape[0]
        U_G = np.zeros((n * self.n_dimensions, 9))

        # Linear terms
        U_G[:n, 0] = 1
        U_G[n:2*n, 1] = 1
        U_G[2*n:3*n, 2] = 1

        # Quadratic terms
        U_G[:n, 3] = 2 * self.gi_reescale * self.dips_position[:, 0]
        U_G[n:2*n, 4] = 2 * self.gi_reescale * self.dips_position[:, 1]
        U_G[2*n:3*n, 5] = 2 * self.gi_reescale * self.dips_position[:, 2]

        U_G[:n, 6] = self.gi_reescale * self.dips_position[:, 1]
        U_G[n:2*n, 6] = self.gi_reescale * self.dips_position[:, 0]

        U_G[:n, 7] = self.gi_reescale * self.dips_position[:, 2]
        U_G[2*n:3*n, 7] = self.gi_reescale * self.dips_position[:, 0]

        U_G[n:2*n, 8] = self.gi_reescale * self.dips_position[:, 2]
        U_G[2*n:3*n, 8] = self.gi_reescale * self.dips_position[:, 1]

        U_I = -np.column_stack([
            self.gi_reescale * (self.rest_layer_points[:, 0] - self.ref_layer_points[:, 0]),
            self.gi_reescale * (self.rest_layer_points[:, 1] - self.ref_layer_points[:, 1]),
            self.gi_reescale * (self.rest_layer_points[:, 2] - self.ref_layer_points[:, 2]),
            self.gi_reescale**2 * (self.rest_layer_points[:, 0]**2 - self.ref_layer_points[:, 0]**2),
            self.gi_reescale**2 * (self.rest_layer_points[:, 1]**2 - self.ref_layer_points[:, 1]**2),
            self.gi_reescale**2 * (self.rest_layer_points[:, 2]**2 - self.ref_layer_points[:, 2]**2),
            self.gi_reescale**2 * (self.rest_layer_points[:, 0]*self.rest_layer_points[:, 1] -
                                   self.ref_layer_points[:, 0]*self.ref_layer_points[:, 1]),
            self.gi_reescale**2 * (self.rest_layer_points[:, 0]*self.rest_layer_points[:, 2] -
                                   self.ref_layer_points[:, 0]*self.ref_layer_points[:, 2]),
            self.gi_reescale**2 * (self.rest_layer_points[:, 1]*self.rest_layer_points[:, 2] -
                                   self.ref_layer_points[:, 1]*self.ref_layer_points[:, 2]),
        ])
        return U_I[:, :self.n_universal_eq_T_op], U_G[:, :self.n_universal_eq_T_op]

    def faults_matrix(self, f_ref=None, f_res=None):
        """Drift matrices for faults"""
        length_of_CG, length_of_CGI, length_of_U_I, length_of_faults = self.matrices_shapes()[:4]

        F_I = (self.fault_drift_at_surface_points_ref -
               self.fault_drift_at_surface_points_rest) + 1e-4
        F_G = np.zeros((length_of_faults, length_of_CG)) + 1e-4

        return F_I, F_G
    
    def covariance_matrix(self):
        """
        Build the universal cokriging covariance matrix (NumPy version).
        """
        # Lengths
        length_of_CG, length_of_CGI, length_of_U_I, length_of_faults, length_of_C = self.matrices_shapes()

        # Individual covariance blocks
        C_G = self.cov_gradients()
        C_I = self.cov_surface_points()
        C_GI = self.cov_interface_gradients()
        U_I, U_G = self.universal_matrix()
        F_I, F_G = self.faults_matrix()

        # Allocate covariance matrix
        C_matrix = np.zeros((length_of_C, length_of_C), dtype=self.dtype)

        # First block row
        C_matrix[0:length_of_CG, 0:length_of_CG] = C_G
        C_matrix[0:length_of_CG, length_of_CG:length_of_CG + length_of_CGI] = C_GI.T
        C_matrix[0:length_of_CG,
                 length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I] = U_G
        C_matrix[0:length_of_CG,
                 length_of_CG + length_of_CGI + length_of_U_I:] = F_G.T

        # Second block row
        C_matrix[length_of_CG:length_of_CG + length_of_CGI, 0:length_of_CG] = C_GI
        C_matrix[length_of_CG:length_of_CG + length_of_CGI,
                 length_of_CG:length_of_CG + length_of_CGI] = C_I
        C_matrix[length_of_CG:length_of_CG + length_of_CGI,
                 length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I] = U_I
        C_matrix[length_of_CG:length_of_CG + length_of_CGI,
                 length_of_CG + length_of_CGI + length_of_U_I:] = F_I.T

        # Third block row (U_G, U_I)
        C_matrix[length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I,
                 0:length_of_CG] = U_G.T
        C_matrix[length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I,
                 length_of_CG:length_of_CG + length_of_CGI] = U_I.T

        # Fourth block row (F_G, F_I)
        C_matrix[length_of_CG + length_of_CGI + length_of_U_I:,
                 0:length_of_CG] = F_G
        C_matrix[length_of_CG + length_of_CGI + length_of_U_I:,
                 length_of_CG:length_of_CG + length_of_CGI] = F_I

        return C_matrix

    def b_vector(self, dip_angles=None, azimuth=None, polarity=None):
        """
        Build the independent vector b for the kriging system (NumPy version).
        """
        length_of_C = self.matrices_shapes()[-1]

        if dip_angles is None:
            dip_angles = self.dip_angles
        if azimuth is None:
            azimuth = self.azimuth
        if polarity is None:
            polarity = self.polarity

        # Cartesian components of dips
        G_x = np.sin(np.deg2rad(dip_angles)) * np.sin(np.deg2rad(azimuth)) * polarity
        G_y = np.sin(np.deg2rad(dip_angles)) * np.cos(np.deg2rad(azimuth)) * polarity
        G_z = np.cos(np.deg2rad(dip_angles)) * polarity

        G = np.concatenate([G_x, G_y, G_z])

        b = np.zeros((length_of_C,), dtype=self.dtype)
        b[0:G.shape[0]] = G
        return b

    def solve_kriging(self, b=None):
        """
        Solve the kriging system using NumPy / SciPy solvers.
        """
        C_matrix = self.covariance_matrix()
        if b is None:
            b = self.b_vector()

        if self.sparse_version:
            C_sparse = csr_matrix(C_matrix)
            b_sparse = b
            DK_parameters = spsolve(C_sparse, b_sparse)
        else:
            DK_parameters = np.linalg.solve(C_matrix, b)

        DK_parameters = DK_parameters.reshape(-1)
        return DK_parameters
    
        # region Evaluate
    def x_to_interpolate(self, grid, verbose=0):
        """
        Add rest and reference points to the interpolation grid.
        Returns:
            np.ndarray: The 3D points of the given grid plus the reference and rest points
        """
        grid_val = np.concatenate(
            [grid, self.rest_layer_points_all, self.ref_layer_points_all],
            axis=0
        )
        return grid_val

    def extend_dual_kriging(self, weights, grid_shape):
        """
        Tile the dual kriging vector to cover all the points to interpolate.
        Returns:
            np.ndarray: Matrix with the DK parameters repeated for all grid points
        """
        DK_parameters = weights
        DK_weights = np.tile(DK_parameters.reshape(-1, 1), (1, grid_shape))
        return DK_weights
    # endregion

    # region Evaluate Geology
    def contribution_gradient_interface(self, grid_val=None, weights=None):
        """
        Contribution of foliations (gradients) at every grid point.
        """
        if weights is None:
            weights = self.extend_dual_kriging()
        if grid_val is None:
            grid_val = self.x_to_interpolate()

        length_of_CG = self.matrices_shapes()[0]

        # Cartesian differences between dips and grid points
        hu_SimPoint = np.vstack([
            (self.dips_position[:, 0][:, None] - grid_val[:, 0]).T,
            (self.dips_position[:, 1][:, None] - grid_val[:, 1]).T,
            (self.dips_position[:, 2][:, None] - grid_val[:, 2]).T,
        ])

        # Euclidean distances
        sed_dips_SimPoint = self.squared_euclidean_distances(
            self.dips_position_tiled, grid_val
        )

        if self.sparse_version:
            cov_aux = csr_matrix(
                self.gi_reescale *
                (-hu_SimPoint *
                 (sed_dips_SimPoint < self.a_T_scalar) *
                 (- self.c_o_T_scalar *
                  ((-14 / self.a_T_scalar ** 2) +
                   105 / 4 * sed_dips_SimPoint / self.a_T_scalar ** 3 -
                   35 / 2 * sed_dips_SimPoint ** 3 / self.a_T_scalar ** 5 +
                   21 / 4 * sed_dips_SimPoint ** 5 / self.a_T_scalar ** 7)))
            )
            sliced_weights = weights[:length_of_CG]
            sigma_0_grad = sliced_weights @ cov_aux.toarray()
        else:
            sigma_0_grad = np.sum(
                (weights[:length_of_CG][:, None] *
                 self.gi_reescale *
                 (-hu_SimPoint *
                  (sed_dips_SimPoint < self.a_T_scalar) *
                  (- self.c_o_T_scalar *
                   ((-14 / self.a_T_scalar ** 2) +
                    105 / 4 * sed_dips_SimPoint / self.a_T_scalar ** 3 -
                    35 / 2 * sed_dips_SimPoint ** 3 / self.a_T_scalar ** 5 +
                    21 / 4 * sed_dips_SimPoint ** 5 / self.a_T_scalar ** 7)))),
                axis=0
            )
        return sigma_0_grad

    def contribution_interface(self, grid_val, weights=None):
        """
        Contribution of surface points at every grid point.
        """
        if weights is None:
            weights = self.compute_weights()

        length_of_CG, length_of_CGI = self.matrices_shapes()[:2]

        # Euclidean distances
        sed_rest_SimPoint = self.squared_euclidean_distances(
            self.rest_layer_points, grid_val
        )
        sed_ref_SimPoint = self.squared_euclidean_distances(
            self.ref_layer_points, grid_val
        )

        if self.sparse_version:
            cov_aux = csr_matrix(
                self.c_o_T_scalar * self.i_reescale * (
                    ((sed_rest_SimPoint < self.a_T_scalar) *
                     (1 - 7 * (sed_rest_SimPoint / self.a_T_scalar) ** 2 +
                      35 / 4 * (sed_rest_SimPoint / self.a_T_scalar) ** 3 -
                      7 / 2 * (sed_rest_SimPoint / self.a_T_scalar) ** 5 +
                      3 / 4 * (sed_rest_SimPoint / self.a_T_scalar) ** 7)) -
                    ((sed_ref_SimPoint < self.a_T_scalar) *
                     (1 - 7 * (sed_ref_SimPoint / self.a_T_scalar) ** 2 +
                      35 / 4 * (sed_ref_SimPoint / self.a_T_scalar) ** 3 -
                      7 / 2 * (sed_ref_SimPoint / self.a_T_scalar) ** 5 +
                      3 / 4 * (sed_ref_SimPoint / self.a_T_scalar) ** 7)))
            )
            weights_sliced = -weights[length_of_CG:length_of_CG + length_of_CGI]
            sigma_0_interf = weights_sliced @ cov_aux.toarray()
        else:
            sigma_0_interf = np.sum(
                -weights[length_of_CG:length_of_CG + length_of_CGI][:, None] *
                (self.c_o_T_scalar * self.i_reescale *
                 ((sed_rest_SimPoint < self.a_T_scalar) *
                  (1 - 7 * (sed_rest_SimPoint / self.a_T_scalar) ** 2 +
                   35 / 4 * (sed_rest_SimPoint / self.a_T_scalar) ** 3 -
                   7 / 2 * (sed_rest_SimPoint / self.a_T_scalar) ** 5 +
                   3 / 4 * (sed_rest_SimPoint / self.a_T_scalar) ** 7) -
                  ((sed_ref_SimPoint < self.a_T_scalar) *
                   (1 - 7 * (sed_ref_SimPoint / self.a_T_scalar) ** 2 +
                    35 / 4 * (sed_ref_SimPoint / self.a_T_scalar) ** 3 -
                    7 / 2 * (sed_ref_SimPoint / self.a_T_scalar) ** 5 +
                    3 / 4 * (sed_ref_SimPoint / self.a_T_scalar) ** 7)))),
                axis=0
            )
        return sigma_0_interf
    
    def contribution_universal_drift(self, grid_val, weights=None):
        """
        Contribution of the universal drift at every interpolation point.
        """
        if weights is None:
            weights = self.compute_weights()

        length_of_CG, length_of_CGI, length_of_U_I, length_of_faults, length_of_C = self.matrices_shapes()

        # Build universal drift matrix
        universal_grid_surface_points_matrix = np.hstack([
            grid_val,
            grid_val ** 2,
            np.stack([
                grid_val[:, 0] * grid_val[:, 1],
                grid_val[:, 0] * grid_val[:, 2],
                grid_val[:, 1] * grid_val[:, 2]
            ], axis=1)
        ]).T

        i_rescale_aux = np.tile(self.gi_reescale, 9)
        i_rescale_aux[:3] = 1
        _aux_magic_term = np.tile(
            i_rescale_aux[:self.n_universal_eq_T_op],
            (grid_val.shape[0], 1)
        ).T

        if self.sparse_version:
            f_0 = weights[length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I] @ (
                self.gi_reescale *
                _aux_magic_term *
                universal_grid_surface_points_matrix[:self.n_universal_eq_T_op]
            )
        else:
            universal_kernel = (
                self.gi_reescale *
                _aux_magic_term *
                universal_grid_surface_points_matrix[:self.n_universal_eq_T_op]
            )
            f_0 = np.sum(
                weights[length_of_CG + length_of_CGI:length_of_CG + length_of_CGI + length_of_U_I][:, None] *
                universal_kernel,
                axis=0
            )
        return f_0

    def contribution_faults(self, weights=None, a=0, b=100000000, f_m=None):
        """
        Contribution of faults to the potential field.
        """
        if weights is None:
            weights = self.compute_weights()

        length_of_CG, length_of_CGI, length_of_U_I, length_of_faults, length_of_C = self.matrices_shapes()

        fault_matrix_selection_non_zero = f_m[:, a:b]

        if self.sparse_version:
            f_1 = weights[length_of_CG + length_of_CGI + length_of_U_I:] @ fault_matrix_selection_non_zero
        else:
            f_1 = np.sum(
                weights[length_of_CG + length_of_CGI + length_of_U_I:, None] *
                fault_matrix_selection_non_zero,
                axis=0
            )
        return f_1

    def scalar_field_loop(self, a, b, Z_x, grid_val, weights, fault_matrix):
        """
        Loop computing partial scalar fields over slices of the grid.
        """
        if self.sparse_version:
            rang = 5
            tiled_weights = self.extend_dual_kriging(weights, rang)
            sigma_0_grad = self.contribution_gradient_interface(grid_val[a:b], tiled_weights)
            sigma_0_interf = self.contribution_interface(grid_val[a:b], tiled_weights)
            f_0 = self.contribution_universal_drift(grid_val[a:b], tiled_weights)
            f_1 = self.contribution_faults(tiled_weights, a, b)

            partial_Z_x = (sigma_0_grad + sigma_0_interf + f_0 + f_1)[0]
        else:
            rang = b - a
            tiled_weights = self.extend_dual_kriging(weights, rang)
            sigma_0_grad = self.contribution_gradient_interface(grid_val[a:b], tiled_weights)
            sigma_0_interf = self.contribution_interface(grid_val[a:b], tiled_weights)
            f_0 = self.contribution_universal_drift(grid_val[a:b], tiled_weights)
            f_1 = self.contribution_faults(tiled_weights, a, b, fault_matrix)

            partial_Z_x = sigma_0_grad + sigma_0_interf + f_0 + f_1

        Z_x[a:b] = partial_Z_x
        return Z_x

    def scalar_field_at_all(self, weights, grid_val, fault_matrix):
        """
        Compute the potential field at all interpolation points (grid + rest + ref).
        """
        grid_shape = grid_val.shape[0]
        Z_x_init = np.zeros(grid_shape)

        # Adaptive chunk size (if memory issues)
        steps = int(5e6 / self.matrices_shapes()[-1])
        slices = np.concatenate([np.arange(0, grid_shape, steps, dtype=np.int64), [grid_shape]])

        if self.sparse_version:
            self.dot_version = True
            Z_x = self.scalar_field_loop(0, grid_shape, Z_x_init, grid_val, weights, fault_matrix)

        elif self.max_speed < 2:
            Z_x = Z_x_init.copy()
            for i in range(len(slices) - 1):
                a, b = slices[i], slices[i+1]
                Z_x = self.scalar_field_loop(a, b, Z_x, grid_val, weights, fault_matrix)
        else:
            tiled_weights = self.extend_dual_kriging(weights, grid_val.shape[0])
            Z_x = self.scalar_field_loop(0, grid_shape, Z_x_init, grid_val, tiled_weights, fault_matrix)

        return Z_x
    
    import numpy as np


    def get_scalar_field_at_surface_points(self, Z_x, npf_op=None):
        if npf_op is None:
            npf_op = self.npf_op
        # Extract scalar field values at surface points
        scalar_field_at_surface_points_values = Z_x[-2 * self.len_points: -self.len_points][npf_op]
        return scalar_field_at_surface_points_values

    def select_finite_faults(self, grid):
        fault_points = np.vstack([self.ref_layer_points[0][None, :], self.rest_layer_points]).T
        ctr = np.mean(fault_points, axis=1)
        x = fault_points - ctr.reshape((-1, 1))
        M = x @ x.T
        U, D, Vt = np.linalg.svd(M)
        V = Vt.T
        rotated_x = grid @ U @ V
        rotated_fault_points = (fault_points.T @ U) @ V
        rotated_ctr = np.mean(rotated_fault_points, axis=0)

        a_radius = (rotated_fault_points[:, 0].max() - rotated_fault_points[:, 0].min()) / 2
        b_radius = (rotated_fault_points[:, 1].max() - rotated_fault_points[:, 1].min()) / 2

        ellipse_factor = ((rotated_x[:, 0] - rotated_ctr[0]) ** 2 / a_radius ** 2 +
                          (rotated_x[:, 1] - rotated_ctr[1]) ** 2 / b_radius ** 2)
        return ellipse_factor

    def compare(self, a, b, slice_init, Z_x, l, n_surface, drift):
        """
        Sigmoid-based thresholding for lithology segmentation.
        """
        n_surface_0 = n_surface[:, slice_init:slice_init + 1]
        n_surface_1 = n_surface[:, slice_init + 1:slice_init + 2]
        drift_slice = drift[:, slice_init:slice_init + 1]

        sigm = (-n_surface_0.reshape((-1, 1)) / (1 + np.exp(-l * (Z_x - a)))
                - n_surface_1.reshape((-1, 1)) / (1 + np.exp(l * (Z_x - b)))
                + drift_slice.reshape((-1, 1)))
        return sigm

    def export_fault_block(self, Z_x, scalar_field_at_surface_points,
                           values_properties_op, finite_faults_sel,
                           slope=50, offset_slope=950):
        """
        Export block model values for faults.
        """
        max_pot = np.max(Z_x)
        min_pot = np.min(Z_x)
        boundary_pad = (max_pot - min_pot) * 0.01

        ellipse_factor_rectified = np.where(finite_faults_sel < 1., finite_faults_sel, 1.)
        sigmoid_slope = (offset_slope
                         - offset_slope * ellipse_factor_rectified ** self.ellipse_factor_exponent
                         + self.not_l)

        scalar_field_iter = np.concatenate([
            np.array([max_pot + boundary_pad]),
            scalar_field_at_surface_points,
            np.array([min_pot - boundary_pad])
        ])

        n_surface_op_float_sigmoid = np.repeat(values_properties_op[[0], :], 2, axis=1)
        n_surface_op_float_sigmoid[:, 1] = -1
        n_surface_op_float_sigmoid[:, -1] = -1
        drift = n_surface_op_float_sigmoid.copy()
        drift[:, 0] = n_surface_op_float_sigmoid[:, 1]

        fault_block = []
        for i in range(0, n_surface_op_float_sigmoid.shape[1], 2):
            a = scalar_field_iter[i]
            b = scalar_field_iter[i + 1]
            sigm = self.compare(a, b, i, Z_x, sigmoid_slope, n_surface_op_float_sigmoid, drift)
            fault_block.append(sigm)

        fault_block = np.sum(np.vstack(fault_block), axis=0)
        return fault_block

    def export_formation_block(self, Z_x, scalar_field_at_surface_points,
                               values_properties_op):
        """
        Export block model values for formations.
        """
        slope = self.sig_slope
        if self.max_speed < 1:
            max_pot = np.max(Z_x)
            min_pot = np.min(Z_x)
            l = slope / (max_pot - min_pot)
        else:
            l = slope

        scalar_field_iter = np.concatenate([
            np.array([0]),
            scalar_field_at_surface_points,
            np.array([0])
        ])

        n_surface_op_float_sigmoid = np.repeat(values_properties_op, 2, axis=1)
        n_surface_op_float_sigmoid[:, 0] = 0
        n_surface_op_float_sigmoid[:, -1] = 0
        drift = n_surface_op_float_sigmoid.copy()
        drift[:, 0] = n_surface_op_float_sigmoid[:, 1]

        formations_block = []
        for i in range(0, n_surface_op_float_sigmoid.shape[1], 2):
            a = scalar_field_iter[i]
            b = scalar_field_iter[i + 1]
            sigm = self.compare(a, b, i, Z_x, l, n_surface_op_float_sigmoid, drift)
            formations_block.append(sigm)

        formations_block = np.sum(np.vstack(formations_block), axis=0)

        if self.gradient:
            ReLU_up = np.where(Z_x < scalar_field_iter[1], 0, -0.01 * (Z_x - scalar_field_iter[1]))
            ReLU_down = np.where(Z_x > scalar_field_iter[-2], 0, 0.01 * np.abs(Z_x - scalar_field_iter[-2]))
            formations_block += ReLU_down + ReLU_up

        return formations_block

    import numpy as np

    def compute_a_series(self,
                            len_i_0=0, len_i_1=None,
                            len_f_0=0, len_f_1=None,
                            len_w_0=0, len_w_1=None,
                            n_form_per_serie_0=0, n_form_per_serie_1=None,
                            u_grade_iter=3,
                            compute_weight_ctr=True,
                            compute_scalar_ctr=True,
                            compute_block_ctr=True,
                            is_finite=False, is_erosion=True,
                            is_onlap=False,
                            n_series=0,
                            range_val=10., c_o=10.,
                            block_matrix=None, weights_vector=None,
                            scalar_field_matrix=None, sfai=None, mask_matrix=None,
                            mask_matrix_f=None, fault_matrix=None, nsle=0, grid=None,
                            shift=None
                            ):
        self.a_T_scalar = range_val
        self.c_o_T_scalar = c_o

        self.number_of_points_per_surface_T_op = self.number_of_points_per_surface_T[
            n_form_per_serie_0:n_form_per_serie_1
        ]

        self.npf_op = self.npf[n_form_per_serie_0:n_form_per_serie_1]
        n_surface_op = self.n_surface[n_form_per_serie_0:n_form_per_serie_1]

        self.dips_position = self.dips_position_all[len_f_0:len_f_1, :]
        self.dips_position_tiled = np.tile(self.dips_position, (self.n_dimensions, 1))

        self.dip_angles = self.dip_angles_all[len_f_0:len_f_1]
        self.azimuth = self.azimuth_all[len_f_0:len_f_1]
        self.polarity = self.polarity_all[len_f_0:len_f_1]

        self.ref_layer_points = self.ref_layer_points_all[len_i_0:len_i_1, :]
        self.rest_layer_points = self.rest_layer_points_all[len_i_0:len_i_1, :]

        self.nugget_effect_scalar_T_op = self.nugget_effect_scalar_T_ref_rest[len_i_0:len_i_1]
        self.nugget_effect_grad_T_op = self.nugget_effect_grad_T[len_f_0*3:len_f_1*3]

        self.n_universal_eq_T_op = u_grade_iter

        x_to_interpolate_shape = grid.shape[0] + 2 * self.len_points

        # Extract faults matrices
        faults_relation_op = self.fault_relation[:, int(n_series)]
        fault_idxs = np.nonzero(faults_relation_op.astype(np.int8))[0]
        fault_matrix_op = fault_matrix[fault_idxs, 0, shift:x_to_interpolate_shape + shift] * self.offset

        self.lenght_of_faults = fault_matrix_op.shape[0]

        interface_loc = grid.shape[0]

        self.fault_drift_at_surface_points_rest = fault_matrix_op[:, interface_loc + len_i_0: interface_loc + len_i_1]
        self.fault_drift_at_surface_points_ref = fault_matrix_op[:, interface_loc + self.len_points + len_i_0:
                                                                interface_loc + self.len_points + len_i_1]

        b = self.b_vector(self.dip_angles, self.azimuth, self.polarity)

        # Compute weights
        if self.sparse_version:
            weights = self.solve_kriging(b)
            Z_x = self.compute_scalar_field(weights, grid)
            weights = weights[0]
        else:
            if compute_weight_ctr:
                weights = self.solve_kriging(b)
            else:
                weights = weights_vector[len_w_0:len_w_1]

            if compute_scalar_ctr:
                Z_x = self.compute_scalar_field(weights, grid, fault_matrix_op)
            else:
                Z_x = scalar_field_matrix[n_series]

        scalar_field_at_surface_points = self.get_scalar_field_at_surface_points(Z_x, self.npf_op)

        # Masks
        if is_erosion:
            mask_e = Z_x > np.min(scalar_field_at_surface_points)
        else:
            mask_e = (~self.is_fault[n_series]) * np.ones_like(Z_x, dtype=bool)

        is_erosion_ = self.is_erosion[:n_series + 1]
        last_erode = np.argmax(np.concatenate(([1], is_erosion_)))

        is_onlap_or_fault = self.is_onlap[n_series] + self.is_fault[n_series]
        nsle = (nsle + is_onlap_or_fault) * is_onlap_or_fault * self.is_onlap[n_series - nsle]
        nsle_op = nsle

        if is_onlap:
            mask_o = Z_x > np.max(scalar_field_at_surface_points)
        else:
            mask_o = mask_matrix[n_series - 1, shift:x_to_interpolate_shape + shift]

        if self.is_fault[n_series]:
            mask_f = Z_x > np.min(scalar_field_at_surface_points)
        else:
            mask_f = np.zeros_like(Z_x, dtype=bool)

        # Compute blocks
        if not self.gradient:
            if compute_block_ctr:
                if is_finite:
                    block = self.compute_fault_block(
                        Z_x, scalar_field_at_surface_points,
                        self.values_properties_op[:, n_form_per_serie_0: n_form_per_serie_1 + 1],
                        n_series, grid
                    )
                else:
                    block = self.compute_formation_block(
                        Z_x, scalar_field_at_surface_points,
                        self.values_properties_op[:, n_form_per_serie_0: n_form_per_serie_1 + 1]
                    )
            else:
                block = block_matrix[n_series, :]
        else:
            if compute_block_ctr:
                block = self.compute_formation_block(
                    Z_x, scalar_field_at_surface_points,
                    self.values_properties_op[:, n_form_per_serie_0: n_form_per_serie_1 + 1]
                )
            else:
                block = block_matrix[n_series, :]

        # Update matrices
        weights_vector[len_w_0:len_w_1] = weights
        scalar_field_matrix[n_series, shift:x_to_interpolate_shape + shift] = Z_x
        block_matrix[n_series, :, shift:x_to_interpolate_shape + shift] = block
        fault_matrix[n_series, :, shift:x_to_interpolate_shape + shift] = block

        mask_matrix[n_series - 1:n_series, shift:x_to_interpolate_shape + shift] = mask_o
        mask_matrix[n_series - nsle_op:n_series, shift:x_to_interpolate_shape + shift] = np.cumprod(
            mask_matrix[n_series - nsle_op:n_series, shift:x_to_interpolate_shape + shift][::-1], axis=0
        )[::-1]
        mask_matrix[n_series, shift:x_to_interpolate_shape + shift] = mask_e

        idx_e = (self.is_fault * ~faults_relation_op.astype(bool))[:self.is_erosion.shape[0]]
        mask_matrix_f[idx_e, shift:x_to_interpolate_shape + shift] = mask_e + mask_f

        sfai[n_series, n_surface_op - 1] = scalar_field_at_surface_points

        return block_matrix, weights_vector, scalar_field_matrix, sfai, mask_matrix, \
            mask_matrix_f, fault_matrix, nsle

        
    import numpy as np

    def compute_forward_gravity_numpy(self, densities=None, pos_density=None):
        assert pos_density is not None or densities is not None, (
            "If a density block is not passed, you need to "
            "specify which interpolated value is density. See :class:`Surface`"
        )

        if densities is None:
            final_model, new_block, new_weights, new_scalar, new_sfai, new_mask = self.compute_series()
            densities = final_model[pos_density, self.lg0:self.lg1]
        else:
            final_model, new_block, new_weights, new_scalar, new_sfai, new_mask = None, None, None, None, None, None

        n_devices = int(densities.shape[0] / self.tz.shape[0])

        tz_rep = np.tile(self.tz, n_devices)

        # density times the component z of gravity
        grav = (densities * tz_rep).reshape((n_devices, -1)).sum(axis=1)

        return final_model, new_block, new_weights, new_scalar, new_sfai, new_mask, grav


    def compute_forward_gravity_pro_numpy(self, densities=None):
        n_devices = int(densities.shape[0] / self.tz.shape[0])

        tz_rep = np.tile(self.tz, n_devices)

        # density times the component z of gravity
        grav = (densities * tz_rep).reshape((n_devices, -1)).sum(axis=1)

        return grav



    import numpy as np

    def compute_forward_magnetics(self, k_vals):
        """
        Compute magnetics

        Args:
            k_vals: Susceptibility values per voxel [-] - varies per device! GemPy

        Returns:
            dT: total magnetic anomaly [nT]
        """

        def magnetic_direction(incl, decl):
            incl_rad = np.deg2rad(incl)
            decl_rad = np.deg2rad(decl)
            x = np.cos(incl_rad) * np.cos(decl_rad)
            y = np.cos(incl_rad) * np.sin(decl_rad)
            z = np.sin(incl_rad)
            return x, y, z

        if 'magnetics' in self.verbose:
            print("Sus. values:", k_vals)

        # induced magnetisation [T]
        J = k_vals * self.B_ext  # [k1dev1, ..., kndevn]

        # components
        dir_x, dir_y, dir_z = magnetic_direction(self.incl, self.decl)
        Jx = dir_x * J
        Jy = dir_y * J
        Jz = dir_z * J

        n_devices = int(k_vals.shape[0] / self.V.shape[1])
        if 'mag_devices' in self.verbose:
            print("n_devices:", n_devices)

        V = np.tile(self.V, (1, n_devices))  # repeat for each device

        # directional magnetic effect on one voxel (3.19)
        Tx = (Jx * V[0, :] + Jy * V[1, :] + Jz * V[2, :]) / (4 * self.pi)
        Ty = (Jx * V[1, :] + Jy * V[3, :] + Jz * V[4, :]) / (4 * self.pi)
        Tz = (Jx * V[2, :] + Jy * V[4, :] + Jz * V[5, :]) / (4 * self.pi)

        T2nT = 1e9  # convert to [nT]

        Tx = np.sum(Tx.reshape((n_devices, -1)), axis=1) * T2nT
        Ty = np.sum(Ty.reshape((n_devices, -1)), axis=1) * T2nT
        Tz = np.sum(Tz.reshape((n_devices, -1)), axis=1) * T2nT

        # total field anomaly in the direction of Earth's main field
        dT = Tx * dir_x + Ty * dir_y + Tz * dir_z
        return dT  # or return Tx, Ty, Tz if you want components too

# from dataclasses import dataclass
# from typing import Annotated

# import numpy as np

# from .encoders.converters import numpy_array_short_validator

# @dataclass
# class GeophysicsInputMagnetics:
#     tz: Annotated[np.ndarray,  numpy_array_short_validator]
#     densities: Annotated[np.ndarray,  numpy_array_short_validator]