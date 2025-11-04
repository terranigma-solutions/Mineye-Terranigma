"""
Configuration parameters for geological model examples.

This module centralizes all configuration parameters for different geological models,
including model extents, resolutions, surface mappings, and visualization settings.
"""
import numpy as np


# ========== THARSIS MODEL PARAMETERS ==========

class TharsisModelConfig:
    """Configuration parameters for the Tharsis geological model."""

    # Model identification
    PROJECT_NAME = 'simple_model'

    # Model extent [min_x, max_x, min_y, max_y, min_z, max_z]
    EXTENT = [-709521, -675558, 4526832, 4551949, -500, 505]

    # Model resolution [nx, ny, nz]
    RESOLUTION = [64, 64, 64]

    # Refinement level
    REFINEMENT = 4

    # Surface mapping
    SURFACE_MAPPING = {
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }

    # Visualization settings
    VIZ_CROSS_SECTION_COUNT = 5
    VIZ_VERTICAL_EXAGGERATION = 6
    VIZ_3D_VERTICAL_EXAGGERATION = 12
    VIZ_CELL_NUMBER = [5]
    VIZ_DIRECTION = ['y']

    # Plot settings
    SHOW_TOPOGRAPHY = True
    SHOW_LITH = True
    SHOW_BOUNDARIES = True
    SHOW_DATA = True
    SHOW_LEGEND = False


# ========== SORICOM MODEL PARAMETERS ==========

class SoricomModelConfig:
    """Configuration parameters for the Soricom geological model."""

    # Model identification
    PROJECT_NAME = 'soricom_model'

    # Model extent [min_x, max_x, min_y, max_y, min_z, max_z]
    # Based on actual data: X: 4442168-4442273, Y: 4588349-4588384, Z: 1558-1660
    # Adding buffer for modeling
    EXTENT = [4442150, 4442300, 4588330, 4588400, 1500, 1700]

    # Model resolution [nx, ny, nz]
    RESOLUTION = [50, 50, 50]

    # Refinement level
    REFINEMENT = 4

    # Surface mapping (fault must be first)
    SURFACE_MAPPING = {
        "Fault_Series": ["Main_Fault"],
        "overlying_lithology": ["overlying lithology"],
        "chromite_lense": ["chromite lense"],
    }

    # Fault configuration
    FAULT_GROUP_INDEX = 0  # Index of the fault in structural groups
    FAULT_RELATIONS_MATRIX = np.array([
        [False, True, True],   # Fault_Series: doesn't affect itself, affects both formations
        [False, False, False], # overlying_lithology: not a fault
        [False, False, False]  # chromite_lense: not a fault
    ])

    # Visualization settings
    VIZ_CROSS_SECTION_COUNT = 5
    VIZ_VERTICAL_EXAGGERATION = 1
    VIZ_CELL_NUMBER = [5]
    VIZ_DIRECTION = ['y']

    # Plot settings
    SHOW_TOPOGRAPHY = False
    SHOW_LITH = True
    SHOW_BOUNDARIES = True
    SHOW_DATA = True
    SHOW_LEGEND = False

