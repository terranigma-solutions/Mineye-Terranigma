import numpy as np

class TharsisModelConfig:
    """Configuration parameters for the Tharsis geological model."""

    PROJECT_NAME = 'simple_model'
    EXTENT = [-709521, -675558, 4526832, 4551949, -500, 505]
    RESOLUTION = [64, 64, 64]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }

    VIZ_CROSS_SECTION_COUNT = 5
    VIZ_VERTICAL_EXAGGERATION = 6
    VIZ_3D_VERTICAL_EXAGGERATION = 12
    VIZ_CELL_NUMBER = [5]
    VIZ_DIRECTION = ['y']
    SHOW_TOPOGRAPHY = True
    SHOW_LITH = True
    SHOW_BOUNDARIES = True
    SHOW_DATA = True
    SHOW_LEGEND = False

class SoricomModelConfig:
    """Configuration parameters for the Soricom geological model."""

    PROJECT_NAME = 'soricom_model'
    EXTENT = [4442150, 4442300, 4588330, 4588400, 1500, 1700]
    RESOLUTION = [50, 50, 50]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "Fault_Series": ["Main_Fault"],
        "overlying_lithology": ["overlying lithology"],
        "chromite_lense": ["chromite lense"],
    }

    FAULT_GROUP_INDEX = 0  # Index of the fault in structural groups
    FAULT_RELATIONS_MATRIX = np.array([
        [False, True, True],   # Fault_Series: doesn't affect itself, affects both formations
        [False, False, False], # overlying_lithology: not a fault
        [False, False, False]  # chromite_lense: not a fault
    ])

    VIZ_CROSS_SECTION_COUNT = 5
    VIZ_VERTICAL_EXAGGERATION = 1
    VIZ_CELL_NUMBER = [5]
    VIZ_DIRECTION = ['y']
    SHOW_TOPOGRAPHY = False
    SHOW_LITH = True
    SHOW_BOUNDARIES = True
    SHOW_DATA = True
    SHOW_LEGEND = False


class SoricomSimpleModelConfig:
    """Configuration parameters for the simple Soricom geological model (without fault)."""
    PROJECT_NAME = 'soricom_simple_model'
    EXTENT = [4442150, 4442300, 4588330, 4588400, 1500, 1700]
    RESOLUTION = [50, 50, 50]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "overlying_lithology": ["overlying lithology"],
        "chromite_lense": ["chromite lense"],
    }

    VIZ_CROSS_SECTION_COUNT = 5
    VIZ_VERTICAL_EXAGGERATION = 1
    VIZ_CELL_NUMBER = [5]
    VIZ_DIRECTION = ['y']
    SHOW_TOPOGRAPHY = False
    SHOW_LITH = True
    SHOW_BOUNDARIES = True
    SHOW_DATA = True
    SHOW_LEGEND = False
