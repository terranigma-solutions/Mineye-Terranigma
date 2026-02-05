import numpy as np

class TharsisModelConfig:
    """Configuration parameters for the Tharsis geological model."""

    PROJECT_NAME = 'simple_model'
    EXTENT = [-707521, -675558, 4526832, 4551949, -500, 505]
    RESOLUTION = [64, 64, 64]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "Tournaisian_Plutonites": ["Tournaisian Plutonites"],
    }

    class TharsisGravityConfig:
        """Configuration parameters for Tharsis gravity forward modeling and inversion."""

        # Model parameters
        DENSITY_PLUTONITES = 2.9  # kg/m³
        DENSITY_SEDIMENTARY_HOST = 2.3  # kg/m³
        GRAVITY_RESOLUTION = 20  # Number of gravity measurement points per axis

        # Centered grid parameters for gravity calculation
        CENTERED_GRID_RESOLUTION = np.array([10, 10, 15])
        CENTERED_GRID_RADIUS = np.array([5000, 5000, 5000])

        # Measurement grid options
        USE_ACTUAL_MEASUREMENT_LOCATIONS = True  # True: use gravity device locations, False: use regular grid

        # Normalization options
        NORMALIZE_DATA = True  # Enable/disable normalization
        NORMALIZATION_METHOD = 'minmax'  # Options: 'zscore', 'minmax', 'mean_center', 'relative'

        # Gravity data field name
        GRAVITY_FIELD_NAME = 'VALU_BOU267'  # Field name in the GeoJSON file (in mGal)

        # Visualization options
        SHOW_FORWARD_MODEL = True
        SHOW_COMPARISON_PLOTS = True
        SAVE_RESULTS = True


    class TharsisDataProcessingConfig:
        """Configuration parameters for Tharsis geological data processing."""

        # Data processing parameters
        DEFAULT_DIP = 10
        DEFAULT_AZIMUTH = 0
        USE_DEFAULT_AZIMUTH = False
        BOUNDARY_TOLERANCE = 800
        FORMATION_ID = 34
        SIMPLIFICATION_LEVEL = 0.9  # 0=no simplification, 1=maximum simplification

        # Shared color scheme for geological formations
        FORMATION_COLORS = {
            'Tournaisian Plutonites': '#e74c3c',        # Red - plutonite
            'Visean Shales': '#3498db',                 # Blue
            'Mid Devonian Siliciclastics': '#2ecc71',   # Green
            'Famennian Siliciclastics': '#f39c12',      # Orange - basement
            'basement': '#8B4513',                      # Brown - default basement
        }

        # Manual orientations to add at specific point IDs (after filtering/simplification)
        MANUAL_ORIENTATIONS_AT_POINTS = [11]  # Point 19 gets filtered out during boundary removal

        # Azimuth flips by orientation ID (180°)
        AZIMUTH_FLIP_BY_ID = {
            0: True,
            1: True,
            2: True,
            3: True,
            4: True,
            5: False,
            6: False
        }

        # Manual dip overrides by orientation ID
        MANUAL_DIP_BY_ID = {}  # Add manual dip values as needed: {3: 60, 7: 30}
        DENSITY_PLUTONITES = 2.9  # kg/m³
        DENSITY_SEDIMENTARY_HOST = 2.3  # kg/m³
        GRAVITY_RESOLUTION = 20  # Number of gravity measurement points per axis

        # Centered grid parameters for gravity calculation
        CENTERED_GRID_RESOLUTION = np.array([10, 10, 15])
        CENTERED_GRID_RADIUS = np.array([5000, 5000, 5000])

        # Measurement grid options
        USE_ACTUAL_MEASUREMENT_LOCATIONS = True  # True: use gravity device locations, False: use regular grid

        # Normalization options
        NORMALIZE_DATA = True  # Enable/disable normalization
        NORMALIZATION_METHOD = 'minmax'  # Options: 'zscore', 'minmax', 'mean_center', 'relative'

        # Gravity data field name
        GRAVITY_FIELD_NAME = 'VALU_BOU267'  # Field name in the GeoJSON file (in mGal)

        # Visualization options
        SHOW_FORWARD_MODEL = True
        SHOW_COMPARISON_PLOTS = True
        SAVE_RESULTS = True




class SoricomModelConfig:
    """Configuration parameters for the Soricom geological model."""

    PROJECT_NAME = 'soricom_model'
    # Updated extent to properly encompass all data points
    # Data ranges: X: 4441902-4442274, Y: 4588250-4588536, Z: 1558-1660
    EXTENT = [4441850, 4442350, 4588200, 4588400, 1500, 1700]
    RESOLUTION = [50, 50, 50]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "Fault_Series": ["Main_Fault"],
        "host_rock": ["host_rock"],
        "chromite_lense": ["chromite lense"],
    }

    FAULT_GROUP_INDEX = 0  # Index of the fault in structural groups
    FAULT_RELATIONS_MATRIX = np.array([
        [False, True, True],   # Fault_Series: doesn't affect itself, affects both formations
        [False, False, False], # host_rock: not a fault
        [False, False, False]  # chromite_lense: not a fault
    ])


class SoricomSimpleModelConfig:
    """Configuration parameters for the simple Soricom geological model (without fault)."""
    PROJECT_NAME = 'soricom_simple_model'
    EXTENT = [4441850, 4442350, 4588200, 4588400, 1700, 1500]
    RESOLUTION = [50, 50, 50]
    REFINEMENT = 4
    SURFACE_MAPPING = {
        "host_rock": ["host_rock"],
        "chromite_lense": ["chromite lense"],
    }


class SoricomErosiveModelConfig:
    PROJECT_NAME = 'soricom_erosive_model'
    # Data ranges: X: 4441902-4442274, Y: 4588250-4588536, Z: 1520-1660 (underlying layer extends deeper)
    EXTENT = [4441850, 4442350, 4588200, 4588400, 1700, 1450]
    RESOLUTION = [50, 50, 50]
    REFINEMENT = 4

    # Surface mapping with chromite lense as erosive intrusion
    # Order: Fault (first), then host rock, then erosive intrusion
    SURFACE_MAPPING = {
        "Fault_Series": ["Main_Fault"],
        "host_rock": ["host_rock"],
        "chromite_lense": ["chromite lense"],  # Erosive intrusion
    }

    FAULT_GROUP_INDEX = 0  # Index of the fault in structural groups
    CHROMITE_GROUP_INDEX = 2  # Index of the chromite lense (erosive) group

    # Fault relations matrix: Fault affects all formations
    FAULT_RELATIONS_MATRIX = np.array([
        [False, True,  True],   # Fault_Series: affects all formations
        [False, False, False],  # host_rock: not a fault
        [False, False, False],  # chromite_lense: not a fault (erosion handled separately)
    ])