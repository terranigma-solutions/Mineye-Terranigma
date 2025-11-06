"""Tests for geological data processing functions."""
import sys
import os
import pandas as pd
from mineye.GeoModel import helper_methods
from mineye.config import paths
from mineye.config.example_parameters import TharsisDataProcessingConfig

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


def test_recreate_tharsis_temp_files():
    """
    Recreate the Tharsis temp files (orientations_mod.csv and points_mod.csv)
    from the original plutonic contact points data.

    This test processes the raw plutonic data and saves the processed files
    to the temp_inputs directory, exactly as in the original workflow.
    """
    print("\n" + "="*60)
    print("RECREATING THARSIS TEMP FILES FROM ORIGINAL DATA")
    print("="*60)

    # Load original full points
    original_contact_points_path = paths.get_plutonic_contact_points_path()
    print(f"\nLoading original data from: {original_contact_points_path}")

    if not os.path.exists(original_contact_points_path):
        raise FileNotFoundError(
            f"Original contact points file not found: {original_contact_points_path}\n"
            f"Please ensure the plutonic data exists in the Plutonic_Data directory."
        )

    original_points = pd.read_csv(original_contact_points_path)
    points_df = original_points.copy()

    print(f"  Loaded {len(points_df)} original contact points")

    # Process the geological data using configuration
    print(f"\nProcessing geological data...")
    print(f"  Default dip: {TharsisDataProcessingConfig.DEFAULT_DIP}°")
    print(f"  Boundary tolerance: {TharsisDataProcessingConfig.BOUNDARY_TOLERANCE}m")
    print(f"  Simplification level: {TharsisDataProcessingConfig.SIMPLIFICATION_LEVEL * 100}%")

    orientations_df, points_df = helper_methods.process_geological_data(
        points_df=points_df,
        default_dip=TharsisDataProcessingConfig.DEFAULT_DIP,
        default_azimuth=TharsisDataProcessingConfig.DEFAULT_AZIMUTH,
        use_default_azimuth=TharsisDataProcessingConfig.USE_DEFAULT_AZIMUTH,
        boundary_tolerance=TharsisDataProcessingConfig.BOUNDARY_TOLERANCE,
        formation_id=TharsisDataProcessingConfig.FORMATION_ID,
        simplification_level=TharsisDataProcessingConfig.SIMPLIFICATION_LEVEL,
        manual_orientations_at_points=TharsisDataProcessingConfig.MANUAL_ORIENTATIONS_AT_POINTS,
        azimuth_flip_by_id=TharsisDataProcessingConfig.AZIMUTH_FLIP_BY_ID,
        manual_dip_by_id=TharsisDataProcessingConfig.MANUAL_DIP_BY_ID
    )

    print(f"\nProcessing complete:")
    print(f"  Generated {len(orientations_df)} orientations")
    print(f"  Retained {len(points_df)} contact points")
    print(f"  Reduction: {len(original_points)} → {len(points_df)} points "
          f"({(1 - len(points_df)/len(original_points))*100:.1f}% reduction)")

    # Get output paths
    mod_or_path = paths.get_orientations_path()
    mod_pts_path = paths.get_points_path()

    # Ensure temp directory exists
    tmp_dir = paths.get_tmp_dir()
    os.makedirs(tmp_dir, exist_ok=True)

    # Save processed data
    print(f"\nSaving processed data:")
    print(f"  Orientations → {mod_or_path}")
    print(f"  Points       → {mod_pts_path}")

    orientations_df.to_csv(mod_or_path, index=False)
    points_df.to_csv(mod_pts_path, index=False)

    print("\n✓ Temp files successfully recreated!")
    print("="*60)

    # Verify outputs
    assert isinstance(orientations_df, pd.DataFrame), "Should return orientations DataFrame"
    assert isinstance(points_df, pd.DataFrame), "Should return points DataFrame"
    assert 'id' in orientations_df.columns, "Orientations should have 'id' column"
    assert 'id' in points_df.columns, "Points should have 'id' column"
    assert len(orientations_df) > 0, "Should have generated some orientations"
    assert len(points_df) > 0, "Should have some points remaining"

    # Verify files were created
    assert os.path.exists(mod_or_path), f"Orientations file should exist: {mod_or_path}"
    assert os.path.exists(mod_pts_path), f"Points file should exist: {mod_pts_path}"

    # Verify file contents
    saved_orientations = pd.read_csv(mod_or_path)
    saved_points = pd.read_csv(mod_pts_path)
    assert len(saved_orientations) == len(orientations_df), "Saved orientations should match"
    assert len(saved_points) == len(points_df), "Saved points should match"

    print(f"\n✓ All verification checks passed!")


def test_process_geological_data_basic():
    """Test that process_geological_data returns valid dataframes."""
    # Create simple test data
    test_points = pd.DataFrame({
        'X': [-700000 + i * 1000 for i in range(20)],
        'Y': [4530000 + i * 1000 for i in range(20)],
        'Z': [300] * 20,
        'formation': ['Tournaisian Plutonites'] * 20,
        'formation_id': [34] * 20
    })

    # Process the data
    orientations_df, points_df = helper_methods.process_geological_data(
        points_df=test_points,
        default_dip=10,
        default_azimuth=0,
        use_default_azimuth=False,
        boundary_tolerance=800,
        formation_id=34,
        simplification_level=0.5,
        manual_orientations_at_points=None,
        azimuth_flip_by_id=None,
        manual_dip_by_id=None
    )

    # Verify outputs
    assert isinstance(orientations_df, pd.DataFrame), "Should return orientations DataFrame"
    assert isinstance(points_df, pd.DataFrame), "Should return points DataFrame"
    assert 'id' in orientations_df.columns, "Orientations should have 'id' column"
    assert 'id' in points_df.columns, "Points should have 'id' column"
    assert len(orientations_df) > 0, "Should have generated some orientations"
    assert len(points_df) > 0, "Should have some points remaining"


def test_process_geological_data_with_manual_adjustments():
    """Test manual adjustments are applied correctly."""
    # Create simple test data
    test_points = pd.DataFrame({
        'X': [-700000 + i * 1000 for i in range(20)],
        'Y': [4530000 + i * 1000 for i in range(20)],
        'Z': [300] * 20,
        'formation': ['Tournaisian Plutonites'] * 20,
        'formation_id': [34] * 20
    })

    # Process with manual adjustments
    orientations_df, points_df = helper_methods.process_geological_data(
        points_df=test_points,
        default_dip=10,
        default_azimuth=0,
        use_default_azimuth=False,
        boundary_tolerance=500,
        formation_id=34,
        simplification_level=0.3,
        manual_orientations_at_points=[0, 5],
        azimuth_flip_by_id={0: True},
        manual_dip_by_id={1: 45}
    )

    # Verify manual adjustments were applied
    assert len(orientations_df) > 0, "Should have orientations after processing"

    # Check if dip override was applied (if ID 1 exists)
    if 1 in orientations_df['id'].values:
        dip_at_id_1 = orientations_df[orientations_df['id'] == 1]['dip'].values[0]
        assert dip_at_id_1 == 45, f"Manual dip should be 45, got {dip_at_id_1}"


def test_process_geological_data_simplification():
    """Test that simplification reduces data appropriately."""
    # Create larger test dataset
    test_points = pd.DataFrame({
        'X': [-700000 + i * 500 for i in range(100)],
        'Y': [4530000 + i * 500 for i in range(100)],
        'Z': [300] * 100,
        'formation': ['Tournaisian Plutonites'] * 100,
        'formation_id': [34] * 100
    })

    # Process with high simplification
    orientations_high, points_high = helper_methods.process_geological_data(
        points_df=test_points.copy(),
        simplification_level=0.9,
        formation_id=34
    )

    # Process with low simplification
    orientations_low, points_low = helper_methods.process_geological_data(
        points_df=test_points.copy(),
        simplification_level=0.1,
        formation_id=34
    )

    # High simplification should result in fewer points
    assert len(points_high) < len(points_low), \
        "High simplification should result in fewer points"
    assert len(orientations_high) < len(orientations_low), \
        "High simplification should result in fewer orientations"

