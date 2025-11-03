"""Tests for geological data processing functions."""
import sys
import os
# Add project root to path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import pandas as pd
from mineye.GeoModel import helper_methods


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


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

