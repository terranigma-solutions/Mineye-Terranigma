import os
import time
import pytest

import gempy as gp
import gempy_viewer as gpv


def test_enmap_offset_gempy(simple_geo_model, topography_dir):
    topography_path = os.path.join(topography_dir, 'topo_reduced_sf0.1.tif')
    gp.set_topography_from_file(
        grid=simple_geo_model.grid,
        filepath=topography_path,
        crop_to_extent=[
                simple_geo_model.grid.extent[0],
                simple_geo_model.grid.extent[2],
                simple_geo_model.grid.extent[1],
                simple_geo_model.grid.extent[3]
        ]
    )

    start_time = time.time()

    simple_geo_model.interpolation_options.evaluation_options.number_octree_levels_surface = 3
    gp.compute_model(simple_geo_model)
    elapsed_time = time.time() - start_time
    print(f"\n⏱️  Model computation time: {elapsed_time:.2f} seconds")


def test_compare_enmap_gempy_labels(simple_geo_model, base_dir):
    """Compare GemPy predicted labels with EnMap extracted labels."""
    import numpy as np
    
    # 1. Load EnMap extracted points and labels
    xyz_path = os.path.join(base_dir, 'central_xyz.npy')
    labels_path = os.path.join(base_dir, 'central_labels.npy')
    
    if not os.path.exists(xyz_path) or not os.path.exists(labels_path):
        pytest.skip("EnMap extracted data not found. Run test_enmap_preprocess.py first.")
        
    xyz_central = np.load(xyz_path)
    labels_enmap = np.load(labels_path)
    
    print(f"\nLoaded {len(xyz_central)} points from EnMap extraction.")
    
    # 2. Set custom grid in GemPy model
    gp.set_custom_grid(simple_geo_model.grid, xyz_central)
    
    # 3. Compute GemPy model
    gp.compute_model(simple_geo_model)
    
    # 4. Get GemPy predicted labels at custom grid points
    # The custom grid results are in solutions.raw_arrays.custom
    labels_gempy = simple_geo_model.solutions.raw_arrays.custom
    
    # 5. Print surface information for mapping verification
    print("\nGemPy Structural Frame:")
    print(simple_geo_model.structural_frame)
    
    # 6. Compare labels
    # We might need to map EnMap labels to GemPy IDs or vice-versa
    # For now, let's just print them and see the alignment
    
    print("\nLabel Comparison:")
    print(f"EnMap labels (unique): {np.unique(labels_enmap)}")
    print(f"GemPy labels (unique): {np.unique(labels_gempy)}")
    print(f"GemPy ID to Name map: {simple_geo_model.structural_frame.element_id_name_map}")
    
    # 7. Find best mapping
    unique_enmap = np.unique(labels_enmap)
    
    mapping = {}
    for e_label in unique_enmap:
        # Find which GemPy label is most frequent for this EnMap label
        mask = (labels_enmap == e_label)
        relevant_gempy = labels_gempy[mask].astype(int)
        if len(relevant_gempy) > 0:
            # Use np.unique to find the most frequent value
            vals, counts = np.unique(relevant_gempy, return_counts=True)
            print(f"EnMap Label {e_label} -> GemPy label distribution: {dict(zip(vals, counts))}")
            best_g_label = vals[np.argmax(counts)]
            mapping[e_label] = best_g_label
            
    print(f"\nBest mapping (EnMap -> GemPy): {mapping}")
    
    # 8. Calculate accuracy with mapping
    mapped_enmap_labels = np.array([mapping[l] for l in labels_enmap])
    matches_mapped = np.sum(mapped_enmap_labels == labels_gempy)
    accuracy_mapped = matches_mapped / len(labels_enmap)
    print(f"Mapped match accuracy: {accuracy_mapped:.2%} ({matches_mapped}/{len(labels_enmap)})")
    
    # TODO: Use this accuracy to define a likelihood for optimization
    
    assert len(labels_gempy) == len(labels_enmap), "Number of labels must match"
