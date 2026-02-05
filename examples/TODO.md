ODO: Improvements for Mineye-Terranigma Examples (01-03)
This TODO list outlines ways to improve the didactic quality, visual appeal, and technical robustness of the example scripts in folders 01_basic_examples, 02_probabilistic_modeling, and 03_segmentation.
📁 01_basic_examples (Foundational Workflows)
•
[ ] 01_simple_tharsis_model.py
◦
Visuals: Add a 2D cross-section plot (using gpv.plot_2d) before the 3D plot to show the internal structure more clearly.
◦
Didactic: Add a small LaTeX block explaining the "Universal Co-Kriging" concept briefly.
◦
Robustness: Standardize the crop_to_extent logic to use a helper function instead of hardcoded values.
•
[ ] 02_complex_tharsis_model.py
◦
Didactic: Enhance the explanation of why we use a two-model workflow for erosive contacts (implicit vs. explicit handling of topological relationships).
◦
Code: Move the FORMATION_COLORS dictionary to a shared configuration or helper to avoid duplication in future examples.
◦
Visuals: Use pyvista to create a "cut-away" view in 3D to show the pluton cutting the layers.
•
[ ] 03_gravity_forward_model.py
◦
Cleanup: Remove commented-out lines like # topography_path = ... (line 328) to keep it clean.
◦
Visuals: Improve the residual plot (Panel 3) by using a diverging colormap (RdBu_r) with a fixed symmetric scale centered at zero.
◦
Theory: Add a note explaining the 1/r² nature of gravity and why the radius in set_centered_grid matters.
📁 02_probabilistic_modeling (Uncertainty & Inversion)
•
[ ] 02_error_propagation.py & 03_error_propagation_dips.py
◦
Fix: Resolve semantic errors regarding PyTorch/NumPy type mismatches (e.g., using .detach().numpy() where needed).
◦
Visuals: In plot_gempy, ensure the colors of the sampled lines match the geological unit they represent (even if transparent).
◦
Didactic: Add a diagram or text explaining "Prior Predictive Sampling" vs. "Posterior Sampling".
•
[ ] 04_gravity_inversion.py & 05_magnetics_inversion.py
◦
Advanced Visuals: Port the "Probability Density Fields" and "Information Entropy" plots from test_gravity_inversion.py. Show entropy in 3D by injecting it into the scalar field.
◦
Diagnostics: Implement the outlier detection logic for stations (Sigma analysis) to show which measurements are "problematic".
◦
Structure: Break 04_gravity_inversion.py into sub-sections using more # %% separators to make the gallery page readable.
◦
Fixes: Correct NUTSConfig arguments and Tensor/ndarray type mismatches in the likelihood functions.
◦
Interactive: Add instructions for downloading pre-computed .nc files.
📁 03_segmentation (Remote Sensing Integration)
•
[ ] 01_enmap_lith_segmentation.py
◦
Standardize: Replace the 'hacky' project root detection with paths.get_base_dir().
◦
Didactic: Add visual "Quicklooks" of MNF bands and explain their importance.
•
[ ] 02_enmap_data_extraction.py
◦
Advanced Sampling: Incorporate the "Boundary-Focused" and "Spatially Reduced" sampling strategies from test_enmap_preprocess.py.
◦
Visuals: Show a side-by-side comparison of raw GeoTIFF vs. sampled points on topography.
◦
Didactic: Explain "Distance Transform" visually with a dedicated subplot.
•
[ ] 03_enmap_gempy_comparison.py
◦
Label Mapping: Port the automated "Best Mapping" logic from test_enmap_likelihood.py to handle different class ID systems.
◦
Visuals: Standardize the 3-panel residual plot (EnMap vs GemPy vs Mismatch) seen in test_enmap_residuals.py.
◦
Metrics: Add a Confusion Matrix and accuracy-by-class reporting.
◦
Benchmarking: Add computation time reporting (e.g., "Model computation time: X seconds") inspired by test_structural_model.py.
🌐 General Improvements (Across All Files)
•
[ ] Path Handling: Replace all occurrences of os.path.join(base_dir, 'examples', 'Data', ...) with direct calls to paths.get_data_dir() or similar.
•
[ ] LaTeX: Ensure all mathematical formulas are wrapped in proper Sphinx math blocks (.. math::) for crisp rendering in documentation.
•
[ ] Captions: Add descriptive captions to every figure generated so they appear nicely in the Sphinx-Gallery.
•
[ ] Requirements: Add a # sphinx_gallery_thumbnail_number = X to every file to ensure the most representative plot is used in the gallery index.