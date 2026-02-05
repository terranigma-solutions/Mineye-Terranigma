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
[ ] 04_gravity_inversion.py
◦
Structure: This file is very long (~1200 lines). Consider breaking it into sub-sections or using more # %% separators to make the generated gallery page more digestible.
◦
Fix: Check line 922 NUTSConfig(...) for the "unexpected argument" error (likely a version mismatch in gempy_probability).
◦
Interactive: Add a note about the RUN_SIMULATION flag and how users can download pre-computed .nc files to save time.
•
[ ] 05_magnetics_inversion.py
◦
Fix: Correct the expected type Tensor, got ndarray error in the likelihood function (line 367).
◦
Visuals: Add a plot showing the IGRF field vector (inclination/declination) to help users visualize the magnetic source of the TMI signal.
◦
Didactic: Explicitly contrast Susceptibility (Magnetic) vs. Density (Gravity) in a small table.
📁 03_segmentation (Remote Sensing Integration)
•
[ ] 01_enmap_lith_segmentation.py
◦
Standardize: Replace the 'hacky' project root detection (lines 63-80) with a standard call to paths.get_base_dir().
◦
Didactic: Add more visual "Quicklooks" of the MNF bands and explain why MNF is preferred over raw bands for geology.
◦
Robustness: Add a check for the existence of the specific EnMap folder and provide a download link or instruction if missing.
•
[ ] 02_enmap_data_extraction.py
◦
Visuals: Add a side-by-side comparison of the raw GeoTIFF vs. the sampled points on top of a topography map.
◦
Didactic: Explain the "Distance Transform" concept visually (show the distance map as a subplot).
•
[ ] 03_enmap_gempy_comparison.py
◦
Visuals: Create a "Confusion Matrix" or "Accuracy by Class" bar chart to provide a more quantitative summary than just a scatter plot.
◦
Integration: Add a concluding section on how this "Likelihood" is actually used in Example 04/05 to drive the inversion.
🌐 General Improvements (Across All Files)
•
[ ] Path Handling: Replace all occurrences of os.path.join(base_dir, 'examples', 'Data', ...) with direct calls to paths.get_data_dir() or similar.
•
[ ] LaTeX: Ensure all mathematical formulas are wrapped in proper Sphinx math blocks (.. math::) for crisp rendering in documentation.
•
[ ] Captions: Add descriptive captions to every figure generated so they appear nicely in the Sphinx-Gallery.
•
[ ] Requirements: Add a # sphinx_gallery_thumbnail_number = X to every file to ensure the most representative plot is used in the gallery index.