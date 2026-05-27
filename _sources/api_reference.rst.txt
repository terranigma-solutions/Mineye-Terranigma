API Reference
=============

Mineye-Terranigma provides a modular API for geological modeling, geophysical
forward computation, probabilistic Bayesian inversion, and GIS data processing.

.. note::
   This API reference focuses on the main modules used in the example galleries.
   For detailed usage, please refer to the corresponding example scripts.

Mineye Package
--------------

The ``mineye`` package is organized into four main sub-packages:

GeoModel — Geological Modeling & Inversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mineye.GeoModel`` contains modules for implicit 3D geological model setup,
geophysical forward modeling, probabilistic inversion, and visualization.

**Core Modules:**

* ``mineye.GeoModel.model_one`` — Tharsis Plutonites model: setup, priors, likelihoods, joint inversion, inference diagnostics, and visualization
* ``mineye.GeoModel.geophysics`` — Forward model computation and quantile alignment of simulated to observed geophysical responses
* ``mineye.GeoModel.helper_methods`` — Data processing utilities (orientation handling, boundary removal, simplification)
* ``mineye.GeoModel.helper_plotter`` — Basic plotting functions for model geometry and geophysical sensors
* ``mineye.GeoModel.plotting.probabilistic_analysis`` — Advanced comparison and uncertainty visualization

**Key Workflow:**

1. ``model_setup.setup_geomodel()`` — Configure the GemPy structural model with PyTorch backend
2. ``geophysics.align_forward_to_observed()`` — Align forward responses to observed anomaly levels
3. ``probabilistic_model.set_priors()`` — Define prior distributions (interface dips, densities, susceptibilities)
4. ``probabilistic_model_likelihoods`` — Generate likelihood functions (gravity diagonal, hierarchical per-station, magnetic, EnMAP categorical)
5. ``joint_probabilistic_model`` — Combine likelihoods for joint inversion with balance diagnostics

**Inference Diagnostics:**

* ``inference_diagnostics.check_mcmc_quality()`` — MCMC convergence assessment
* ``inference_diagnostics.check_likelihood_balance()`` — Verify no single dataset dominates the joint posterior

**Visualization:**

* ``visualization.generate_gravity_uncertainty_plots()`` — Observed vs. forward comparison
* ``visualization.probability_fields_for()`` — Probability density maps and information entropy
* ``visualization.gempy_viz()`` — 3D model rendering


BayesianSegmentation — EO Data Processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mineye.BayesianSegmentation`` integrates hyperspectral satellite data with
Bayesian segmentation workflows.

**Core Modules:**

* ``full_workflow`` — End-to-end BaySeg integration pipeline
* ``EnMap_feature_extraction`` — EnMAP band parsing, wavelength extraction, and feature stacking
* ``prepare_data`` — GeoTIFF cropping and spatial sampling


GisHelpers — GIS Utilities
~~~~~~~~~~~~~~~~~~~~~~~~~~

``mineye.GisHelpers`` provides geospatial data processing utilities.

**Core Modules:**

* ``crs_utils`` — Coordinate reference system handling
* ``extractPointsFromMap`` — Point extraction from raster layers
* ``raster2obj`` — Raster to structured object conversion
* ``shapefiles2numpy`` — Shapefile to NumPy array conversion


Config — Configuration Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``mineye.config`` handles path resolution and model parameter configuration.

**Core Modules:**

* ``paths`` — Centralized path resolution for data, geophysical, topography, and model directories
* ``example_parameters`` — ``TharsisModelConfig`` and ``SoricomModelConfig`` dataclasses with extent, refinement, surface mapping, and formation colors