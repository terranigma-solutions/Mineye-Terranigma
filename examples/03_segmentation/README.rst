Hyperspectral Segmentation
==========================

This gallery demonstrates the lithological segmentation workflow using EnMap hyperspectral data and Bayesian inference.

Gallery Overview
----------------

These examples cover the complete pipeline from raw hyperspectral cubes to georeferenced lithological maps with associated uncertainty. The workflow integrates remote sensing feature extraction with probabilistic classification.

Key Concepts
------------

* **Feature Extraction**: Dimensionality reduction (MNF) and extraction of geologically relevant spectral features (e.g., absorption depths).
* **Bayesian Segmentation**: Using spatial priors and spectral likelihoods to classify surface units.
* **Geological Integration**: Comparing surface segmentation results with 3D geological model predictions.
* **Uncertainty Quantification**: Estimating classification confidence and entropy across the scene.

Prerequisites
-------------

* Basic understanding of hyperspectral remote sensing.
* Familiarity with Bayesian classification concepts.
* Installed packages: ``rasterio``, ``scikit-image``, ``scipy``, in addition to the core Mineye dependencies.

Workflow Stages
---------------

1. **Preprocessing & Feature Extraction**: Preparing the hyperspectral data for classification.
2. **Segmentation**: Running the Bayesian engine to produce lithological labels.
3. **Validation**: Evaluating results against known geological structures and extracting points for further modeling.

.. note::
   The examples in this gallery require EnMap L2A products. Paths to data should be configured in ``mineye.config.paths`` or provided as environment variables as shown in the scripts.
