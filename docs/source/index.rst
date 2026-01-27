.. Mineye-Terranigma documentation master file

==============================
Mineye-Terranigma Documentation
==============================

Overview
========

``Mineye-Terranigma`` is a research project for **probabilistic geological modeling**
and **geophysical inversion**. The project combines state-of-the-art implicit geological
modeling with Bayesian inference techniques.

**Core Capabilities:**

* **3D Geological Modeling**: Build implicit 3D models from structural data using `GemPy <https://www.gempy.org>`_
* **Geophysical Forward Modeling**: Compute gravity and magnetic responses from geological models
* **Probabilistic Inference**: Quantify uncertainty using `Pyro <https://pyro.ai>`_ and PyTorch
* **Satellite Image Segmentation**: Bayesian classification of remote sensing data
* **Joint Inversion**: Integrate multiple data sources for improved geological interpretations

----

Getting Started
===============

.. toctree::
   :maxdepth: 1

   installation

----

.. _basic-examples:

Basic Examples
==============

Foundational workflows for 3D geological modeling and geophysical forward modeling.
These examples introduce core concepts using real data from the **Tharsis mining district** in Spain.

.. include:: examples_basic/index.rst
   :start-after: .. thumbnail-parent-div-open
   :end-before: .. thumbnail-parent-div-close

.. toctree::
   :maxdepth: 2
   :hidden:

   examples_basic/index

----

.. _probabilistic-examples:

Probabilistic Modeling
======================

Advanced examples demonstrating uncertainty quantification, error propagation,
and Bayesian inversion techniques for geological and geophysical applications.

.. include:: examples_probabilistic/index.rst
   :start-after: .. thumbnail-parent-div-open
   :end-before: .. thumbnail-parent-div-close

.. toctree::
   :maxdepth: 2
   :hidden:

   examples_probabilistic/index

----

API Reference
=============

Detailed documentation of the Mineye-Terranigma Python API.

.. toctree::
   :maxdepth: 2

   api_reference

----

Key Features
============

Bayesian Inference with Pyro
----------------------------

The project uses PyTorch and Pyro for probabilistic programming, enabling:

* **Variational Inference (VI)**: Fast approximate posterior estimation
* **Hamiltonian Monte Carlo (HMC/NUTS)**: Accurate sampling for complex posteriors
* **GPU Acceleration**: Scalable inference for large-scale inversions

Geological Modeling with GemPy
------------------------------

Building on the GemPy framework, the project supports:

* **Implicit Surface Modeling**: Continuous 3D geological surfaces from sparse data
* **Structural Complexity**: Faults, unconformities, and erosive contacts
* **Forward Modeling**: Gravity and magnetic field computations
* **Uncertainty Quantification**: Probabilistic geological interpretations

----

Indices and Tables
==================

* :ref:`genindex`
* :ref:`search`
