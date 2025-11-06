Installation
============

Prerequisites
-------------

``Mineye-Terranigma`` requires Python 3.10 or higher.

Installation from Source
-------------------------

1. Clone the repository:

``$ git clone https://github.com/leguark/Mineye-Terranigma.git``

``$ cd Mineye-Terranigma``

2. Install dependencies:

``$ pip install -r requirements.txt``

3. Install the package in development mode:

``$ pip install -e .``

Dependencies
------------

Main dependencies include:

* arviz - Bayesian analysis and visualization
* torch - PyTorch for deep learning and automatic differentiation
* gempy, gempy_probability, gempy_viewer - Geological modeling
* matplotlib, numpy, pandas, scipy - Scientific computing
* rasterio, geopandas - Geospatial data processing
* scikit-image, scikit-learn - Machine learning and image processing

PyTorch with CUDA (Optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For GPU acceleration, install PyTorch with CUDA support:

``$ pip install torch --index-url https://download.pytorch.org/whl/cu118``

Visit https://pytorch.org/get-started/locally/ for platform-specific instructions.

Development Installation
------------------------

For development, also install testing and documentation dependencies:

``$ pip install pytest sphinx sphinx-gallery myst-parser``
