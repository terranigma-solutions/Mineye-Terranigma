"""
Bayesian Magnetic Inversion: Soricom Chromite Lens
===================================================
This tutorial demonstrates Bayesian inversion of Total Magnetic Intensity (TMI)
data using the Soricom fault structural model. It extends the magnetic inversion
methodology to a more complex geological setting with a fault and a chromite
lens hosted in ultramafic rocks.

**What Makes the Soricom Model Different?**

Unlike the simple Tharsis plutonite model (:ref:`sphx_glr_02_probabilistic_modeling_05_magnetics_inversion.py`),
the Soricom model features:

1. **A fault**: The Main_Fault truncates all formations (fault-first structural
   frame ordering)
2. **A chromite lens**: A thin, high-susceptibility target embedded in host rock
3. **Smaller extent**: ~500 m × 350 m area at UTM zone 34N coordinates
4. **Higher resolution**: Octree refinement level 5 on a much smaller domain

**Data and Preprocessing**

TMI data is extracted from a merged B1B2 raster at the Soricom prospect.
The raw measurements are in **nanoTesla (nT)** and include Earth's background
magnetic field (~47,500 nT at this location). Since the GemPy forward model
computes **TMI anomalies** (deviations from the background field), we subtract
the IGRF (International Geomagnetic Reference Field) intensity before inversion:

.. math::

    TMI_{\text{anomaly}} = TMI_{\text{measured}} - IGRF_{\text{intensity}}


For comprehensive theory on Bayesian inversion, MCMC, and hierarchical
likelihoods, see :ref:`sphx_glr_02_probabilistic_modeling_04_gravity_inversion.py`.
"""

