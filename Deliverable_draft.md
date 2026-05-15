# 1. Probabilistic Framework Modernization

## 1.1 Migration to PyTorch & Probabilistic API

As part of Task 2.2, the probabilistic geomodelling framework of GemPy was comprehensively modernized to support
scalable Bayesian inference and the integration of Earth Observation (EO)-derived constraints. This modernization was
necessary to align the probabilistic modelling capabilities with the new GemPy v3 engine and to enable robust,
extensible workflows for data fusion in geological modelling.

A complete refactoring of the probabilistic framework was carried out, most notably the migration from the legacy Theano
backend to the modern tensor-based PyTorch framework. This transition significantly improved performance, flexibility,
and long-term maintainability of the codebase. Deprecated dependencies were removed, and the overall structure of the
package was reorganized to support future extensions and integration into larger processing environments, including EO
processing pipelines and the LiquidEarth/IPOP platform.

To enable efficient and user-friendly probabilistic geomodelling, a high-level API was designed to simplify the setup,
execution, and analysis of Bayesian workflows. The API abstracts the underlying complexity of probabilistic inference,
allowing users to define models, priors, and likelihoods in a consistent and structured manner while maintaining full
flexibility for advanced use cases.

A central outcome of this effort is the introduction of a unified data model, enabling consistent handling of
heterogeneous data sources. Core model elements such as scalar fields, geological surfaces, and observation datasets are
now handled through unified data structures. This ensures consistency across workflows and facilitates the integration
of geological inputs (interface points, orientations), geophysical observations, and EO-derived constraints within a
single probabilistic framework.

## 1.2 Performance & Numerical Stability

A systematic analysis of the previous probabilistic implementation identified several sources of numerical instability
and performance limitations, primarily the use of non-vectorized operations and inconsistent parameter handling.

To address these limitations, critical components of the framework were refactored and implemented using fully
vectorized PyTorch operations. This resulted in more stable and efficient execution of probabilistic inference
workflows, particularly during the probabilistic learning phase and iterative Bayesian methods. Enhanced methods for
scalar field-based segmentation of geological features (such as faults and layers) were introduced, which are critical
for stable model convergence.

Memory usage was further optimized through improved data handling strategies, enabling more efficient storage and
processing of intermediate results. This is highly relevant for large-scale applications where memory constraints become
a limiting factor. Furthermore, the migration to the PyTorch backend prepares the framework for native GPU acceleration,
allowing computationally intensive operations to be executed on modern hardware architectures, providing a clear pathway
for future performance gains in demanding MINEYE use cases.

# 2. Forward Modeling Operators

Forward modelling operators establish a direct link between geological models and geophysical observables, forming the
basis for probabilistic inversion. They translate the implicit geological representation into physically meaningful
property fields, enabling the simulation of measurable geophysical responses. Within this work, forward models were
implemented for potential field methods (gravity and magnetics), which can be computed efficiently from volumetric
property distributions without solving partial differential equations. This efficiency makes them exceptionally well-
suited for repeated forward evaluations required in probabilistic workflows.

## 2.1 Gravity Formulation

A gravity forward modelling approach was implemented based on density contrasts between geological units. The implicit
formulation allows a continuous representation of subsurface structures, which are subsequently translated into
discretized volumetric density fields for geophysical computation. Lithological model outputs are converted into density
distributions, supporting spatially varying densities within individual units to represent intra-unit heterogeneity and
improve realism.

The gravity signal is computed using a discretized volumetric representation of the model domain. The forward simulation
can be mathematically expressed as:

$$ \mathbf{d}_{sim} = \mathbf{G} oldsymbol{ho} $$

where $\mathbf{d}_{sim}$ represents the simulated gravity anomalies, $\mathbf{G}$ is the sensitivity matrix (or forward
operator) accounting for the geometry of the discretized cells, and $oldsymbol{ ho}$ is the density contrast
distribution derived from the implicit geological model. This approach provides the flexibility to handle complex
geometries while maintaining computational efficiency on regional-scale datasets and irregular observation grids. The
forward model is tightly integrated into the probabilistic framework via a likelihood function, comparing simulated
gravity anomalies directly to observed data.

## 2.2 Magnetics (TMI) Formulation

A forward modelling approach for magnetic data was implemented based on the distribution of magnetic susceptibility
within geological units. Geological model outputs are mapped to susceptibility values to construct volumetric property
models.

Magnetic anomalies (such as Total Magnetic Intensity - TMI) are computed from this volumetric model, accounting
primarily for induced magnetization driven by the ambient magnetic field. Where applicable, the formulation can be
extended to include remanent magnetization for specific geological settings. Similar to gravity, the implementation
supports multi-scale datasets and accommodates spatially heterogeneous distributions, optimized for efficient repeated
evaluation.

## 2.3 Categorical Likelihoods & Ordinal Probabilities (EO-derived data)

To integrate surface lithological information into the probabilistic geomodelling framework, categorical likelihood
functions were developed. This approach allows for the direct incorporation of EO-derived classifications, such as
lithological classes obtained from hyperspectral segmentation workflows, as constraints within the Bayesian inference.

The formulation distinguishes between categorical likelihoods for discrete class assignments and ordinal probabilities
for ordered geological relationships (e.g., stratigraphic sequences). Crucially, to account for uncertainty in EO
classification results, **soft probabilistic representations** are used instead of hard labels. Class membership
probabilities mapped onto the geological model domain define a categorical likelihood, ensuring that classification
uncertainty is coherently propagated from EO data processing into the resulting geological models.

*(Note: Within MINEYE, workflows often derive these categorical constraints using Bayesian segmentation approaches like BaySeg. While BaySeg serves as a powerful demonstration for generating these soft probabilities from hyperspectral data, the implemented probabilistic API is fully agnostic and can incorporate any advanced, site-specific, or machine learning-driven lithological classifications.)*

# 3. Bayesian Inference Strategy

Bayesian formulations were developed to enable the simultaneous integration of geological, geophysical, and EO-derived
information within a unified probabilistic framework. Instead of yielding a single deterministic model, this approach
explicitly captures uncertainties in data and model assumptions, sampling posterior distributions to generate ensembles
of plausible geological realizations.

## 3.1 Inference and Sampling Methodology

The inversion process is driven by combining prior geological knowledge with likelihood functions defined by the misfit
between observed data and the simulated responses derived from the implicit geological model. The foundational Bayesian
formulation is:

$$ P(\mathbf{m} | \mathbf{d}_{obs}) \propto P(\mathbf{d}_{obs} | \mathbf{m}) P(\mathbf{m}) $$

where $\mathbf{m}$ represents the geological model parameters (e.g., interface positions, densities), $\mathbf{d}_{obs}$
represents the observed geophysical or EO data, $P(\mathbf{m})$ is the prior distribution, and $P(\mathbf{d}_{obs} |
\mathbf{m})$ is the likelihood function. For potential field geophysics (gravity and magnetics), the likelihood is
defined as a Multivariate Normal distribution governed by a covariance matrix $\mathbf{\Sigma}$ that captures
measurement noise and model spatial correlation:

$$ P(\mathbf{d}_{obs} | \mathbf{m}) = \mathcal{N}(\mathbf{d}_{sim}(\mathbf{m}), \mathbf{\Sigma}) $$

Inference is carried out using advanced Markov Chain Monte Carlo (MCMC) techniques. Specifically, the framework utilizes
the **No-U-Turn Sampler (NUTS)** implemented via the **Pyro** probabilistic programming library. During the sampling
process (typically involving ~1000 samples following an initial warmup phase), each iteration systematically updates the
geological parameters, recomputes the implicit geological model, and evaluates the forward response against the combined
likelihoods. This rigorous approach ensures the robust handling of measurement noise and data uncertainty, leading to
progressively improved model realizations.

## 3.2 Joint Inversion Framework

A framework for Bayesian joint inversion was established to simultaneously evaluate multiple data types. A joint
likelihood function is constructed by mathematically integrating the individual contributions from gravity data
($\mathbf{d}_{grav}$), magnetic data ($\mathbf{d}_{mag}$), and EO-derived constraints ($\mathbf{d}_{EO}$). Assuming
conditional independence between the different observation types given the geological model, the joint posterior is
formulated as:

$$ P(\mathbf{m} | \mathbf{d}_{grav}, \mathbf{d}_{mag}, \mathbf{d}_{EO}) \propto P(\mathbf{d}_{grav} | \mathbf{m}) P(\mathbf{d}_{mag} | \mathbf{m}) P(\mathbf{d}_{EO} | \mathbf{m}) P(\mathbf{m}) $$

By demanding that the updated geological parameters simultaneously satisfy all available data sources, the joint
inversion significantly improves model constraints and reduces the inherent ambiguity associated with single-data
inversions (e.g., the non-uniqueness of potential fields).

## 3.3 Uncertainty Estimation

Uncertainty across the resulting ensemble of geological realizations is quantified using multiple complementary
measures, including variance fields, information entropy maps, and ensemble-based statistics. These spatially resolved
metrics allow for the identification of well-constrained regions strongly supported by data, versus zones of high
uncertainty requiring alternative modelling assumptions or further data acquisition, thereby providing a foundation for
risk-aware decision-making.

# 4. Application Case Study: Tharsis AOI 1

To demonstrate the capabilities of the modernized framework, a series of Bayesian inversions were performed on the
Tharsis AOI 1 region. The inversions build upon an initial GemPy structural model of the area, targeting key geometrical
parameters (such as interface positions) and physical properties (e.g., density contrasts and magnetic
susceptibilities). Prior distributions were defined based on the initial conceptual model and available petrophysical
information.

## 4.1 Gravity Data Constraints

A Bayesian gravity inversion was performed using available regional gravity observations to constrain the large-scale
structural model.

**Outcomes:** The inversion demonstrated that gravity data provides meaningful constraints on large-scale density
contrasts and deep interface geometries, despite limited direct subsurface information. Posterior analysis revealed a
marked reduction of uncertainty in regions well-covered by the gravity data. However, as expected with potential field
methods, ambiguity remained in areas where different deep geometric configurations produced similar surface gravity
responses.

*[Figure: initial model vs constrained model?]*
*[Figure: constrained model, information entropy]*

## 4.2 Magnetic (TMI) Data Constraints

A subsequent inversion utilized Total Magnetic Intensity (TMI) data to further refine the subsurface interpretations,
focusing on the estimation of susceptibility distributions and their structural boundaries.

**Outcomes:** The TMI data provided critical additional constraints, particularly in areas exhibiting significant
susceptibility contrasts. Compared to the gravity inversion, the magnetic data offered higher sensitivity to specific
lithological variations and shallower structural features. Improved constraints were observed in areas with a strong
magnetic signal, effectively complementing the gravity-based results and highlighting regions requiring joint
interpretation.

## 4.3 EnMap EO-derived Constraints

To evaluate the integration of surface data, a Bayesian inversion was performed using hyperspectral data from the EnMAP
satellite. Probabilistic categorical constraints were derived from the imagery using a Bayesian segmentation workflow
(e.g., BaySeg), which included spectral preprocessing, masking, and clustering into lithological classes with soft
membership probabilities.

**Outcomes:** The likelihood evaluated the agreement between the model-predicted surface lithologies and the EnMap-
derived class probabilities. The results demonstrated that EO constraints provide highly valuable, high-resolution
information regarding surface geology. This substantially improved the spatial consistency of the geological model near
the surface and reduced structural ambiguity where subsurface data was sparse. Naturally, the influence of the EO
constraints decreased with depth.

## 4.4 Joint Inversion Outcomes

*(Note: The Joint interpretation capitalizes on the combined framework detailed in Section 3.2)*

To maximize the reduction of uncertainty, a Joint Inversion was conducted combining both the Gravity observations and
the EnMap EO-derived categorical constraints within a single unified likelihood function.

**Outcomes:** The joint interpretation capitalized on the complementary strengths of both datasets. The EnMap
hyperspectral constraints anchored the near-surface lithological boundaries with high confidence, preventing the
potential field inversion from introducing geologically unreasonable structures at the surface. Simultaneously, the
gravity data constrained the deeper volumetric density distributions and interface geometries where the EO data lacked
penetration. This combined evidence approach resulted in a more robust characterization of the subsurface, significantly
reducing the non-uniqueness and overall spatial uncertainty compared to utilizing either dataset in isolation.

*[FIGURE: Lith Segmentation, derived classes]*
