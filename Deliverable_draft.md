# Introduction
Within the MINEYE project, probabilistic methods are developed to enable the integration of Earth Observation (EO) data into the construction and updating of 3D geological models. This requires a framework that is both robust and computationally efficient, while remaining flexible and accessible for practical application in exploration workflows.

The methodological developments build upon the scientific open-source packages GemPy [1] and BaySeg [2]. GemPy provides the foundation for implicit geological modelling and probabilistic inversion, while BaySeg enables location aware Bayesian segmentation of arbitrary data in 2D and 3D. Although previous work has demonstrated the potential of probabilistic inversion using GemPy, these approaches have largely remained exploratory and limited in scalability and maturity.

In this work, both frameworks were substantially modernized and extended to meet the requirements of MINEYE. This includes improvements in numerical stability, performance, data integration, and workflow automation. In addition, new formulations were developed to consistently incorporate EO-derived constraints, geophysical data, and uncertainty estimates into a unified probabilistic modelling framework. Together, these developments establish the basis for general, scalable application of probabilistic geomodelling and data fusion within the project.

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

$$ \mathbf{d}_{sim} = \mathbf{G}\boldsymbol{\rho} $$

where $\mathbf{d}_{sim}$ represents the simulated gravity anomalies, $\mathbf{G}$ is the sensitivity matrix (or forward
operator) accounting for the geometry of the discretized cells, and $\boldsymbol{\rho}$ is the density contrast
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

In the implemented EnMAP workflow, this mapping is performed through a differentiable ordinal-probability formulation:
continuous GemPy scalar fields at EO observation points are transformed into class probabilities using sigmoid
transitions across geological boundaries. A temperature parameter controls transition sharpness (from soft transitions
to near-discrete boundaries), preserving differentiability required by NUTS while still representing classification
uncertainty near contacts.

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
\mathbf{m})$ is the likelihood function. For potential field geophysics (gravity and magnetics), this can be written in
generic form as a multivariate Gaussian likelihood:

$$ P(\mathbf{d}_{obs} | \mathbf{m}) = \mathcal{N}(\mathbf{d}_{sim}(\mathbf{m}), \mathbf{\Sigma}) $$

Inference is carried out using advanced Markov Chain Monte Carlo (MCMC) techniques. Specifically, the framework utilizes
the **No-U-Turn Sampler (NUTS)** implemented via the **Pyro** probabilistic programming library. During the sampling
process (typically involving ~1000 samples following an initial warmup phase), each iteration systematically updates the
geological parameters, recomputes the implicit geological model, and evaluates the forward response against the combined
likelihoods. This rigorous approach ensures the robust handling of measurement noise and data uncertainty, leading to
progressively improved model realizations.

To improve inference robustness and transparency, the workflow explicitly includes prior predictive checks (to verify
that prior assumptions can plausibly reproduce the observed data domain) and posterior predictive checks (to assess fit
quality, residual structure, and potential model misspecification after conditioning on observations).

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

The practical implementation also incorporates a dedicated likelihood-balance diagnostic before full sampling to verify
that no single dataset numerically dominates the joint posterior. This step is especially important when combining
continuous potential-field signals with categorical EO constraints that have different scales and noise structures.

At model level, the workflow supports multi-grid evaluation within one probabilistic run (e.g., centered grids for
gravity responses and custom grids at EO pixel coordinates), allowing both data types to be enforced consistently on
their native observation supports.

## 3.3 Uncertainty Estimation

Uncertainty across the resulting ensemble of geological realizations is quantified using multiple complementary
measures, including variance fields, information entropy maps, and ensemble-based statistics. These spatially resolved
metrics allow for the identification of well-constrained regions strongly supported by data, versus zones of high
uncertainty requiring alternative modelling assumptions or further data acquisition, thereby providing a foundation for
risk-aware decision-making.

For geophysical datasets, uncertainty treatment is further strengthened by hierarchical per-station noise modelling,
where station-specific noise levels are inferred jointly with geological parameters. This enables automatic detection of
stations with anomalously high inferred noise (potential outliers or locally unmodelled complexity), reducing their
undue influence on posterior structure.

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

From an implementation perspective, the gravity workflow also includes alignment of forward responses to observed anomaly
levels (regional-residual handling via reference alignment) and hierarchical station-noise estimation, both of which
improve numerical stability and robustness against heterogeneous data quality.

Key geological takeaway: gravity primarily constrains deep and regional-scale geometry, while residual ambiguity remains
where multiple deep configurations can explain similar anomaly patterns.

*[Figure suggestion: Prior vs posterior gravity response maps at station locations. Panels show observed gravity, posterior mean forward gravity, residuals (observed minus predicted), and correlation plot with 1:1 line.]*
*[Figure suggestion: Posterior geological uncertainty summary with probability density fields and information entropy map; low entropy highlights well-constrained structures, high entropy highlights ambiguous zones.]*

## 4.2 Magnetic (TMI) Data Constraints

A subsequent inversion utilized Total Magnetic Intensity (TMI) data to further refine the subsurface interpretations,
focusing on the estimation of susceptibility distributions and their structural boundaries.

**Outcomes:** The TMI data provided critical additional constraints, particularly in areas exhibiting significant
susceptibility contrasts. Compared to the gravity inversion, the magnetic data offered higher sensitivity to specific
lithological variations and shallower structural features. Improved constraints were observed in areas with a strong
magnetic signal, effectively complementing the gravity-based results and highlighting regions requiring joint
interpretation.

The probabilistic magnetic inversion mirrors the gravity workflow with hierarchical per-station noise treatment and
posterior predictive residual diagnostics, enabling direct identification of magnetically inconsistent stations and
improved confidence assessment of susceptibility-driven interpretations.

Key geological takeaway: magnetic inversion adds sharper constraints on susceptibility-driven and relatively shallower
lithological boundaries, complementing the broader depth sensitivity of gravity.

## 4.3 EnMAP EO-derived Constraints

To evaluate the integration of surface data, a Bayesian inversion was performed using hyperspectral data from the EnMAP
satellite. Probabilistic categorical constraints were derived from the imagery using a Bayesian segmentation workflow
(e.g., BaySeg), which included spectral preprocessing, masking, and clustering into lithological classes with soft
membership probabilities.

**Outcomes:** The likelihood evaluated the agreement between the model-predicted surface lithologies and the EnMAP-
derived class probabilities. The results demonstrated that EO constraints provide highly valuable, high-resolution
information regarding surface geology. This substantially improved the spatial consistency of the geological model near
the surface and reduced structural ambiguity where subsurface data was sparse. Naturally, the influence of the EO
constraints decreased with depth.

The EnMAP inversion examples further show that using EO observations on custom surface grids and differentiable ordinal
likelihoods allows direct gradient-based assimilation of classification information without collapsing uncertainty into
hard labels.

Key geological takeaway: EO constraints are most powerful for anchoring near-surface contacts and lithological
organization, but their direct control decreases with depth.

## 4.4 Joint Inversion Outcomes

*(Note: The Joint interpretation capitalizes on the combined framework detailed in Section 3.2)*

To maximize the reduction of uncertainty, a Joint Inversion was conducted combining both the Gravity observations and
the EnMAP EO-derived categorical constraints within a single unified likelihood function.

**Outcomes:** The joint interpretation capitalized on the complementary strengths of both datasets. The EnMAP
hyperspectral constraints anchored the near-surface lithological boundaries with high confidence, preventing the
potential field inversion from introducing geologically unreasonable structures at the surface. Simultaneously, the
gravity data constrained the deeper volumetric density distributions and interface geometries where the EO data lacked
penetration. This combined evidence approach resulted in a more robust characterization of the subsurface, significantly
reducing the non-uniqueness and overall spatial uncertainty compared to utilizing either dataset in isolation.

Operationally, this combined workflow benefits from explicit likelihood-balance checks and shared-parameter inference
across multi-grid observations, which improves confidence that both deep geophysical structure and surface EO evidence
are jointly honored in the posterior ensemble.

Key geological takeaway: joint inversion delivers the most balanced structural interpretation by combining shallow
surface control from EnMAP with deeper volumetric constraints from gravity.

Current limitations: uncertainty remains higher where observation coverage is sparse, potential-field non-uniqueness is
not fully eliminated at depth, and results remain conditional on the assumed structural model parameterization and prior
ranges.

*[Figure suggestion: Joint inversion summary panel showing EO class probability map at surface, gravity fit improvement, and posterior entropy reduction relative to single-data inversions.]*

# 5. Application Case Study: Soricom Chromite Lens

To further evaluate the flexibility of the framework in structurally complex, mineral-exploration settings,
a Bayesian magnetic inversion was performed on the Soricom prospect, targeting a thin chromite lens hosted
within faulted ultramafic rocks. Unlike the regional-scale Tharsis model, the Soricom case represents a
deposit-scale problem where the primary uncertainties concern the precise position, thickness, and geometry
of a mineralized lens rather than broad regional interface geometries.

## 5.1 Geological Setting and Model Setup

The Soricom structural model covers a compact domain of approximately 500 m × 350 m at octree refinement level 5,
providing the high spatial resolution needed to resolve a thin chromite lens. The stratigraphy includes a
Main_Fault that truncates all formations (fault-first structural frame ordering), a host ultramafic rock,
and the chromite lens as a discrete, high-susceptibility target. Four geological units are defined, with
susceptibility values dominated by the strong magnetic contrast of the chromite lens (~0.5 SI) against the
weakly magnetic host rock (~10⁻⁴ SI).

Total Magnetic Intensity (TMI) observations were extracted from a merged aeromagnetic raster and reduced to
20 representative stations via random subsampling. Following the established preprocessing pipeline, the IGRF
intensity (~47,500 nT) was subtracted from raw TMI values to obtain anomaly-level measurements consistent
with the forward model output.

## 5.2 Direct Position Priors for Lens Geometry

A key methodological distinction for the Soricom case is the choice to parameterize uncertainty directly in
terms of **surface point Z positions** rather than the orientation or dip priors used in the Tharsis
study. This design choice reflects the geological nature of the problem at deposit scale: for a thin,
high-contrast lens body hosted within a faulted sequence, the primary unknowns are:

- **Lens position**: Where exactly does the chromite lens sit within the host rock?
- **Lens thickness**: How thick is the lens, and does it pinch or swell laterally?
- **Host rock interface depth**: At what depth does the upper boundary of the ultramafic host occur?

Orientation-based priors (e.g., Normal distributions on dip angles) are well-suited for regional studies
where interface attitudes are the dominant source of geometric uncertainty. However, for a mineralized
lens where contacts are sub-parallel, the critical question is *where* the body is located — not *at what
angle*. By directly perturbing the Z coordinates of individual surface points, the inversion can
independently explore translations of the entire lens body as well as local variations in lens thickness
along its strike. This approach also preserves full differentiability for gradient-based inference via NUTS,
since Z-coordinate perturbations propagate naturally through the implicit interpolation.

The prior structure is organized as a two-tier Z-shift model with 10 parameters:

| Parameter | Geological target | Prior | Interpretation |
|-----------|-------------------|-------|----------------|
| `z_shifts[0]` | All host rock points (12 points) | Normal(0, 15 m) | Bulk vertical shift of the host rock upper contact |
| `z_shifts[1]…z_shifts[9]` | Individual chromite lens points (9 points) | Normal(0, 15 m) | Independent vertical perturbation per lens surface point |

The use of **independent priors** for each lens surface point allows the ensemble to express a
richer range of lens geometries than a single-thickness or single-position parameter would permit:
uniform thickness changes (all lens points shift together in the posterior), lateral thickness
variations (some points shift more than others), and coherent translations (lens and host shift in the
same direction). This parameterization naturally captures the geological expectation that the lens may
pinch out, swell, or shift vertically within the host, while remaining compact enough (10 parameters)
to be efficiently explored by MCMC.

In combination with the position priors, LogNormal priors are placed on the susceptibility values for all
four units, with the chromite lens centred at a high susceptibility (~0.5 SI) reflecting its strong
magnetic signature relative to the weakly magnetic host.

## 5.3 Inference and Results

Inference was performed using the No-U-Turn Sampler with 200 warmup and 200 sampling steps. A custom
`set_interp_input_fn` applies sampled Z shifts to the correct surface-point index ranges
(indices 1–12 for host rock, 13–21 for the lens) at each MCMC iteration, recomputing the implicit model
and TMI forward response. An equivalent `update_model_fn` is used during post-inference visualization and
probability field generation to consistently apply posterior samples to the geological model.

**Outcomes:** Posterior analysis of the Z-shift parameters revealed significant learning from the magnetic
data, with posterior distributions substantially narrower than the priors and shifted away from zero.
The magnetic observations provided meaningful constraints on both the bulk host rock depth and the
individual lens point positions. Where multiple lens points showed consistent posterior shifts, the data
favoured a particular lens thickness and vertical position; where posterior uncertainty remained elevated,
the magnetic response was less diagnostic of local lens geometry — an interpretable outcome that directly
informs where additional data (e.g., drilling) would be most valuable.

Posterior predictive checks demonstrated a good match between observed TMI anomalies and the forward
model ensemble mean, confirming that the Z-position parameterization, combined with susceptibility
inference, can explain the magnetic observations without introducing spurious structural complexity.

Probability density fields and information entropy maps were computed for both prior and posterior
ensembles, revealing a marked reduction in lithological ambiguity within the lens–host rock transition
zone after conditioning on magnetic data. The entropy reduction was most pronounced directly beneath
the lens and in the immediate vicinity of the observation stations, where the TMI gradient is most
sensitive to lens thickness and depth. The probability fields provide a spatially explicit map of
where the chromite lens is most likely to be found and how its thickness varies across the prospect.

Key geological takeaway: for deposit-scale targets characterized by thin, high-contrast bodies,
direct position priors on surface-point Z coordinates offer a more natural and interpretable
uncertainty parameterization than orientation or dip priors. This approach explicitly quantifies the
posterior uncertainty in lens position and thickness — directly translatable into actionable information
for resource estimation, drill targeting, and risk assessment.

Current limitations: the inversion uses a fixed (non-hierarchical) noise model (σ = 150 nT);
the 20-station subset may not fully capture spatial TMI variability across the prospect; and lens
geometry uncertainty remains conditional on the assumed number and distribution of surface points.
Relaxing the point count, exploring transdimensional lens representations, or incorporating remanent
magnetization are natural extensions for future work.

*[Figure suggestion: Prior vs posterior Z-shift parameter marginal distributions, showing narrowing and shifting of individual lens point positions.]*
*[Figure suggestion: Posterior probability density fields and information entropy cross-sections through the chromite lens, highlighting uncertainty reduction after magnetic conditioning.]*
*[Figure suggestion: TMI observation vs forward prediction correlation plot with posterior ensemble spread and 1:1 reference line.]*
