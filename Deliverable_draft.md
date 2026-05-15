General Platform modernization: GemPy

As part of Task 2.2, the probabilistic geomodelling framework of GemPy was comprehensively modernized to support scalable Bayesian inference and the integration of Earth Observation (EO)-derived constraints. This modernization was necessary to align the probabilistic modelling capabilities with the new GemPy v3 engine and to enable robust, extensible workflows for data fusion in geological modelling.

A complete refactoring of the probabilistic framework was carried out, including the migration from the legacy Theano backend to the modern tensor-based PyTorch framework. This transition significantly improved performance, flexibility, and long-term maintainability of the codebase. Deprecated dependencies were removed, and the overall structure of the package was reorganized to support future extensions and integration into larger processing environments .

A central outcome of this effort is the introduction of a unified data model, enabling consistent handling of heterogeneous data sources. This includes geological inputs such as interface points and orientations, EO-derived constraints, and geophysical observations. The unified representation forms the basis for joint probabilistic modelling and inversion workflows.

In parallel, a modular probabilistic API was developed to standardize the definition and execution of Bayesian workflows. This API enables users to define prior distributions and likelihood functions, execute inference procedures, and access posterior diagnostics, including analysis and visualization of uncertainty. The framework was further extended with new functionality for posterior evaluation and monitoring during inference.

Significant improvements were achieved in numerical stability, particularly in the probabilistic learning phase. This includes enhanced methods for scalar field-based segmentation of geological features such as faults and layers, which are critical for stable model convergence  . Additionally, automated test coverage was substantially expanded to ensure robustness and reproducibility of probabilistic workflows.

Finally, the modernized framework was prepared for integration into external systems and workflows, including EO processing pipelines and the LiquidEarth/IPOP platform. This ensures that the developed probabilistic methods can be deployed within the broader MINEYE ecosystem and used in end-to-end data fusion and modelling applications.



1.1.1 Probabilistic API and Data Structures

To enable efficient and user-friendly probabilistic geomodelling, a high-level API was designed to simplify the setup, execution, and analysis of Bayesian workflows. The API abstracts the underlying complexity of probabilistic inference, allowing users to define models, priors, and likelihoods in a consistent and structured manner while maintaining full flexibility for advanced use cases.

A key component of this development is the standardization of internal data representations. Core model elements such as scalar fields, geological surfaces, and observation datasets are now handled through unified data structures. This ensures consistency across workflows and facilitates the integration of heterogeneous data sources within a single probabilistic framework.

In addition, dedicated interfaces were implemented to support the incorporation of external data into the modelling process. This includes mechanisms for adding EO-derived constraints as probabilistic inputs, enabling their direct use within likelihood formulations. Furthermore, flexible connectors were developed to link external datasets such as raster data products and segmentation outputs, allowing seamless integration of results from upstream EO processing workflows.

Finally, interoperability with the broader Python scientific and machine learning ecosystem was significantly improved. The API is designed to integrate naturally with common libraries for numerical computing, data handling, and machine learning, facilitating the development of hybrid workflows that combine probabilistic modelling with modern data-driven approaches.



1.1.2 Numerical Stability and Performance Improvements

A systematic analysis of the previous probabilistic implementation identified several sources of numerical instability and performance limitations. Key issues included the use of non-vectorized operations and inconsistent parameter handling, both of which negatively affected convergence behaviour and computational efficiency.

To address these limitations, critical components of the framework were refactored and implemented using fully vectorized operations. This resulted in more stable and efficient execution of probabilistic inference workflows, particularly in the context of iterative Bayesian methods. In parallel, parameter handling was standardized across the codebase, reducing ambiguity and ensuring consistent behaviour during model evaluation and sampling.

These improvements led to a significantly enhanced convergence behaviour in Bayesian inference workflows. Numerical stability during the sampling process was increased, reducing the likelihood of divergence and improving the reliability of posterior estimates.

Memory usage was further optimized through improved data handling strategies, enabling more efficient storage and processing of intermediate results. This is particularly relevant for large-scale applications where memory constraints can become a limiting factor.

The migration to the PyTorch backend also prepares the framework for GPU acceleration, allowing computationally intensive operations to be executed on modern hardware architectures. This provides a clear pathway for further performance gains in future developments.

Overall, these enhancements substantially increase the robustness of the framework when applied to demanding use cases, which are central to the objectives of the MINEYE project 

Forward Modeling implementations

Forward modelling operators were developed to establish a direct link between geological models and geophysical observables, forming the basis for probabilistic inversion. These operators translate the implicit geological representation into physically meaningful property fields, enabling the simulation of measurable geophysical responses.

Within this work, forward models were implemented for gravity and magnetic data, both belonging to the class of potential field methods. In this context, potential fields are governed by scalar potentials and can be computed efficiently from volumetric property distributions without the need to solve partial differential equations. This makes them particularly well-suited for probabilistic workflows that require repeated forward evaluations.

The geological scalar field representation is coupled to physical property distributions, such as density and magnetic susceptibility, allowing consistent propagation of geological information into geophysical space. The resulting forward responses are integrated into the probabilistic framework as likelihood components, enabling quantitative comparison between simulated and observed geophysical data.

Other geophysical methods, such as electrical resistivity, electromagnetic, or seismic approaches, require the numerical solution of differential equations and are therefore significantly more computationally demanding. While these methods are in principle compatible with probabilistic inversion, their computational cost currently limits their applicability in large-scale Bayesian workflows.

Overall, the implemented forward models provide a robust and efficient foundation for inversion and joint inversion workflows, where multiple data sources can be combined to improve model accuracy and reduce uncertainty.



Forward Gravity Formulation

A gravity forward modelling approach was implemented based on density contrasts between geological units represented within an implicit geological model. The implicit formulation allows continuous representation of subsurface structures, which are subsequently translated into volumetric density fields for geophysical computation.

Lithological model outputs are converted into density distributions, including support for spatially varying densities within individual units. This enables the representation of intra-unit heterogeneity and improves the realism of simulated gravity responses.

The gravity signal is computed using a discretized volumetric representation of the model domain. This approach provides the flexibility to handle complex geometries while maintaining computational efficiency. The implementation is designed to operate on regional-scale datasets and supports irregular observation grids, reflecting typical acquisition conditions in gravity surveys.

The forward model is tightly integrated into the probabilistic framework through the likelihood function, where simulated gravity anomalies are compared to observed data. This enables direct incorporation of geophysical information into Bayesian inference workflows.

The formulation is optimized for repeated forward evaluations, making it suitable for probabilistic inversion schemes that require large numbers of model realizations during sampling.



Forward Magnetics Formulation

A forward modelling approach for magnetic data was implemented based on the distribution of magnetic susceptibility within geological units. Geological model outputs are mapped to susceptibility values, enabling the construction of volumetric property models from the implicit geological representation.

Magnetic anomalies are computed from this volumetric model, accounting primarily for induced magnetization driven by the ambient magnetic field. Where applicable, the formulation can be extended to include remanent magnetization, allowing for more realistic representation of magnetic behaviour in specific geological settings.

The forward model is integrated into the probabilistic framework as a likelihood term, enabling comparison between simulated and observed magnetic anomalies. This allows magnetic data to directly constrain geological model realizations within Bayesian inference workflows.

The implementation supports multi-scale magnetic datasets and accommodates spatially heterogeneous susceptibility distributions. Similar to the gravity formulation, it is designed for efficient repeated evaluation, ensuring suitability for computationally intensive probabilistic inversion and sampling procedures.



Bayesian Formulations

Bayesian formulations were developed to enable the consistent integration of geological, geophysical, and EO-derived information within a unified probabilistic framework. This approach allows explicit representation and propagation of uncertainty, while incorporating multiple data sources as constraints through likelihood functions. The implemented formulations cover uncertainty estimation, individual inversion of gravity and magnetic data, integration of categorical EO-derived information, and joint inversion strategies. Together, these components provide a flexible foundation for data-driven geological modelling and uncertainty-aware interpretation within the MINEYE workflow.

Uncertainty Estimation in Geological Models

Geological structures are represented as probabilistic entities within the modelling framework, allowing uncertainties in data and model assumptions to be explicitly captured. Instead of a single deterministic model, posterior distributions are sampled to generate ensembles of plausible geological realizations.

Uncertainty is quantified using multiple complementary measures, including variance fields, entropy maps, and ensemble-based statistics. These metrics provide spatially resolved insight into model confidence and variability across the domain.

The analysis enables the identification of well-constrained regions, where data strongly supports the model, as well as zones of high uncertainty, which may require additional data acquisition or alternative modelling assumptions. Visualization tools were developed to support interpretation, including monitoring of posterior evolution during sampling and assessment of uncertainty propagation through the model.

This probabilistic representation of geological models provides a foundation for risk-aware decision-making in exploration contexts, where uncertainty plays a critical role in planning and evaluation.



Bayesian Gravity Inversion

The gravity forward model was integrated into a Bayesian inference framework to enable inversion of geological models based on observed gravity data. A likelihood function is defined from the misfit between observed gravity anomalies and the simulated response derived from the geological model.

Within this framework, posterior distributions of geological parameters are estimated, conditioned on the available gravity observations. This allows key model components, including layer geometries and density distributions, to be updated during inference.

Measurement noise and data uncertainty are explicitly accounted for in the likelihood formulation, ensuring robust handling of real-world gravity datasets. The inversion proceeds iteratively, refining the geological model through repeated forward evaluations and comparison with observed data, leading to progressively improved model realizations.



Bayesian Magnetic Inversion

Magnetic forward modelling was incorporated into the Bayesian inference framework to enable inversion of geological models using magnetic observations. The likelihood function is defined based on the misfit between observed magnetic anomalies and the simulated response derived from the susceptibility model.

The inversion process estimates posterior distributions of both susceptibility fields and geological structures, allowing magnetic data to directly constrain subsurface interpretations. Uncertainty in the observations is explicitly handled through noise modelling in the likelihood, while the inherent ambiguity of magnetic data is addressed through the probabilistic formulation.

The framework supports the integration of multiple magnetic datasets, enabling joint use of data with different resolutions and acquisition characteristics. Compared to deterministic approaches, this leads to improved constraints on subsurface structures and provides a more comprehensive characterization of uncertainty in the resulting models.



Categorical Likelihood and Ordinal Probabilities

Likelihood functions were developed to incorporate categorical geological information into the probabilistic modelling framework. This enables the direct use of EO-derived classifications, such as lithological classes obtained from segmentation workflows, as constraints within Bayesian inference.

The approach distinguishes between categorical likelihoods for discrete class assignments and ordinal probabilities for ordered geological relationships, such as stratigraphic sequences. EO segmentation outputs are mapped onto the geological model domain and translated into probabilistic constraints that inform the likelihood evaluation.

To account for uncertainty in classification results, soft probabilistic representations are used instead of hard labels. This allows class membership probabilities to be propagated through the inference process, ensuring that uncertainty in EO-derived inputs is consistently reflected in the resulting geological models.

The formulation enables the integration of a wide range of data sources, including remote sensing products and interpreted geological maps. In particular, uncertainty estimates derived from Bayesian segmentation methods such as BaySeg can be directly incorporated, allowing for coherent propagation of uncertainty from EO data processing into geological modelling.

Bayesian Joint inversion

A framework for Bayesian joint inversion was developed to enable the simultaneous integration of multiple data types within a unified probabilistic formulation. This includes the combination of gravity data, magnetic data, and EO-derived constraints within a single inversion workflow.

A joint likelihood function is constructed by integrating the individual contributions from each observation type. This allows all available data sources to consistently inform the inference process, ensuring that geological model updates are driven by the combined evidence.

The use of complementary datasets improves model constraints and reduces ambiguity compared to single-data inversions. While individual data types may be non-unique or ambiguous, their joint interpretation provides a more robust and consistent characterization of subsurface structures.

The framework enables the consistent update of geological parameters across all data sources during inference, ensuring coherence between geological, geophysical, and EO-derived information. It provides the foundation for multi-physics inversion approaches and supports integrated workflows that combine EO data processing with geophysical modelling and inversion. 

Bayesian Gravity Inversion

A Bayesian gravity inversion was performed for Tharsis AOI 1 to constrain the geological model using available gravity observations. The inversion builds on the initial GemPy model and integrates the gravity forward formulation described in Section 1.3 as a likelihood component within the probabilistic framework.

The inversion targets both geometrical parameters (e.g. interface positions) and physical properties (density contrasts between units). Prior distributions were defined based on the initial model and available petrophysical information, while the likelihood was formulated from the misfit between observed gravity data and simulated responses.

Inference was carried out using a sampling-based approach, generating an ensemble of model realizations (~200-500 samples). Each iteration involves updating model parameters, recomputing the implicit geological model, and evaluating the forward gravity response. The resulting posterior distribution captures the range of geological configurations consistent with the observed data.

The inversion demonstrates that gravity data provides meaningful constraints on large-scale density contrasts and interface geometries, despite the limited subsurface information available for this area. At the same time, ambiguity remains in regions where different model configurations produce similar gravity responses, reflecting the non-uniqueness of potential field inversion.

Posterior analysis shows a reduction of uncertainty in regions well covered by gravity data, while poorly constrained areas remain sensitive to prior assumptions. Overall, the results highlight the value of integrating gravity data within a probabilistic framework and provide a basis for further refinement through joint inversion with additional data sources.

[Figure: initial model vs constrained model?]
[Figure: constrained model, information entropy]

Bayesian Magnetic Inversion (TMI)

A Bayesian magnetic inversion was carried out for Tharsis AOI 1 to further constrain the geological model using magnetic observations. The inversion builds on the initial GemPy model and incorporates the magnetic forward formulation as a likelihood component within the probabilistic framework.

The inversion is based on Total Magnetic Intensity (TMI) data, which represents the measured magnetic field anomalies. The likelihood is defined through the misfit between observed TMI and the simulated magnetic response derived from the susceptibility model. The formulation primarily considers induced magnetization, consistent with the available data and modelling assumptions.

The inversion focuses on the estimation of susceptibility distributions and their relationship to geological structures. Prior information is derived from the initial model and available petrophysical constraints. Inference was performed using a sampling-based approach, generating an ensemble of model realizations (~200-500 samples). Each iteration updates model parameters, recomputes the implicit geological model, and evaluates the corresponding magnetic response.

The results indicate that TMI data provides additional constraints on subsurface structures, particularly in areas with significant susceptibility contrasts. Compared to gravity, magnetic data can offer higher sensitivity to lithological variations, although ambiguity remains due to non-uniqueness and trade-offs between geometry and property distributions.

Posterior analysis shows improved constraint in areas with strong magnetic signal, while uncertainty persists in regions with weak or ambiguous responses. The inversion complements the gravity-based results and provides a basis for integrated interpretation within joint inversion workflows.

Bayesian EnMap Inversion: Categorical Likelihood and Ordinal Probabilities

A Bayesian inversion incorporating EO-derived constraints was performed for Tharsis AOI 1 using hyperspectral data from EnMAP. The objective is to integrate surface lithological information into the probabilistic geomodelling framework through categorical likelihoods.

Categorical constraints were derived from hyperspectral imagery using a Bayesian segmentation approach based on BaySeg. The workflow includes spectral preprocessing (e.g. band selection and smoothing), masking of non-relevant pixels, and clustering of spectral signatures into lithological classes. The segmentation produces probabilistic class assignments rather than hard labels, allowing uncertainty in classification to be explicitly represented. These class probabilities are mapped onto the geological model domain and used to define a categorical likelihood, linking EO-derived surface information to subsurface model realizations.

It should be noted that the implemented segmentation serves primarily as a demonstration of the workflow. The approach is not limited to BaySeg and can incorporate more advanced or site-specific lithological classification methods developed within MINEYE, including alternative machine learning or expert-driven interpretations.



The inversion is performed by sampling model realizations (~200-500 samples), where each iteration updates geological parameters, recomputes the implicit model, and evaluates the agreement between model-derived lithologies and EO-based class probabilities. This allows EO-derived information to directly influence the posterior distribution of geological structures.

The results demonstrate that EO constraints provide valuable information on surface geology, improving model consistency near the surface and helping to reduce ambiguity in areas with limited subsurface data. At the same time, the influence of EO data decreases with depth, highlighting the importance of combining EO constraints with geophysical data in joint inversion workflows.

Overall, the approach illustrates how probabilistic EO-derived classifications can be integrated into geological modelling, enabling consistent propagation of uncertainty from remote sensing data into subsurface interpretations.



Bayesian Joint Inversion: Gravity and EnMap

A Bayesian inversion incorporating EO-derived constraints was performed for Tharsis AOI 1 using hyperspectral data from EnMAP. The objective is to integrate surface lithological information into the probabilistic geomodelling framework through categorical likelihoods.

Categorical constraints were derived from hyperspectral imagery using a Bayesian segmentation workflow based on BaySeg. The processing includes extraction of relevant spectral bands from EnMAP data, filtering of noisy wavelength regions, and application of smoothing and normalization techniques. Additional masking steps are applied to remove vegetation, clouds, and water using quality layers, ensuring that only geologically meaningful pixels are considered. The resulting spectral dataset is segmented into lithological classes, producing probabilistic class memberships for each pixel. (Add: Sampling method?)

These probabilistic classifications are spatially aligned with the geological model and mapped onto the model domain. Instead of using hard class assignments, soft probabilities are used to define a categorical likelihood, allowing uncertainty in the EO-derived segmentation to be propagated into the inversion. The likelihood evaluates the agreement between model-predicted lithologies at the surface and the EO-derived class probabilities.



It should be noted that the BaySeg-based segmentation is used here as a demonstration of the workflow. The approach is designed to be flexible and can incorporate more advanced or site-specific lithological classifications developed within MINEYE, including alternative machine learning methods or expert-derived interpretations.

The inversion is performed using a sampling-based approach (~200-500 realizations), where each iteration updates geological parameters, recomputes the implicit model, and evaluates consistency with the EO-derived constraints. The results show that EO data provides strong constraints on near-surface geology, improving model consistency in areas with limited subsurface information.

At the same time, the influence of EO constraints decreases with depth, highlighting the importance of combining EO-derived information with geophysical data in joint inversion workflows. Overall, this approach demonstrates how probabilistic EO-derived classifications can be effectively integrated into geological modelling, enabling consistent uncertainty propagation from remote sensing data to subsurface interpretations.



[FIGURE: Lith Segmentation, derived classes] 

