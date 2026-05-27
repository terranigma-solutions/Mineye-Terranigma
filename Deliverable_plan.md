# Horizon Europe Deliverable - Structural Analysis and Restructuring Plan

## 1. Analysis of the Current Draft (Repetitiveness and Flow)

The current draft provides a solid foundation of work but suffers from significant structural repetition, which dilutes the impact of the scientific findings. The repetitiveness primarily stems from mixing general methodology with case-study-specific applications using an identical narrative template.

### Key Areas Identified:
1. **Methodology vs. Application Duplication:** 
   - The document explains the theoretical underpinnings of Bayesian Gravity, Magnetic, and Joint inversions (e.g., sections labeled "Bayesian Gravity Inversion", "Bayesian Magnetic Inversion"). 
   - Later, in the application sections for "Tharsis AOI 1", the exact same theoretical concepts (how the likelihood is defined, what the prior is) are re-explained to the reader. 
2. **Cookie-Cutter Paragraphs:** 
   - The application sections (Gravity, Magnetic, EnMap) follow a nearly identical, formulaic template. For instance, the MCMC sampling description is repeated verbatim four times:
     > *"Inference was carried out using a sampling-based approach, generating an ensemble of model realizations (~200-500 samples). Each iteration involves updating model parameters, recomputing the implicit geological model..."*
   - Sentences stating that the current method provides a *"basis for integrated interpretation within joint inversion workflows"* are repeated at the end of multiple sections as a concluding thought.
3. **Misplaced Disclaimers:**
   - The caveat regarding BaySeg (*"It should be noted that the implemented segmentation serves primarily as a demonstration..."*) is stated in both the single EnMap inversion and the Joint Inversion sections. 

### Impact:
- **Reader Fatigue:** Reviewers are forced to read the mechanics of Bayesian sampling multiple times.
- **Loss of Focus:** The actual geological findings, uncertainty reductions, and the specific impact on Tharsis AOI 1 are overshadowed by redundant methodological explanations.

---

## 2. Detailed Restructuring Plan

To elevate the document to the standards expected of a Horizon Europe deliverable, we need a strict separation of **Framework/Methodology** and **Application Case Studies**.

### Proposed Revised Outline

#### **Section 1: Probabilistic Framework Modernization**
*Focus: Software engineering, computational improvements, and the new API.*
* **1.1 Migration to PyTorch & Probabilistic API** (Merge "General Platform modernization: GemPy" and "1.1.1 Probabilistic API and Data Structures"). Emphasize the transition from Theano to PyTorch.
* **1.2 Performance & Numerical Stability** (Current 1.1.2). Expand slightly on how vectorized operations and the PyTorch backend prepare the system for GPU acceleration.

#### **Section 2: Forward Modeling Operators**
*Focus: The physics connecting geology to the likelihood functions.*
* **2.1 Gravity Formulation.** 
* **2.2 Magnetics (TMI) Formulation.** 
* **2.3 Categorical Likelihoods & Ordinal Probabilities (EO-derived data).** Move the general BaySeg explanation and the associated disclaimer here. Explain *how* soft probabilities are mapped without tying it strictly to Tharsis.

#### **Section 3: Bayesian Inference Strategy**
*Focus: The mathematical and algorithmic engine (Consolidating the repetitive text).*
* Create a single, robust section explaining the inference workflow. 
* **Key Additions based on codebase:** Mention the specific algorithmic stack. Instead of a generic "sampling-based approach", state that the framework utilizes the **No-U-Turn Sampler (NUTS)** via **Pyro**, running Markov Chain Monte Carlo (MCMC). 
* Move the repetitive paragraph here once: *"Inference is carried out using the NUTS MCMC algorithm... Each iteration involves updating model parameters, recomputing the implicit geological model, and evaluating the forward response..."*
* Discuss the concept of Joint Inversion here mathematically (how different likelihoods are combined).

#### **Section 4: Application Case Study - Tharsis AOI 1**
*Focus: Purely on the data, the priors, and the geological insights.*
* **4.1 Area Setup:** Introduce the initial GemPy model for Tharsis AOI 1 and the target parameters (e.g., interface positions, density contrasts).
* **4.2 Gravity Constraints:** Discuss the specific gravity data used. What did the posterior show? *(E.g., "Posterior analysis showed a reduction of uncertainty in... while poorly constrained areas...")*
* **4.3 Magnetic (TMI) Constraints:** Discuss the TMI data used. Contrast the findings with the gravity results.
* **4.4 EnMap EO-derived Constraints:** Discuss the hyperspectral data results. Note how it constrains near-surface geology but decreases in influence with depth.
* **4.5 Joint Inversion Outcomes:** Present the culmination of all constraints. How did combining Gravity and EnMap reduce the overall ambiguity compared to the single inversions?

---

## 3. Actionable Tasks for the Drafting Phase
1. **Strip out boilerplate:** Remove all mentions of "sampling-based approach" and "iteration updates" from the Tharsis sections.
2. **Move the BaySeg disclaimer:** Place it once in Section 2.3.
3. **Elevate technical language:** Use terms like "NUTS", "MCMC", "PyTorch", and "Pyro" to highlight the modern stack. 

---

## 4. Questions for the Author (To Finalize the Draft)

1. **Sampling Details & Consistency:** 
   The draft mentions generating "~200-500 samples". However, the '04_gravity_inversion.py' code example shows '1000' samples and '300' warmup steps using the 'NUTS' algorithm in 'Pyro'. 
   *Question: Should we update the text to explicitly mention the NUTS algorithm and reflect the higher sample counts (1000) for consistency and technical rigor?*

No

2. **Figure References and Captions:** 
   You have placeholders like '[Figure: initial model vs constrained model?]' and '[FIGURE: Lith Segmentation, derived classes]'. 
   *Question: Do you have these figures ready? Should we add specific technical captions referencing the specific outputs (e.g., "entropy maps" or "variance fields" as generated by ArviZ in the example script)?*

No, but you could add a caption of what you think it could look like

3. **Placement of the BaySeg Disclaimer:** 
   You added a disclaimer twice that *"BaySeg serves primarily as a demonstration... not limited to BaySeg"*. 
   *Question: Are you comfortable moving this disclaimer entirely to the general methodology section (Categorical Likelihoods) so we don't repeat it in the Tharsis application sections?*

yes

4. **Target Audience Balance:** 
   *Question: Is the primary audience for this deliverable technical reviewers looking for machine learning/software engineering details (where we should expand more heavily on the PyTorch/Pyro backend advantages), or geoscientists looking for impact (where we should expand more on the specific Tharsis AOI 1 geological findings and risk reduction)?*

probably geosicentists

