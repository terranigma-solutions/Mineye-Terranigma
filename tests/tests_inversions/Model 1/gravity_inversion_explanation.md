### Planning structure of the explanation
- Multiple steps on the construction of a probabilistic model:
  - Set up initial model: this is the base model (need more explanation): model_setup.py/setup_geomodel
  - Constructing the Probabilistic Model
      1. Define priors: Follow explanation
      2. Define deterministic parameters: Follow explanation 
      3. Define likelihood functions: follow explanation with example from generate_multigravity_likelihood_hierarchical_per_station
      4. Define set_prior: Basically the functions used to pass the priors to gempy (probabilistic_model.py)


- Prior predictive model:
    - test_run_predictive
    - This gives us the range of the forward observations:
      - The interesting thing about this model is since we are simulating 20 observations in every iteration,
       we end up with models that can explain some of the observations but not others. In this simplify model 
       that is obvious when we plot the priors in test_run_predictive_analysis
        1) This tells us two things, the models should be more complex (either by perturbing more of the gempy
         parameters or adding things like variability of density within a layer)
        2) Adding more data (like magnetics) will not help us too much.
        3) But we will be able to detect outliers that our model cannot explain buy running test_run_predictive_analysis on 
         the posterior
         

- Inference: _pyro_runner.py

- Explaining posterior:
    - Since there are no model that fits all the observations, the posterior is concentrating
     around the maximum number of observations. This leave the outliers (that may be more interesting)
     without being explained by this simplified model.
    - In a sense it is almost a regression where the "polynomial" is a complex gravity fw