# crywolf 1.0.0

First release, accompanying the simulation study of Provete et al.
(fpcalib research compendium, v0.1.0, doi:10.5281/zenodo.23379100).

## Models

* `fp_occu()` fits the multiple-methods false-positive occupancy model of
  Miller et al. (2011) in Stan, with any number of occupancy, detection, and
  false-positive covariates, an optional recorder random intercept, and
  configurable priors (`fp_priors()`). `fp = NULL` gives the standard
  occupancy model.
* Observation-confirmation design of Chambert et al. (2015): argument
  `confirmed` of `fp_data()` and `frac_confirmed` of `simulate_fp()`.
* Methods `print`, `summary`, `coef`, `predict`, `loo`, and `as_draws_df`.
* `fp_ppc()`: posterior predictive check with a Bayesian p-value.

## Simulation and calibration

* `simulate_fp()` generates homogeneous, unstructured, and structured
  false-positive regimes; with the same seed it returns the same datasets
  as the simulator of the fpcalib compendium.
* `fp_calibrate()` and `fp_verification_curve()` report the coverage and
  bias of credible intervals for a design, with a `plot()` method and a
  progress bar that shows the time remaining (needs 'progressr'). They fit
  the recorder random intercept only when the simulated design has one,
  which avoids divergent transitions.

## Data and documentation

* `curiango` dataset: acoustic monitoring of the Common Pauraque in the
  Pantanal (Piatti et al. 2026).
* Articles: introduction, simulation, calibration, a case study with the
  `curiango` data, coming from `unmarked::occuFP()`, and the Stan model.
