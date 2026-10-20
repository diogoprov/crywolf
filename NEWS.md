# crywolf 1.3.0

* Faster likelihood, with the same posterior: the unoccupied branch of each
  site is computed from its totals of automated occasions and detections,
  the occupied branch uses one `binomial_logit_glm_lpmf()` call per site
  when detection has covariates, and with `random_fp = FALSE` the site
  effects enter as data rather than as parameters fixed at zero. On the
  curiango data each gradient evaluation is 4 times faster than in version
  1.2.0 with one thread, and a short run of the sampler with four chains
  took 2.3 times less per leapfrog step with one thread and 2.6 times less
  with four. The log density, its gradient and `log_lik` are unchanged up
  to floating-point rounding. Executables compiled by earlier versions are
  rebuilt once.

# crywolf 1.2.0

* `fp_occu(random_fp = TRUE)` adds a site random intercept to the logit of
  the false-positive rate, for rates that vary among sites in ways the
  covariates do not capture; its prior scale is `fp_priors(fp_site_sd = )`.
  It is off by default: it slows the fits and mixes less well, and it does
  not replace covariates when false positives follow an environmental
  gradient (see `?fp_occu`). With the option off, the draws are identical
  to those of crywolf 1.1.0.
* `fp_calibrate_fit()` checks the calibration of candidate models on data
  simulated from a fitted model, which keep the design of the real data
  (sites, covariates, occasions, recorders, and verification).
* `simulate()` method for `fp_occu` fits, from the posterior means or from
  random posterior draws.
* The summary of a calibration prints every occupancy slope, not only `occ_x`.

# crywolf 1.1.0

* The likelihood is summed over sites with Stan's `reduce_sum()`, and the
  model is compiled with threading, so `fp_occu(threads_per_chain = )` can
  split each chain among several cores. On the curiango data (97,806
  automated occasions) a fit with four chains of four threads each took
  16.5 minutes instead of 70, and each gradient evaluation was 2.7 times faster than with the
  previous model in a short benchmark on the same computer. The log
  density and its gradient are unchanged up to floating-point rounding, and
  the automated occasions use `y * lp - n * log1p_exp(lp)`, an equivalent
  form of the binomial log-likelihood with half the transcendental
  functions. Executables compiled by earlier versions are rebuilt once.

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
