# crywolf 1.4.0

## Plots and site-level estimates

* `plot()` for fits draws the posterior intervals of the coefficients, one
  panel per process, or sampler diagnostics (`type = "trace"`, `"rank"`,
  `"pairs"`, via 'bayesplot', with divergent transitions marked).
* `fp_effects()` gives the occupancy probability or the false-positive rate
  along one covariate, with the others at their means, and a `transform`
  argument for covariates standardized before fitting; `plot()` draws it.
* `fp_sites()` returns, for every site, the occupancy probability, the
  posterior probability that it is occupied given its data, and its
  false-positive rate; `plot()` ranks the sites or, with coordinates, maps
  them.
* `plot()` for `fp_ppc()` results: replicated against observed discrepancy,
  or the distribution of detections per site.
* `fp_prior_predictive()` simulates data from the priors with the design of
  the real data and compares simple summaries with the observed ones.
* `plot()` for calibrations can show verification effort as the expected
  number of verified sites or occasions, and mark a target coverage.

## Reading classifier output

* `fp_read_birdnet()` reads BirdNET selection tables (Raven format).
* `fp_read_validation()` reads the clip verdicts of the
  birdnet-validation-sampler.
* `fp_aggregate()` turns detections into automated occasions at a
  confidence threshold, and verdicts into counts for the
  observation-confirmation design, warning when a checked occasion has
  detections that were not reviewed.
* `fp_data()` warns about sites without automated occasions and about
  automated data that are all detections or all non-detections.

## Other changes

* Fits keep every draw in memory, so a fit saved with `saveRDS()` works in
  later sessions, after CmdStan's temporary files are gone.

# crywolf 1.3.0

## Faster likelihood, same posterior

The Stan program evaluates the likelihood up to three times faster, with
the log density, its gradient and `log_lik` unchanged up to floating-point
rounding (checked in nine designs; scripts and results in `dev/profile/`
of the repository). Executables compiled by earlier versions are rebuilt
once, and CmdStan 2.34 or later is now required.

* The unoccupied branch of each site is computed from its totals of
  automated occasions and detections, because the false-positive rate is
  constant within a site.
* When `det` has covariates, the occupied branch uses one
  `binomial_logit_glm_lpmf()` call per site, which evaluates the linear
  predictor, the log-likelihood and its gradient in one pass instead of
  creating an autodiff variable per occasion.
* With `random_fp = FALSE` the site effects on the false-positive rate enter
  as data rather than as parameters fixed at zero.

The gain grows with the number of occasion-level rows per site. Time per
gradient evaluation with one thread, same parameter values:

| Data | Rows | 1.2.0 | 1.3.0 | Gain |
|---|---|---|---|---|
| 300 sites x 20 occasions, `det = ~ 1` (rows collapsed) | 300 | 0.064 ms | 0.056 ms | 1.1x |
| 300 sites x 20 occasions, `det = ~ bio + aci`, recorder effect | 6,000 | 0.58 ms | 0.23 ms | 2.5x |
| `curiango`, `det = ~ bio + aci`, recorder effect | 97,806 | 8.5 ms | 2.6 ms | 3.3x |

In a short run of the sampler on the `curiango` data (four chains), each
leapfrog step took 2.3 times less with one thread and 2.6 times less with
four threads per chain.

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
