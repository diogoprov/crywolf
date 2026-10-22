# Changelog

## crywolf (development version)

### Bug fixes

- [`print()`](https://rdrr.io/r/base/print.html) of an `fp_data` object
  counted the sites known to be occupied (from verified detections and
  confirmed true observations together) on the line of verified sites,
  so data with only confirmations showed “0 verified sites … 24 with a
  confirmed detection”. Verified sites, confirmed occasions, and sites
  known to be occupied now have their own lines.
- With `parallel = TRUE`, the workers of ‘furrr’ use the L’Ecuyer-CMRG
  generator, so `simulate_fp(seed = s)` inside
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md),
  [`fp_calibrate_fit()`](https://provetelab.org/crywolf/reference/fp_calibrate_fit.md),
  and [`fp_sbc()`](https://provetelab.org/crywolf/reference/fp_sbc.md)
  simulated other datasets than a sequential run with the same seeds.
  Every function that takes a `seed` now sets R’s default generators
  along with it, so the two modes give identical results. In a session
  with the default generators nothing changes:
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  still reproduces the datasets of the fpcalib compendium exactly.
- Two fits made within the same minute after the same
  [`set.seed()`](https://rdrr.io/r/base/Random.html) could share
  CmdStan’s output file names, which ‘cmdstanr’ draws from R’s random
  numbers; discarding the first fit then deleted the files of the
  second, and reading its draws failed with “File does not exist” or
  “Supplied CSV file is corrupt”.
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) now
  gives every fit its own file names, without using R’s random numbers.

### Tests

- The precomputed results of the vignettes carry the MD5 of the Stan
  program that made them, and a vignette shows a note when its results
  come from another version of the program. `dev/precompute-all.R`
  refreshes them all; `dev/restamp-precomputed.R` marks them current
  without recomputing, only after the posterior regression test passes.
  The MD5 is of the whole file, so even an edit to a comment in the Stan
  program shows the note until the results are restamped.
- A posterior regression test checks, in four designs, that the
  log-likelihood at saved parameter values is unchanged and that
  posterior means stay within Monte Carlo error of saved references
  (`tests/testthat/test-regression.R`; references made by
  `dev/make-regression-refs.R`).

### New articles

- “From BirdNET output to an occupancy model”: BirdNET selection tables
  and ‘sheepdog’ verdicts to the data of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).
- “Choosing a classifier threshold”:
  [`fp_threshold_sweep()`](https://provetelab.org/crywolf/reference/fp_threshold_sweep.md)
  on the `curiango` data.
- “Reporting a crywolf analysis”: the items of the Bayesian Analysis
  Reporting Guidelines, with a template for methods and results.
- “Detections that come in runs”: the Markov term on the `curiango`
  data.
- “False-positive rates that vary among sites and occasions”:
  `random_fp` and occasion-level covariates on simulated data.
- “Planning a study”: decisions to make before recording starts.
- “Performance”: computing time over the versions, threads, and what did
  not help.

Results that need long fits are precomputed by
`dev/precompute-markov.R`, `dev/precompute-threshold.R`, and
`dev/precompute-fp-variation.R`.

### Installation

- ‘sheepdog’, suggested for reading clip reviews, is now installed from
  <https://diogoprov.r-universe.dev> (listed in
  `Additional_repositories`) instead of from GitHub through `Remotes`.
  The README says how.

## crywolf 1.6.0

### Speed

- `fp_occu(threads_per_chain = "auto")` uses the physical cores divided
  among the parallel chains. On the `curiango` data, 150 iterations of
  one chain took 55 s with one thread, 16 s with four, and 10 s with
  eight. The default stays at one thread, because with more the order of
  the partial sums varies from run to run and a fit is no longer
  reproduced exactly by its `seed` (the posterior is the same). A fit
  whose chains took over a minute on one core now says so and suggests
  `"auto"`.
- [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  documents that two chains per fit give an effective sample size of
  about 2000 and, with the chains of a fit run one after another as in a
  parallel calibration, half the cost of four.
- Measured and not adopted (`dev/profile/rejected-2026-10-21.md`): a
  dense mass matrix, sampling the false-positive intercept with the gap,
  a centred recorder effect, and `reduce_sum_static()`. Considered and
  not pursued, from the same measurements: initialization by Pathfinder
  (warm-up is 40% of the time) and a QR decomposition (no posterior
  correlation among covariate effects above 0.5). The confirmation
  design keeps its current likelihood path, for lack of a large dataset
  to measure it on.
- With `random_fp = TRUE`, 3 of 20 fits of the unstructured regime had
  an R-hat above 1.01 with the default settings; a warm-up of 2000 (2
  of 20) or `adapt_delta = 0.99` (3 of 20) took longer for no clear
  gain, so the defaults stay.

### Data

- [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)
  reads `validation_by_clip.csv` from ‘sheepdog’ 2.0.0 and later
  (English column names and verdicts: `site`, `source_audio`,
  `clip_start`, `verdict = "positive"`, …) as well as the legacy
  `validacao_por_clipe.csv` of the sampler up to 1.2.0. The format is
  detected from the columns, and both give the same output, which now
  keeps `occasion_id` when the file has it.

## crywolf 1.5.0

### Model options

All three are off by default; with them off, the draws are identical to
those of crywolf 1.4.0.

- `fp_occu(ver = ~ ...)` gives the verified method covariates of its own
  (from the verified table or the sites), reported as `verified_*`.
- False-positive formulas can use covariates measured on each occasion
  (columns of the automated table), so the rate varies among the
  occasions of a site.
- `fp_occu(markov = "time")` adds the detection on the previous occasion
  of the same site as a covariate of true detection (`det_y_prev`), for
  species recorded in runs. The previous occasion must be one time step
  earlier: after a gap (a missing occasion) a new run starts.
  [`simulate()`](https://rdrr.io/r/stats/simulate.html) simulates the
  runs in order, and
  [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md)
  checks them one step ahead.

`dev/validate-model-options.R` checks the identity with 1.4.0 in three
designs, compares the R likelihood used by
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md),
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md),
and [`simulate()`](https://rdrr.io/r/stats/simulate.html) with the Stan
`log_lik` in six combinations of the options, and checks in simulations
that each option’s effect is recovered; its output is in
`dev/validate-model-options.txt`.

### Checking and comparing models

- [`fp_sbc()`](https://provetelab.org/crywolf/reference/fp_sbc.md):
  simulation-based calibration of the computation (Talts et al. 2018;
  Modrák et al. 2025), with the log-likelihood of the data as a test
  quantity besides the parameters, and
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).
- `loo(reloo = TRUE)` refits the model without each site whose Pareto k
  exceeds a threshold and uses the exact leave-one-out value for it. A
  site that is the only one with its recorder, or with a level of a
  factor, keeps its approximate value, with a warning.
- [`fp_loo_weights()`](https://provetelab.org/crywolf/reference/fp_loo_weights.md):
  stacking or pseudo-BMA weights of several fits.
- [`fp_gap_prior()`](https://provetelab.org/crywolf/reference/fp_gap_prior.md)
  turns the precision measured on a random sample of reviewed detections
  into a prior for the gap.
- `predict(process = "fp")` predicts for the automated occasions by
  default when the false-positive rate has occasion covariates.
- [`fp_threshold_sweep()`](https://provetelab.org/crywolf/reference/fp_threshold_sweep.md)
  fits the same model at several classifier thresholds and plots the
  estimates against the threshold.

## crywolf 1.4.0

### Plots and site-level estimates

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for fits
  draws the posterior intervals of the coefficients, one panel per
  process, or sampler diagnostics (`type = "trace"`, `"rank"`,
  `"pairs"`, via ‘bayesplot’, with divergent transitions marked).
- [`fp_effects()`](https://provetelab.org/crywolf/reference/fp_effects.md)
  gives the occupancy probability or the false-positive rate along one
  covariate, with the others at their means, and a `transform` argument
  for covariates standardized before fitting;
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.
- [`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)
  returns, for every site, the occupancy probability, the posterior
  probability that it is occupied given its data, and its false-positive
  rate; [`plot()`](https://rdrr.io/r/graphics/plot.default.html) ranks
  the sites or, with coordinates, maps them.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md)
  results: replicated against observed discrepancy, or the distribution
  of detections per site.
- [`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)
  simulates data from the priors with the design of the real data and
  compares simple summaries with the observed ones.
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) for
  calibrations can show verification effort as the expected number of
  verified sites or occasions, and mark a target coverage.

### Reading classifier output

- [`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md)
  reads BirdNET selection tables (Raven format).
- [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)
  reads the clip verdicts of the birdnet-validation-sampler.
- [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
  turns detections into automated occasions at a confidence threshold,
  and verdicts into counts for the observation-confirmation design,
  warning when a checked occasion has detections that were not reviewed.
- [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
  warns about sites without automated occasions and about automated data
  that are all detections or all non-detections.

### Other changes

- Fits keep every draw in memory, so a fit saved with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) works in later
  sessions, after CmdStan’s temporary files are gone.

## crywolf 1.3.0

### Faster likelihood, same posterior

The Stan program evaluates the likelihood up to three times faster, with
the log density, its gradient and `log_lik` unchanged up to
floating-point rounding (checked in nine designs; scripts and results in
`dev/profile/` of the repository). Executables compiled by earlier
versions are rebuilt once, and CmdStan 2.34 or later is now required.

- The unoccupied branch of each site is computed from its totals of
  automated occasions and detections, because the false-positive rate is
  constant within a site.
- When `det` has covariates, the occupied branch uses one
  `binomial_logit_glm_lpmf()` call per site, which evaluates the linear
  predictor, the log-likelihood and its gradient in one pass instead of
  creating an autodiff variable per occasion.
- With `random_fp = FALSE` the site effects on the false-positive rate
  enter as data rather than as parameters fixed at zero.

The gain grows with the number of occasion-level rows per site. Time per
gradient evaluation with one thread, same parameter values:

| Data | Rows | 1.2.0 | 1.3.0 | Gain |
|----|----|----|----|----|
| 300 sites x 20 occasions, `det = ~ 1` (rows collapsed) | 300 | 0.064 ms | 0.056 ms | 1.1x |
| 300 sites x 20 occasions, `det = ~ bio + aci`, recorder effect | 6,000 | 0.58 ms | 0.23 ms | 2.5x |
| `curiango`, `det = ~ bio + aci`, recorder effect | 97,806 | 8.5 ms | 2.6 ms | 3.3x |

In a short run of the sampler on the `curiango` data (four chains), each
leapfrog step took 2.3 times less with one thread and 2.6 times less
with four threads per chain.

## crywolf 1.2.0

- `fp_occu(random_fp = TRUE)` adds a site random intercept to the logit
  of the false-positive rate, for rates that vary among sites in ways
  the covariates do not capture; its prior scale is
  `fp_priors(fp_site_sd = )`. It is off by default: it slows the fits
  and mixes less well, and it does not replace covariates when false
  positives follow an environmental gradient (see
  [`?fp_occu`](https://provetelab.org/crywolf/reference/fp_occu.md)).
  With the option off, the draws are identical to those of crywolf
  1.1.0.
- [`fp_calibrate_fit()`](https://provetelab.org/crywolf/reference/fp_calibrate_fit.md)
  checks the calibration of candidate models on data simulated from a
  fitted model, which keep the design of the real data (sites,
  covariates, occasions, recorders, and verification).
- [`simulate()`](https://rdrr.io/r/stats/simulate.html) method for
  `fp_occu` fits, from the posterior means or from random posterior
  draws.
- The summary of a calibration prints every occupancy slope, not only
  `occ_x`.

## crywolf 1.1.0

- The likelihood is summed over sites with Stan’s `reduce_sum()`, and
  the model is compiled with threading, so
  `fp_occu(threads_per_chain = )` can split each chain among several
  cores. On the curiango data (97,806 automated occasions) a fit with
  four chains of four threads each took 16.5 minutes instead of 70, and
  each gradient evaluation was 2.7 times faster than with the previous
  model in a short benchmark on the same computer. The log density and
  its gradient are unchanged up to floating-point rounding, and the
  automated occasions use `y * lp - n * log1p_exp(lp)`, an equivalent
  form of the binomial log-likelihood with half the transcendental
  functions. Executables compiled by earlier versions are rebuilt once.

## crywolf 1.0.0

First release, accompanying the simulation study of Provete et al.
(fpcalib research compendium, v0.1.0, <doi:10.5281/zenodo.23379100>).

### Models

- [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  fits the multiple-methods false-positive occupancy model of Miller et
  al. (2011) in Stan, with any number of occupancy, detection, and
  false-positive covariates, an optional recorder random intercept, and
  configurable priors
  ([`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)).
  `fp = NULL` gives the standard occupancy model.
- Observation-confirmation design of Chambert et al. (2015): argument
  `confirmed` of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) and
  `frac_confirmed` of
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).
- Methods `print`, `summary`, `coef`, `predict`, `loo`, and
  `as_draws_df`.
- [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md):
  posterior predictive check with a Bayesian p-value.

### Simulation and calibration

- [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  generates homogeneous, unstructured, and structured false-positive
  regimes; with the same seed it returns the same datasets as the
  simulator of the fpcalib compendium.
- [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  and
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)
  report the coverage and bias of credible intervals for a design, with
  a [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method and
  a progress bar that shows the time remaining (needs ‘progressr’). They
  fit the recorder random intercept only when the simulated design has
  one, which avoids divergent transitions.

### Data and documentation

- `curiango` dataset: acoustic monitoring of the Common Pauraque in the
  Pantanal (Piatti et al. 2026).
- Articles: introduction, simulation, calibration, a case study with the
  `curiango` data, coming from
  [`unmarked::occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html),
  and the Stan model.
