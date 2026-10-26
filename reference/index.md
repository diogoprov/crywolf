# Package index

## Data

Assemble automated and verified detections and site covariates.

- [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) :
  Assemble data for a false-positive occupancy model
- [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md)
  : Assemble data from detection matrices
- [`fp_to_unmarked()`](https://provetelab.org/crywolf/reference/fp_to_unmarked.md)
  : Convert data for unmarked::occuFP()

## Reading classifier output

Read BirdNET selection tables and reviewed clips, and build occasions.

- [`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md)
  : Read BirdNET selection tables
- [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)
  : Read the verdicts of a clip review
- [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md)
  : Read an export of the sheepdog package
- [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
  : Detections and confirmations per occasion
- [`fp_threshold_sweep()`](https://provetelab.org/crywolf/reference/fp_threshold_sweep.md)
  [`plot(`*`<fp_threshold_sweep>`*`)`](https://provetelab.org/crywolf/reference/fp_threshold_sweep.md)
  : Fit the model at several classifier thresholds

## Fitting

Fit the model and inspect, summarize, predict, and compare fits.

- [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) :
  Fit a Bayesian false-positive occupancy model
- [`fp_occu_mc()`](https://provetelab.org/crywolf/reference/fp_occu_mc.md)
  [`summary(`*`<fp_occu_mc>`*`)`](https://provetelab.org/crywolf/reference/fp_occu_mc.md)
  [`print(`*`<summary.fp_occu_mc>`*`)`](https://provetelab.org/crywolf/reference/fp_occu_mc.md)
  [`print(`*`<fp_occu_mc>`*`)`](https://provetelab.org/crywolf/reference/fp_occu_mc.md)
  [`as_draws_df(`*`<fp_occu_mc>`*`)`](https://provetelab.org/crywolf/reference/fp_occu_mc.md)
  : Fit the model to draws of the detection histories
- [`fp_draws_data()`](https://provetelab.org/crywolf/reference/fp_draws_data.md)
  : Detection histories over the draws of a classifier calibration
- [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)
  : Prior settings for fp_occu()
- [`print(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`summary(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`print(`*`<summary.fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`coef(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`predict(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`loo(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  [`as_draws_df(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  : Methods for fp_occu fits
- [`crywolf-output`](https://provetelab.org/crywolf/reference/crywolf-output.md)
  [`output-names`](https://provetelab.org/crywolf/reference/crywolf-output.md)
  : Output names and interpretation
- [`crywolf-api`](https://provetelab.org/crywolf/reference/crywolf-api.md)
  : Stability of the crywolf interface
- [`fp_loo_weights()`](https://provetelab.org/crywolf/reference/fp_loo_weights.md)
  : Model weights by stacking
- [`fp_gap_prior()`](https://provetelab.org/crywolf/reference/fp_gap_prior.md)
  : A prior for the gap from validated detections
- [`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)
  [`plot(`*`<fp_prior_predictive>`*`)`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)
  : Prior predictive check
- [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) :
  Posterior predictive check for a false-positive occupancy model
- [`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md)
  : Compile (or load) the Stan model
- [`fp_stan_file()`](https://provetelab.org/crywolf/reference/fp_stan_file.md)
  : Path to the Stan program

## Plots and site-level estimates

Coefficients, diagnostics, marginal effects, and occupancy of each site.

- [`plot(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/plot.fp_occu.md)
  : Plot a fitted model
- [`fp_effects()`](https://provetelab.org/crywolf/reference/fp_effects.md)
  [`plot(`*`<fp_effects>`*`)`](https://provetelab.org/crywolf/reference/fp_effects.md)
  : Marginal effect of a covariate
- [`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)
  [`plot(`*`<fp_sites>`*`)`](https://provetelab.org/crywolf/reference/fp_sites.md)
  : Site-level estimates
- [`plot(`*`<fp_ppc>`*`)`](https://provetelab.org/crywolf/reference/plot.fp_ppc.md)
  : Plot a posterior predictive check

## Simulation and calibration

Simulate false-positive regimes and check interval calibration for a
design.

- [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  : Simulate data for a false-positive occupancy model
- [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  : Simulation-based calibration check for a study design
- [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)
  : Coverage as a function of verification effort
- [`fp_calibrate_fit()`](https://provetelab.org/crywolf/reference/fp_calibrate_fit.md)
  : Calibration of a fitted model on data like the real ones
- [`fp_sbc()`](https://provetelab.org/crywolf/reference/fp_sbc.md)
  [`plot(`*`<fp_sbc>`*`)`](https://provetelab.org/crywolf/reference/fp_sbc.md)
  : Simulation-based calibration
- [`simulate(`*`<fp_occu>`*`)`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
  : Simulate data from a fitted model
- [`plot(`*`<fp_calibration>`*`)`](https://provetelab.org/crywolf/reference/plot.fp_calibration.md)
  : Plot the results of a calibration check

## Example data

Passive acoustic monitoring of the Common Pauraque in the Pantanal.

- [`curiango`](https://provetelab.org/crywolf/reference/curiango.md) :
  Passive acoustic monitoring of the Common Pauraque in the Pantanal

## Package

- [`crywolf`](https://provetelab.org/crywolf/reference/crywolf-package.md)
  [`crywolf-package`](https://provetelab.org/crywolf/reference/crywolf-package.md)
  : crywolf: Bayesian false-positive occupancy models with calibration
  checks
- [`reexports`](https://provetelab.org/crywolf/reference/reexports.md)
  [`loo`](https://provetelab.org/crywolf/reference/reexports.md)
  [`as_draws_df`](https://provetelab.org/crywolf/reference/reexports.md)
  : Objects exported from other packages
