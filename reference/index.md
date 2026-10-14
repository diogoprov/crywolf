# Package index

## Data

Assemble automated and verified detections and site covariates.

- [`fp_data()`](https://diogoprov.github.io/crywolf/reference/fp_data.md)
  : Assemble data for a false-positive occupancy model
- [`fp_data_wide()`](https://diogoprov.github.io/crywolf/reference/fp_data_wide.md)
  : Assemble data from detection matrices

## Fitting

Fit the model and inspect, summarize, predict, and compare fits.

- [`fp_occu()`](https://diogoprov.github.io/crywolf/reference/fp_occu.md)
  : Fit a Bayesian false-positive occupancy model
- [`fp_priors()`](https://diogoprov.github.io/crywolf/reference/fp_priors.md)
  : Prior settings for fp_occu()
- [`print(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`summary(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`print(`*`<summary.fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`coef(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`predict(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`loo(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  [`as_draws_df(`*`<fp_occu>`*`)`](https://diogoprov.github.io/crywolf/reference/fp_occu-methods.md)
  : Methods for fp_occu fits
- [`fp_ppc()`](https://diogoprov.github.io/crywolf/reference/fp_ppc.md)
  : Posterior predictive check for a false-positive occupancy model
- [`fp_compile()`](https://diogoprov.github.io/crywolf/reference/fp_compile.md)
  : Compile (or load) the Stan model
- [`fp_stan_file()`](https://diogoprov.github.io/crywolf/reference/fp_stan_file.md)
  : Path to the Stan program

## Simulation and calibration

Simulate false-positive regimes and check interval calibration for a
design.

- [`simulate_fp()`](https://diogoprov.github.io/crywolf/reference/simulate_fp.md)
  : Simulate data for a false-positive occupancy model
- [`fp_calibrate()`](https://diogoprov.github.io/crywolf/reference/fp_calibrate.md)
  : Simulation-based calibration check for a study design
- [`fp_verification_curve()`](https://diogoprov.github.io/crywolf/reference/fp_verification_curve.md)
  : Coverage as a function of verification effort

## Package

- [`crywolf`](https://diogoprov.github.io/crywolf/reference/crywolf-package.md)
  [`crywolf-package`](https://diogoprov.github.io/crywolf/reference/crywolf-package.md)
  : crywolf: Bayesian false-positive occupancy models with calibration
  checks
- [`reexports`](https://diogoprov.github.io/crywolf/reference/reexports.md)
  [`loo`](https://diogoprov.github.io/crywolf/reference/reexports.md)
  [`as_draws_df`](https://diogoprov.github.io/crywolf/reference/reexports.md)
  : Objects exported from other packages
