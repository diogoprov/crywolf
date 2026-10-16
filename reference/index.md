# Package index

## Data

Assemble automated and verified detections and site covariates.

- [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) :
  Assemble data for a false-positive occupancy model
- [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md)
  : Assemble data from detection matrices

## Fitting

Fit the model and inspect, summarize, predict, and compare fits.

- [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) :
  Fit a Bayesian false-positive occupancy model
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
- [`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) :
  Posterior predictive check for a false-positive occupancy model
- [`fp_compile()`](https://provetelab.org/crywolf/reference/fp_compile.md)
  : Compile (or load) the Stan model
- [`fp_stan_file()`](https://provetelab.org/crywolf/reference/fp_stan_file.md)
  : Path to the Stan program

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
