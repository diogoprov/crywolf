# crywolf 0.0.0.9000

* `fp_ppc()`: posterior predictive check with a Bayesian p-value.
* `plot()` method for `fp_calibration` objects.
* Observation-confirmation design of Chambert et al. (2015): argument
  `confirmed` of `fp_data()` and `frac_confirmed` of `simulate_fp()`.
* `curiango` dataset (Piatti et al. 2026).
* Progress: `fp_occu()` reports chain progress in interactive sessions
  (`refresh`), and `fp_calibrate()` and `fp_verification_curve()` show a
  progress bar with the estimated time remaining, also in parallel (needs
  'progressr').

* First version: `fp_data()`, `fp_data_wide()`, `fp_occu()` with methods,
  `fp_priors()`, `simulate_fp()`, `fp_calibrate()`, and
  `fp_verification_curve()`.
* The Stan model generalizes the one in the fpcalib research compendium:
  any number of occupancy, detection, and false-positive covariates, an
  optional recorder random intercept, and configurable priors.
