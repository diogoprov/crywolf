# crywolf 0.0.0.9000

* First version: `fp_data()`, `fp_data_wide()`, `fp_occu()` with methods,
  `fp_priors()`, `simulate_fp()`, `fp_calibrate()`, and
  `fp_verification_curve()`.
* The Stan model generalizes the one in the fpcalib research compendium:
  any number of occupancy, detection, and false-positive covariates, an
  optional recorder random intercept, and configurable priors.
