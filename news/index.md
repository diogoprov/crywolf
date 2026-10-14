# Changelog

## crywolf 0.0.0.9000

- First version:
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
  [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  with methods,
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
  and
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md).
- The Stan model generalizes the one in the fpcalib research compendium:
  any number of occupancy, detection, and false-positive covariates, an
  optional recorder random intercept, and configurable priors.
