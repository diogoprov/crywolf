# Changelog

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
