# Fit a Bayesian false-positive occupancy model

Fits the multiple-methods false-positive occupancy model of Miller et
al. (2011) in Stan. Detections from an automated method may be false
positives, whereas detections from a verified subset are certain. The
occupancy probability depends on site covariates, true detection on
observation (or site) covariates and an optional recorder random
intercept, and the false-positive rate is either constant or a function
of site covariates.

## Usage

``` r
fp_occu(
  data,
  occ = ~1,
  det = ~1,
  fp = ~1,
  random_recorder = !is.null(data$recorder),
  priors = fp_priors(),
  chains = 4,
  parallel_chains = chains,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = NULL,
  refresh = 0,
  ...
)
```

## Arguments

- data:

  An `fp_data` object from
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
  [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
  or
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).

- occ:

  One-sided formula for the logit of occupancy, using columns of
  `data$sites`.

- det:

  One-sided formula for the logit of the true-detection probability of
  the automated method. It may use columns of `data$auto` and of
  `data$sites`. It must keep the intercept.

- fp:

  One-sided formula for the logit of the false-positive rate, using
  columns of `data$sites`. `~ 1` (the default) fits a constant rate;
  `~ w` lets the rate change along `w`. Use `NULL` to fit a standard
  occupancy model without false positives. The intercept is fixed by the
  model (see Details) and cannot be removed.

- random_recorder:

  Logical. Add a recorder random intercept to true detection? Defaults
  to `TRUE` when `data` has a recorder column.

- priors:

  Prior settings from
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

- chains, iter_warmup, iter_sampling, adapt_delta, seed, refresh:

  Passed to the `$sample()` method of 'cmdstanr'.

- parallel_chains:

  Number of chains run in parallel.

- ...:

  Further arguments passed to `$sample()`.

## Value

An object of class `fp_occu`, with methods for
[`base::print()`](https://rdrr.io/r/base/print.html),
[`base::summary()`](https://rdrr.io/r/base/summary.html),
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::predict()`](https://rdrr.io/r/stats/predict.html),
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html), and
[`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html).
The 'cmdstanr' fit is in element `fit`.

## Details

**Identifiability.** The intercept of the false-positive rate equals the
intercept of true detection minus a positive gap, \\\gamma_0 =
\alpha_1 - \Delta\\, so that at the reference values of the covariates
(all zero) a true detection is more likely than a false positive. This
removes the label switching between the occupied and unoccupied states.
Centre the covariates, so that the reference is an average site;
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
already returns them standardized.

**Verified data.** The verified detections anchor the false-positive
rate. With none, the model still runs, but the rate is informed only by
the priors and the gap, and the fit should not be trusted.

**When to model the false-positive rate.** In the simulations behind
this package, a constant rate kept nominal coverage of the occupancy
slope when false positives were homogeneous or varied at random among
sites, but not when they tracked an environmental gradient correlated
with the occupancy covariate; verifying more sites did not fix it.
Placing that gradient, or the occupancy covariate itself, in `fp` did.
Use
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
to check your own design.

**Parameter names.** Coefficients are labelled by process and column of
the design matrix: `occ_*` (occupancy), `det_*` (true detection), `fp_*`
(false positive, where `fp_(Intercept)` is \\\gamma_0\\), plus `fp_gap`
(\\\Delta\\), `sd_recorder`, `p_verified` (detection probability of the
verified method), `p11_ref` and `p10_ref` (true- and false-positive
probabilities at the reference), and `p10_bar` (false-positive
probability averaged over the sites).

## References

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
[doi:10.1890/10-1396.1](https://doi.org/10.1890/10-1396.1)

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 42)

# constant false-positive rate
m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1)
# false-positive rate that tracks the occupancy covariate
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x)

coef(m0)["occ_x"]; coef(m1)["occ_x"]   # truth: -0.8
loo::loo_compare(loo(m0), loo(m1))
} # }
```
