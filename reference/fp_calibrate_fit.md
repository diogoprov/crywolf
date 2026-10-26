# Calibration of a fitted model on data like the real ones

Asks whether the credible intervals of a model would cover the truth in
data like the ones it was fitted to. The function simulates datasets
from the fit with
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md),
so that they keep the design of the real data (sites, covariates,
occasions, verification, recorders) and follow the fitted model, refits
each candidate model to every dataset, and scores the coverage and bias
of their intervals, as
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
does for designed simulations.

## Usage

``` r
fp_calibrate_fit(
  fit,
  models = list(fitted = fit$formulas$fp),
  reps = 50,
  truth = c("mean", "posterior"),
  prob = 0.95,
  base_seed = 1L,
  parallel = FALSE,
  fit_args = list(chains = 4, parallel_chains = 1, iter_warmup = 1000, iter_sampling =
    1000, adapt_delta = 0.95)
)
```

## Arguments

- fit:

  An `fp_occu` object: the generating model.

- models:

  Named list of false-positive formulas (or `NULL` for no false
  positives) to refit. The default refits the generating model. The
  occupancy and detection formulas, the recorder effect, and the priors
  are those of `fit`.

- reps:

  Number of simulated datasets.

- truth:

  `"mean"` (default) or `"posterior"`; see
  [`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md).

- prob:

  Width of the credible intervals.

- base_seed:

  Seeds are `base_seed + 1, ..., base_seed + reps`.

- parallel, fit_args:

  As in
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md).

## Value

An object of class `fp_calibration`, as from
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
with a single design cell.

## Details

The scored quantities are the occupancy coefficients, the detection
coefficients, the intercept of the verified method, the mean
false-positive rate over sites (`p10_bar`, when both the generating and
the refitted model have false positives), and the slopes of the
false-positive rate in the refitted model. A slope that is not in the
generating model has a true value of zero.

The most useful comparison pits the model you fitted against a simpler
one. For example, if a model with `fp = ~ temp + rh` estimated clear
effects of the covariates on the false-positive rate, refitting
`fp = ~ 1` to data simulated from it shows how much a constant rate
would bias the occupancy effects at your sample size and verification
effort. This answers the question for the fitted structure only: a
dependence of the false-positive rate on something you did not measure
is not part of the generating model.

Each replicate refits every model, so the cost is that of `reps` times
the number of models fits of the size of the real data. Use `parallel`,
and `threads_per_chain` in `fit_args` for large datasets.

## See also

[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 1)
m_x <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
cal <- fp_calibrate_fit(m_x, models = list(constant = ~ 1, fp_on_x = ~ x),
                        reps = 20)
summary(cal)
plot(cal, parameter = "occ_x")
} # }
```
