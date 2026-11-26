# Simulation-based calibration check for a study design

Simulates many datasets from a design with
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
fits one or more false-positive models to each with
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), and
scores how often the credible intervals contain the true values and how
biased the posterior means are. With `dynamic = TRUE`, the designs are
multi-season: they are simulated with
[`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md)
and fitted with
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md).
Use it before collecting data, to ask whether the correction will be
calibrated for the number of sites, occasions, and verified sites you
can afford, and whether a constant false-positive rate is enough under
the confusion you suspect.

## Usage

``` r
fp_calibrate(
  design,
  models = list(constant = ~1),
  occ = ~x,
  det = ~1,
  reps = 50,
  prob = 0.95,
  base_seed = 1L,
  parallel = FALSE,
  fit_args = list(chains = 4, parallel_chains = 1, iter_warmup = 1000, iter_sampling =
    1000, adapt_delta = 0.95),
  dynamic = FALSE,
  psi1 = ~x,
  col = ~1,
  ext = ~1
)
```

## Arguments

- design:

  A data frame (or list) whose columns are arguments of
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  (of
  [`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md)
  with `dynamic = TRUE`); each row is one cell of the design. Arguments
  not given take the defaults of the simulator.

- models:

  Named list of formulas for the false-positive rate, one model per
  element, as in the `fp` argument of
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  The default fits only the constant-rate model.

- occ, det:

  Formulas for occupancy and true detection, shared by all models.
  Simulated data have the occupancy covariate `x` and the detection
  covariates `bio` and `aci`. With `dynamic = TRUE`, `occ` is not used,
  and the data have no detection covariates.

- reps:

  Number of replicate datasets per cell. Each costs one Stan fit per
  model. Fifty replicates give a Monte Carlo standard error of about
  0.03 for a coverage of 0.95.

- prob:

  Width of the credible intervals scored.

- base_seed:

  Integer. Replicate `j` (counting over cells and replicates) uses seed
  `base_seed + j`, so results are reproducible and every model sees
  exactly the same datasets.

- parallel:

  Logical. Fit replicates in parallel with 'furrr'. Set the plan first,
  e.g. `future::plan(future::multisession, workers = 8)`, and keep
  `parallel_chains = 1` in `fit_args`. In interactive sessions with
  'progressr' installed, a progress bar shows the estimated time
  remaining, also when fitting in parallel.

- fit_args:

  List of further arguments for
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  Unless set here, `random_recorder` is `TRUE` only when the design
  simulates a recorder effect (`sd_recorder > 0`): fitting one whose
  true standard deviation is zero produces divergent transitions. With
  the chains of each fit run one after another, the cost is proportional
  to `chains`: `chains = 2` takes half the time of four and still gives
  an effective sample size of about 2000 per parameter in designs like
  the default one, ample for the 95% intervals that are scored. Each fit
  runs on one thread per chain (`threads_per_chain = 1`, unlike
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
  whose default is `"auto"`): `"auto"` divides the cores among the
  chains of one fit, not among parallel workers, and one thread keeps
  every replicate reproducible by its seed.

- dynamic:

  Logical. `TRUE` calibrates a multi-season design: the columns of
  `design` are arguments of
  [`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md),
  and each dataset is fitted with
  [`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md).

- psi1, col, ext:

  Formulas for initial occupancy, colonization and extinction, shared by
  all models, used only with `dynamic = TRUE`. Simulated dynamic data
  have the site covariate `x` and the site-season covariate `rain` (see
  [`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md)).
  Include the covariates of the design: with `b_col` not zero, for
  instance, use `col = ~ rain`.

## Value

A tibble of class `fp_calibration`, one row per replicate, model, and
scored quantity, with the truth, the posterior mean and interval,
whether the interval covered the truth, and sampler diagnostics. Use
[`base::summary()`](https://rdrr.io/r/base/summary.html) for coverage,
its Monte Carlo standard error, bias, and RMSE per cell and model.

## Details

The scored quantities are the occupancy intercept and slope
(`occ_(Intercept)`, `occ_x`), the false-positive rate averaged over
sites (`p10_bar`, against the realized mean rate), the detection
probability of the verified method (`p_verified`), and, for models that
include the confusing gradient `w`, its effect on the false-positive
rate (`fp_w`).

**Dynamic designs.** With `dynamic = TRUE` the scored quantities are the
intercepts of initial occupancy, colonization and extinction
(`psi1_(Intercept)`, `col_(Intercept)`, `ext_(Intercept)`), their slopes
when the formulas include the covariates of the simulator (`psi1_x`,
`col_rain`, `ext_x`), the false-positive rate (`p10_bar`, against the
constant `p10` of the design) and the detection probability of the
verified method (`p_verified`). Colonization and extinction are hard to
tell from false positives in seasons without verified data, so designs
with few verified site-seasons are the ones worth checking.

## See also

[`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
[`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md)

## Examples

``` r
# \donttest{
design <- data.frame(p10 = 0.05, conf_effect = c(0, 1.5), conf_cor = 0.7)
cal <- fp_calibrate(design,
                    models = list(constant = ~ 1, fp_on_x = ~ x),
                    occ = ~ x, reps = 5,
                    fit_args = list(chains = 2, threads_per_chain = 1))
summary(cal)   # use 50 or more replicates in practice
#> Calibration summary (95% intervals)
#> 
#> Occupancy slope:
#>  cell    model  p10 conf_effect conf_cor n failed coverage mc_se   bias  rmse
#>     1 constant 0.05         0.0      0.7 5      0      1.0  0.00 -0.014 0.076
#>     1  fp_on_x 0.05         0.0      0.7 5      0      1.0  0.00 -0.014 0.075
#>     2 constant 0.05         1.5      0.7 5      0      0.0  0.00  0.442 0.444
#>     2  fp_on_x 0.05         1.5      0.7 5      0      0.2  0.18  0.315 0.316
#>  prop_divergent   verdict
#>               0        ok
#>               0        ok
#>               0 collapsed
#>               0 collapsed
#> 
#> All quantities: inspect the returned tibble.

# a dynamic design: four seasons, colonization driven by rain
cal_dyn <- fp_calibrate(data.frame(nseason = 4, b_col = 1),
                        dynamic = TRUE, col = ~ rain, reps = 5,
                        fit_args = list(chains = 2, threads_per_chain = 1))
summary(cal_dyn)
#> Calibration summary (95% intervals)
#> 
#> Slopes of initial occupancy, colonization and extinction:
#>  cell    model parameter nseason b_col n failed coverage mc_se   bias rmse
#>     1 constant  col_rain       4     1 5      0        1     0 -0.054 0.38
#>     1 constant    psi1_x       4     1 5      0        1     0  0.011 0.23
#>  prop_divergent verdict
#>               0      ok
#>               0      ok
#> 
#> All quantities: inspect the returned tibble.
# }
```
