# Simulation-based calibration check for a study design

Simulates many datasets from a design with
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
fits one or more false-positive models to each with
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), and
scores how often the credible intervals contain the true values and how
biased the posterior means are. Use it before collecting data, to ask
whether the correction will be calibrated for the number of sites,
occasions, and verified sites you can afford, and whether a constant
false-positive rate is enough under the confusion you suspect.

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
    1000, adapt_delta = 0.95)
)
```

## Arguments

- design:

  A data frame (or list) whose columns are arguments of
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md);
  each row is one cell of the design. Arguments not given take the
  defaults of
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).

- models:

  Named list of formulas for the false-positive rate, one model per
  element, as in the `fp` argument of
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  The default fits only the constant-rate model.

- occ, det:

  Formulas for occupancy and true detection, shared by all models.
  Simulated data have the occupancy covariate `x` and the detection
  covariates `bio` and `aci`.

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

## See also

[`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
design <- data.frame(p10 = 0.05, conf_effect = c(0, 1.5), conf_cor = 0.7)
cal <- fp_calibrate(design,
                    models = list(constant = ~ 1, fp_on_x = ~ x),
                    occ = ~ x, reps = 50)
summary(cal)
} # }
```
