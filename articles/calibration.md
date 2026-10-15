# Checking calibration and how much to verify

A 95% credible interval should contain the true value in about 95% of
datasets. When the model is misspecified, it may not, and nothing in a
single fit will warn you. Simulation-based calibration answers the
question for a given design: simulate many datasets from known
parameters, fit the model to each, and count how often the intervals
cover the truth.

## The workflow

1.  Describe the design you can afford: number of sites, automated
    occasions, fraction of sites verified, and plausible detection and
    false-positive rates.
2.  Describe the false-positive structure you suspect: homogeneous,
    random among sites, or tracking an environmental gradient.
3.  Run
    [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
    with the candidate models.
4.  Read the coverage and bias of the quantity you care about, usually
    the occupancy slope.

``` r

library(crywolf)
fp_compile()   # compile once before starting parallel workers

future::plan(future::multisession, workers = 8)

design <- expand.grid(
  nsite = 300, nocc = 20, p10 = 0.05, frac_verified = 0.3,
  conf_effect = c(0, 0.5, 1.5), conf_cor = 0.7
)

cal <- fp_calibrate(
  design,
  models = list(constant = ~ 1, fp_on_x = ~ x),
  occ = ~ x, reps = 50, parallel = TRUE
)
summary(cal)
```

The code above is not run when this vignette is built; the output below
comes from the same call, saved in advance.
[`summary()`](https://rdrr.io/r/base/summary.html) prints, for the
occupancy slope, the coverage of the 95% intervals, its Monte Carlo
standard error, the bias, and a verdict:

    #> Calibration summary (95% intervals)
    #> 
    #> Occupancy slope:
    #>  cell    model nsite nocc  p10 frac_verified conf_effect conf_cor  n failed
    #>     1 constant   300   20 0.05           0.3         0.0      0.7 50      0
    #>     1  fp_on_x   300   20 0.05           0.3         0.0      0.7 50      0
    #>     2 constant   300   20 0.05           0.3         0.5      0.7 50      0
    #>     2  fp_on_x   300   20 0.05           0.3         0.5      0.7 50      0
    #>     3 constant   300   20 0.05           0.3         1.5      0.7 50      0
    #>     3  fp_on_x   300   20 0.05           0.3         1.5      0.7 50      0
    #>  coverage mc_se    bias rmse prop_divergent   verdict
    #>      1.00 0.000 0.06986 0.19           0.66        ok
    #>      0.98 0.020 0.00032 0.14           0.52        ok
    #>      0.96 0.028 0.07085 0.21           0.58        ok
    #>      0.96 0.028 0.01103 0.15           0.70        ok
    #>      0.22 0.059 0.42449 0.45           0.60 collapsed
    #>      0.62 0.069 0.24085 0.29           0.68 collapsed
    #> 
    #> All quantities: inspect the returned tibble.

![Coverage of the 95% interval of the occupancy slope against the
strength of the confusion, for the constant-rate model and for a
false-positive rate that depends on the occupancy covariate. Bars span
two Monte Carlo standard errors; the dashed line is the nominal
coverage.](calibration_files/figure-html/plot-confusion-1.png)

Coverage of the 95% interval of the occupancy slope against the strength
of the confusion, for the constant-rate model and for a false-positive
rate that depends on the occupancy covariate. Bars span two Monte Carlo
standard errors; the dashed line is the nominal coverage.

Each replicate costs one Stan fit per model, so this design needs 3
cells × 50 replicates × 2 models = 300 fits. With 50 replicates, the
Monte Carlo standard error of a coverage of 0.95 is about 0.03, enough
to tell a calibrated model from a broken one; raise `reps` only to pin a
borderline value. Every model is fitted to the same datasets, so
comparisons between models are paired.

The column `prop_divergent` is the share of fits with at least one
divergent transition. In this run most fits had a few divergences
(median 1 of 4000 draws), so check `max_rhat` and `divergences` in the
returned tibble before trusting a cell whose fits had many.

## How much to verify

Verification anchors the false-positive rate. It is natural to think
that verifying more sites will rescue a model that is not calibrated,
and
[`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)
tests that directly:

``` r

fp_verification_curve(
  nsite = 300, nocc = 20, p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
  frac_verified = c(0.1, 0.3, 0.5, 0.7, 0.9),
  models = list(constant = ~ 1, fp_on_x = ~ x),
  reps = 50, parallel = TRUE
)
```

    #> Calibration summary (95% intervals)
    #> 
    #> Occupancy slope:
    #>  cell    model frac_verified nsite nocc  p10 conf_effect conf_cor  n failed
    #>     1 constant           0.1   300   20 0.05         1.5      0.7 50      0
    #>     1  fp_on_x           0.1   300   20 0.05         1.5      0.7 50      0
    #>     2 constant           0.3   300   20 0.05         1.5      0.7 50      0
    #>     2  fp_on_x           0.3   300   20 0.05         1.5      0.7 50      0
    #>     3 constant           0.5   300   20 0.05         1.5      0.7 50      0
    #>     3  fp_on_x           0.5   300   20 0.05         1.5      0.7 50      0
    #>     4 constant           0.7   300   20 0.05         1.5      0.7 50      0
    #>     4  fp_on_x           0.7   300   20 0.05         1.5      0.7 50      0
    #>     5 constant           0.9   300   20 0.05         1.5      0.7 50      0
    #>     5  fp_on_x           0.9   300   20 0.05         1.5      0.7 50      0
    #>  coverage mc_se bias rmse prop_divergent   verdict
    #>      0.06 0.034 0.55 0.57           0.82 collapsed
    #>      0.38 0.069 0.39 0.43           0.82 collapsed
    #>      0.12 0.046 0.45 0.48           0.64 collapsed
    #>      0.48 0.071 0.29 0.32           0.68 collapsed
    #>      0.22 0.059 0.40 0.42           0.58 collapsed
    #>      0.64 0.068 0.21 0.26           0.68 collapsed
    #>      0.10 0.042 0.40 0.42           0.62 collapsed
    #>      0.64 0.068 0.21 0.24           0.58 collapsed
    #>      0.24 0.060 0.39 0.42           0.64 collapsed
    #>      0.64 0.068 0.21 0.25           0.60 collapsed
    #> 
    #> All quantities: inspect the returned tibble.

![Coverage of the 95% interval of the occupancy slope against the
fraction of sites verified, under structured confusion (effect 1.5,
correlation
0.7).](calibration_files/figure-html/plot-verification-1.png)

Coverage of the 95% interval of the occupancy slope against the fraction
of sites verified, under structured confusion (effect 1.5, correlation
0.7).

## What to expect

In the simulation study that motivated this package (Provete et al., in
preparation), at a base false-positive rate of 0.05:

- When false positives were homogeneous, the constant-rate model kept
  nominal coverage of the occupancy slope. When they varied at random
  among sites, the estimate of the false-positive rate became
  unreliable, but coverage of the slope fell only to 0.81 (across all
  cells of that regime).
- When the false-positive rate tracked a gradient correlated with the
  occupancy covariate (effect 1.5 on the logit scale, correlation 0.7),
  coverage of the slope fell to 0.08–0.20, and verifying up to 90% of
  the sites did not restore it.
- Letting the false-positive rate depend on the occupancy covariate
  (`fp = ~ x`) recovered part of the coverage (0.48–0.66) and improved
  with verification; with the true gradient (`fp = ~ w`), coverage
  returned to 0.90–0.97.
- Weak confusion (effects of 0.5 or less) left the constant-rate model
  calibrated at every correlation, including when false positives
  depended on the occupancy covariate itself.

So, if
[`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md)
shows the constant-rate model stuck below nominal coverage, verifying
more is not the remedy. Model the false-positive rate instead, and
verify enough to tell its effect apart from the occupancy effect.

## Plotting

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
coverage of the occupancy slope against the design column that varies
(here, the fraction of sites verified), one colour per model, with
intervals of two Monte Carlo standard errors and a dashed line at the
nominal level. `type = "bias"` plots the bias instead, and `parameter`
selects another scored quantity. The result is a ‘ggplot2’ object that
can be modified as usual.

``` r

plot(cal)
plot(cal, type = "bias")
plot(cal, parameter = "p10_bar")
```
