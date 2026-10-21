# Marginal effect of a covariate

Posterior of the occupancy probability (or the false-positive rate)
along one covariate, with the other covariates of that process held at
their means (numeric) or first level (factors) over the sites.

## Usage

``` r
fp_effects(
  object,
  covariate,
  process = c("occ", "fp"),
  values = NULL,
  transform = NULL,
  prob = 0.95
)

# S3 method for class 'fp_effects'
plot(x, xlab = x$covariate[1], ...)
```

## Arguments

- object:

  An `fp_occu` object.

- covariate:

  Name of a covariate of the process.

- process:

  `"occ"` (default) or `"fp"`.

- values:

  Values of the covariate, on the scale used in the model. The default
  spans its observed range in 100 steps, or takes every level of a
  factor.

- transform:

  Optional function that maps model-scale values to the scale shown on
  the x axis, for covariates standardized or transformed before fitting;
  for a covariate standardized with
  [`base::scale()`](https://rdrr.io/r/base/scale.html),
  `function(z) z * sd(raw) + mean(raw)`.

- prob:

  Width of the credible interval.

- x:

  An `fp_effects` object.

- xlab:

  Label of the x axis.

- ...:

  Unused.

## Value

A data frame of class `fp_effects` with the covariate value (`value`,
model scale), the value shown on the plot (`x`), and the posterior mean
and interval of the probability;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws it.

## See also

[`predict.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
plot(fp_effects(m, "x"))
plot(fp_effects(m, "x", process = "fp"))
} # }
```
