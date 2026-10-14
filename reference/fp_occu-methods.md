# Methods for fp_occu fits

Summaries, point estimates, predictions, and model comparison for
objects returned by
[`fp_occu()`](https://diogoprov.github.io/crywolf/reference/fp_occu.md).

## Usage

``` r
# S3 method for class 'fp_occu'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu'
summary(object, prob = 0.95, ...)

# S3 method for class 'summary.fp_occu'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu'
coef(object, type = c("all", "occ", "det", "fp"), ...)

# S3 method for class 'fp_occu'
predict(
  object,
  newdata,
  process = c("occ", "fp"),
  prob = 0.95,
  summary = TRUE,
  ...
)

# S3 method for class 'fp_occu'
loo(x, ...)

# S3 method for class 'fp_occu'
as_draws_df(x, ...)
```

## Arguments

- x, object:

  An `fp_occu` object.

- digits:

  Number of digits to print.

- ...:

  Unused, or passed on (see each method).

- prob:

  Width of the equal-tailed credible interval.

- type:

  Which coefficients to return: `"all"` (default), `"occ"`, `"det"`, or
  `"fp"`.

- newdata:

  Data frame with the covariates of the requested process. If missing,
  the sites used in the fit.

- process:

  `"occ"` for the occupancy probability or `"fp"` for the false-positive
  rate, both at the site level.

- summary:

  Logical. Return a summary per row (`TRUE`) or the matrix of posterior
  draws (rows are draws, columns are rows of `newdata`)?

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.fp_occu` whose element `table` is a tibble with the
posterior mean, standard deviation, credible interval, R-hat, and
effective sample sizes of every parameter.

[`coef()`](https://rdrr.io/r/stats/coef.html) returns a named vector of
posterior means.

[`predict()`](https://rdrr.io/r/stats/predict.html) returns a data frame
with the posterior mean and credible interval of the probability for
each row of `newdata`, or a matrix of draws.

[`loo()`](https://mc-stan.org/loo/reference/loo.html) returns a
`psis_loo` object from 'loo'. The pointwise log-likelihood is computed
per site, with the latent state marginalized, so the comparison is
leave-one-site-out.

[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
returns the posterior draws with the readable parameter names.
