# Methods for fp_occu fits

Summaries, point estimates, predictions, and model comparison for
objects returned by
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

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
loo(x, reloo = FALSE, k_threshold = 0.7, seed = NULL, ...)

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
  the sites used in the fit; for a false-positive rate with covariates
  measured on each occasion, the automated occasions, with the
  covariates of their sites.

- process:

  `"occ"` for the occupancy probability or `"fp"` for the false-positive
  rate, both at the site level.

- summary:

  Logical. Return a summary per row (`TRUE`) or the matrix of posterior
  draws (rows are draws, columns are rows of `newdata`)?

- reloo:

  Logical. Refit the model without each site whose Pareto \\k\\ exceeds
  `k_threshold`, and replace the approximate leave-one-out value of that
  site by the exact one. Each such site costs one fit, with the sampler
  settings of `x`. A site that is the only one with its recorder, or
  with a level of a factor covariate, cannot be predicted from a fit
  without it; it keeps its approximate value, with a warning.

- k_threshold:

  Pareto \\k\\ above which a site is refitted.

- seed:

  Seed of the refits.

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
leave-one-site-out. Pareto \\k\\ above 0.7 marks sites where the
approximation is unreliable; `reloo = TRUE` computes their values
exactly, and their \\k\\ is then set to 0.

[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
returns the posterior draws with the readable parameter names.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
m                                   # posterior summary and diagnostics
summary(m, prob = 0.9)$table
coef(m, "occ")
predict(m, newdata = data.frame(x = c(-1, 0, 1)), process = "occ")
loo(m)
draws <- as_draws_df(m)
} # }
```
