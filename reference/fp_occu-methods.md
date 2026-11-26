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
loo(
  x,
  reloo = FALSE,
  k_threshold = 0.7,
  seed = NULL,
  integrate = TRUE,
  nodes = NULL,
  ...
)

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

- integrate:

  Logical. With a site random effect on the false-positive rate
  (`random_fp = TRUE`), on true detection (`random_det = TRUE`), or
  both, integrate each site's effects out of its log-likelihood (by
  Gauss-Hermite quadrature) before the leave-one-out computation. The
  default, `TRUE`, is the right choice for leaving a site out: a
  left-out site's effects are not known from its own data. With `FALSE`,
  the log-likelihood is conditional on the fitted site effects, which
  gives many sites with high Pareto \\k\\. Without a site random effect,
  `integrate` has no effect.

- nodes:

  Number of quadrature nodes for `integrate`, per site effect. The
  default, `NULL`, uses 100 with one site effect and 30 per effect with
  both (a grid of 900 nodes). With one effect and 100 nodes, the
  log-likelihood of a site was within 6e-4 of its value with 400 nodes
  in a test with 200 sites. The cost grows with the number of nodes in
  the grid, `nodes` with one effect and `nodes^2` with two.

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
# \donttest{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
             chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 11.5 seconds.
#> Chain 2 finished in 11.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 11.6 seconds.
#> Total execution time: 11.7 seconds.
#> 
m                                   # posterior summary and diagnostics
#> Bayesian false-positive occupancy model (multiple-methods design)
#> 
#>   occupancy:       ~x 
#>   true detection:  ~1   
#>   false positive:  ~1  
#>   data:            200 sites, 2000 automated and 138 verified occasions
#> 
#>              variable    mean      sd  lower   upper rhat ess_bulk ess_tail
#>       occ_(Intercept)  0.0727 0.15269 -0.235  0.3758    1     2590     1178
#>                 occ_x -0.5195 0.15355 -0.813 -0.2350    1     2788     1473
#>       det_(Intercept)  0.3772 0.06972  0.241  0.5166    1     1996     1254
#>  verified_(Intercept)  0.5934 0.23311  0.138  1.0420    1     1938     1153
#>            p_verified  0.6424 0.05294  0.535  0.7392    1     1938     1153
#>               p11_ref  0.5931 0.01681  0.560  0.6263    1     1996     1254
#>        fp_(Intercept) -2.8455 0.16104 -3.178 -2.5467    1     2839     1490
#>                fp_gap  3.2227 0.16532  2.914  3.5753    1     2695     1488
#>               p10_ref  0.0555 0.00833  0.040  0.0727    1     2839     1490
#>               p10_bar  0.0555 0.00833  0.040  0.0727    1     2839     1490
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.003, min bulk/tail ESS 1938/1153
summary(m, prob = 0.9)$table
#> # A tibble: 10 × 8
#>    variable                mean      sd   lower   upper  rhat ess_bulk ess_tail
#>    <chr>                  <dbl>   <dbl>   <dbl>   <dbl> <dbl>    <dbl>    <dbl>
#>  1 occ_(Intercept)       0.0727 0.153   -0.186   0.324  1.00     2590.    1178.
#>  2 occ_x                -0.520  0.154   -0.776  -0.269  1.00     2788.    1473.
#>  3 det_(Intercept)       0.377  0.0697   0.260   0.492  1.00     1996.    1254.
#>  4 verified_(Intercept)  0.593  0.233    0.218   0.988  1.00     1938.    1153.
#>  5 p_verified            0.642  0.0529   0.554   0.729  1.00     1938.    1153.
#>  6 p11_ref               0.593  0.0168   0.565   0.621  1.00     1996.    1254.
#>  7 fp_(Intercept)       -2.85   0.161   -3.12   -2.59   1.000    2839.    1490.
#>  8 fp_gap                3.22   0.165    2.96    3.51   1.00     2695.    1488.
#>  9 p10_ref               0.0555 0.00833  0.0421  0.0696 1.000    2839.    1490.
#> 10 p10_bar               0.0555 0.00833  0.0421  0.0696 1.000    2839.    1490.
coef(m, "occ")
#> occ_(Intercept)           occ_x 
#>      0.07265452     -0.51954884 
predict(m, newdata = data.frame(x = c(-1, 0, 1)), process = "occ")
#>    x  estimate     lower     upper
#> 1 -1 0.6424340 0.5427749 0.7340527
#> 2  0 0.5180488 0.4415028 0.5928556
#> 3  1 0.3913724 0.2912191 0.4978346
loo(m)
#> 
#> Computed from 2000 by 200 log-likelihood matrix.
#> 
#>          Estimate   SE
#> elpd_loo  -1083.0 43.4
#> p_loo         5.3  0.4
#> looic      2165.9 86.8
#> ------
#> MCSE of elpd_loo is 0.0.
#> MCSE and ESS estimates assume MCMC draws (r_eff in [0.9, 1.6]).
#> 
#> All Pareto k estimates are good (k < 0.7).
#> See help('pareto-k-diagnostic') for details.
draws <- as_draws_df(m)
# }
```
