# Methods for Royle-Nichols fits

[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html) for fits of
[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md).
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md) and
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) also
accept them.

[`loo()`](https://mc-stan.org/loo/reference/loo.html) computes the
leave-one-site-out predictive density from the marginal log-likelihood
of each site (abundance summed over 0..`n_max`), with the site effects
integrated out (see `integrate`); `reloo = TRUE` refits
[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md)
without each site with a high Pareto \\k\\.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) keeps the design
of the data (sites, covariates, occasions, recorders, which occasions
were verified or confirmed) and draws new abundances \\N_i \sim
\mathrm{Poisson}(\lambda_i)\\ and detections from the fitted model, as
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
does for
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md); with
`truth = "mean"`, new site effects on log abundance are drawn from a
normal distribution with the posterior mean of `sd_abund_site`.

## Usage

``` r
# S3 method for class 'fp_occu_rn'
loo(
  x,
  reloo = FALSE,
  k_threshold = 0.7,
  seed = NULL,
  integrate = TRUE,
  nodes = NULL,
  ...
)

# S3 method for class 'fp_occu_rn'
simulate(object, nsim = 1, seed = NULL, truth = c("mean", "posterior"), ...)

# S3 method for class 'fp_occu_rn'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_rn'
summary(object, prob = 0.95, ...)

# S3 method for class 'fp_occu_rn'
as_draws_df(x, ...)
```

## Arguments

- x, object:

  A fit of
  [`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md).

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

  Logical. With site random effects (`random_det`, `random_abund` or
  `random_fp`), integrate each site's effects out of its log-likelihood
  by Gauss-Hermite quadrature before the leave-one-out computation, as
  [`loo.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md)
  does. The default number of nodes per effect is 100 with one effect,
  30 with two and 12 with three; each node of the grid costs one pass
  over the abundances 0..`n_max`, so the time grows with `nodes`,
  `n_max` and the number of draws. In a test with 40 sites, the
  log-likelihood of a site moved by at most 0.005 from 30 to 60 nodes
  with two effects, and by at most 0.07 (median 1e-7) from 12 to 20
  nodes with three, where 20 nodes took about five times as long. Raise
  `nodes` when models are close; `integrate = FALSE` is quicker but
  gives many high Pareto \\k\\.

- nodes:

  Number of quadrature nodes for `integrate`, per site effect. The
  default, `NULL`, uses 100 with one site effect and 30 per effect with
  both (a grid of 900 nodes). With one effect and 100 nodes, the
  log-likelihood of a site was within 6e-4 of its value with 400 nodes
  in a test with 200 sites. The cost grows with the number of nodes in
  the grid, `nodes` with one effect and `nodes^2` with two.

- ...:

  Not used, except by
  [`loo()`](https://mc-stan.org/loo/reference/loo.html), which passes
  them to [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html).

- nsim:

  Number of datasets.

- truth:

  Parameter values to simulate from: `"mean"` (default) uses the
  posterior means for every dataset; `"posterior"` uses a different
  random posterior draw for each dataset, which carries the uncertainty
  of the fit into the simulations.

- digits:

  Digits to print.

- prob:

  Probability mass of the equal-tailed credible intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`summary.fp_occu` (see
[`summary.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md));
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
the posterior draws with readable names;
[`loo()`](https://mc-stan.org/loo/reference/loo.html) a `psis_loo`
object; [`simulate()`](https://rdrr.io/r/stats/simulate.html) a list of
`nsim` elements, each with `data` and `truth`, as from
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md),
with the coefficients of log abundance named `abund_*` and the
abundances `N` in `truth`.
