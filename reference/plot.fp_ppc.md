# Plot a posterior predictive check

Plot a posterior predictive check

## Usage

``` r
# S3 method for class 'fp_ppc'
plot(x, type = c("discrepancy", "histogram", "counts"), ...)
```

## Arguments

- x:

  An `fp_ppc` object.

- type:

  `"discrepancy"` (default) plots the replicated against the observed
  discrepancy for every draw, with the line of equality; the Bayesian
  p-value is the share of points above it. `"histogram"` shows the same
  as the distribution of their difference, replicated minus observed,
  with the p-value as the share above zero. `"counts"` compares the
  distribution of detections per site with the replicates, with
  [`bayesplot::ppc_bars()`](https://mc-stan.org/bayesplot/reference/PPC-discrete.html).

- ...:

  Passed to
  [`bayesplot::ppc_bars()`](https://mc-stan.org/bayesplot/reference/PPC-discrete.html).

## Value

A 'ggplot2' object.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
p <- fp_ppc(m, ndraws = 500, seed = 1)
plot(p)
plot(p, type = "counts")
} # }
```
