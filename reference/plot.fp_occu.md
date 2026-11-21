# Plot a fitted model

Posterior intervals of the coefficients, or sampler diagnostics.

## Usage

``` r
# S3 method for class 'fp_occu'
plot(
  x,
  type = c("intervals", "trace", "rank", "pairs"),
  pars = NULL,
  prob = 0.95,
  ...
)
```

## Arguments

- x:

  An `fp_occu` object.

- type:

  `"intervals"` (default) draws the posterior median with 50% and `prob`
  intervals of every coefficient on the logit scale, one panel per
  process, and the standard deviations of the random effects. `"trace"`,
  `"rank"`, and `"pairs"` draw sampler diagnostics with 'bayesplot':
  trace plots and rank histograms by chain, and pairs plots with
  divergent transitions marked.

- pars:

  Output names of the parameters for the diagnostics. The default takes
  the occupancy slopes and the parameters that are hardest to sample in
  this model: the false-positive intercept, the gap, and the standard
  deviations of the random effects.

- prob:

  Width of the outer interval for `type = "intervals"`.

- ...:

  Passed to the 'bayesplot' function for the diagnostics.

## Value

A 'ggplot2' object (for `"pairs"`, a 'bayesplot' grid).

## See also

[`fp_effects()`](https://provetelab.org/crywolf/reference/fp_effects.md),
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md),
[`plot.fp_ppc()`](https://provetelab.org/crywolf/reference/plot.fp_ppc.md)

## Examples

``` r
# \donttest{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE,
             chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 18.4 seconds.
#> Chain 2 finished in 18.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 18.4 seconds.
#> Total execution time: 18.5 seconds.
#> 
plot(m)

plot(m, type = "trace")
#> No divergences to plot.

plot(m, type = "pairs", pars = c("fp_(Intercept)", "fp_gap"))
#> Error in suggested_package("gridExtra"): Please install the gridExtra package to use this function.
# }
```
