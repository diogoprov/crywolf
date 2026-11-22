# Simulate data from a fitted dynamic model

Generates new detections from an
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
fit, keeping the structure of the original data (sites, seasons,
occasions and their covariates, which site-seasons were verified and
which occasions were checked), as
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
does for a single season. The latent states are drawn anew over seasons,
from initial occupancy, colonization and extinction, and the detections
given them.

## Usage

``` r
# S3 method for class 'fp_occu_dyn'
simulate(object, nsim = 1, seed = NULL, truth = c("mean", "posterior"), ...)
```

## Arguments

- object:

  An `fp_occu_dyn` object.

- nsim:

  Number of datasets.

- seed:

  Optional integer seed.

- truth:

  Parameter values to simulate from: `"mean"` (default) uses the
  posterior means for every dataset; `"posterior"` uses a different
  random posterior draw for each dataset, which carries the uncertainty
  of the fit into the simulations.

- ...:

  Unused.

## Value

A list of `nsim` elements, each a list with `data` (an `fp_data` object
with the structure of the original data) and `truth` (the generating
values named as the output of
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
in `pars`, the latent states `z` as a matrix with sites in rows and
seasons in columns, the false-positive rate of each site-season
`p10_site`, and its mean `p10_bar`).

## See also

[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md),
[`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md)

## Examples

``` r
# \donttest{
s <- simulate_fp_dyn(nsite = 100, nseason = 4, seed = 1)
m <- fp_occu_dyn(s$data, col = ~ rain, ext = ~ x, random_recorder = FALSE,
                 chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 2 finished in 17.0 seconds.
#> Chain 1 finished in 17.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 17.2 seconds.
#> Total execution time: 17.4 seconds.
#> 
new <- simulate(m, nsim = 2, seed = 1)
new[[1]]$truth$pars
#>     psi1_(Intercept)      col_(Intercept)             col_rain 
#>           0.06125464          -1.43877557           0.02254640 
#>      ext_(Intercept)                ext_x      det_(Intercept) 
#>          -1.75610400           0.48620397           0.44871645 
#> verified_(Intercept)       fp_(Intercept) 
#>           0.80990599          -2.86552348 
# }
```
