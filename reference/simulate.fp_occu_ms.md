# Simulate data from a fitted multispecies model

Generates new detections from an
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
fit, keeping the structure of the original data (sites, occasions and
their covariates, recorders, which occasions were verified or checked),
as
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
does for one species. **In development**, as
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).

## Usage

``` r
# S3 method for class 'fp_occu_ms'
simulate(object, nsim = 1, seed = NULL, truth = c("mean", "posterior"), ...)
```

## Arguments

- object:

  An `fp_occu_ms` object.

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
with the structure and the species column of the original data) and
`truth`: `pars`, the generating values of each species named as in
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
with the species in brackets (`occ_x[species]`,
`det_(Intercept)[species]`, `fp_(Intercept)[species]`, ...); `z`, the
latent states, a matrix of sites by species; `p10_site`, the
false-positive rate of each site and species, a matrix of the same
shape; and `p10_bar`, its mean over sites for each species. All the
species of a dataset take their values from the same posterior draw (or
from the posterior means).

## See also

[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
[`simulate_fp_ms()`](https://provetelab.org/crywolf/reference/simulate_fp_ms.md)

## Examples

``` r
# \donttest{
s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE,
                chains = 2, threads_per_chain = 1)
#> Parametrization (4 species at 60 sites): centered = FALSE, centered_occ = FALSE; see ?fp_occu_ms.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 99.7 seconds.
#> Chain 2 finished in 106.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 103.2 seconds.
#> Total execution time: 106.8 seconds.
#> 
#> Warning: Largest R-hat is 1.013 (> 1.01); run longer chains.
new <- simulate(m, nsim = 2, seed = 1)
new[[1]]$truth$pars
#>      occ_(Intercept)[sp1]                occ_x[sp1]      det_(Intercept)[sp1] 
#>               0.440132908              -0.104310404              -0.798581731 
#>              det_bio[sp1] verified_(Intercept)[sp1]       fp_(Intercept)[sp1] 
#>               0.009752225               0.365604973              -3.119346576 
#>      occ_(Intercept)[sp2]                occ_x[sp2]      det_(Intercept)[sp2] 
#>              -0.256445254               0.169746561              -0.202209963 
#>              det_bio[sp2] verified_(Intercept)[sp2]       fp_(Intercept)[sp2] 
#>               0.009752225               0.365604973              -2.402080349 
#>      occ_(Intercept)[sp3]                occ_x[sp3]      det_(Intercept)[sp3] 
#>              -0.274142055              -0.069720873               0.656982628 
#>              det_bio[sp3] verified_(Intercept)[sp3]       fp_(Intercept)[sp3] 
#>               0.009752225               0.365604973              -1.845296964 
#>      occ_(Intercept)[sp4]                occ_x[sp4]      det_(Intercept)[sp4] 
#>              -0.041047102               0.008685685              -0.078655141 
#>              det_bio[sp4] verified_(Intercept)[sp4]       fp_(Intercept)[sp4] 
#>               0.009752225               0.365604973              -5.663891029 
# }
```
