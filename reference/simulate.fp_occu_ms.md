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
if (FALSE) { # \dontrun{
s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE)
new <- simulate(m, nsim = 2, seed = 1)
new[[1]]$truth$pars
} # }
```
