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
if (FALSE) { # \dontrun{
s <- simulate_fp_dyn(nsite = 100, nseason = 4, seed = 1)
m <- fp_occu_dyn(s$data, col = ~ rain, ext = ~ x, random_recorder = FALSE)
new <- simulate(m, nsim = 2, seed = 1)
new[[1]]$truth$pars
} # }
```
