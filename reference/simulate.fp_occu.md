# Simulate data from a fitted model

Generates new detections from an `fp_occu` fit, keeping everything else
of the original data: the sites and their covariates, the automated
occasions with their observation covariates and recorders, which sites
and occasions were verified, and which occasions were checked under the
observation-confirmation design. Only the latent occupancy states and
the detections are drawn anew, from the fitted model.

## Usage

``` r
# S3 method for class 'fp_occu'
simulate(object, nsim = 1, seed = NULL, truth = c("mean", "posterior"), ...)
```

## Arguments

- object:

  An `fp_occu` object.

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
values on the scale of the output of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), the
latent states `z`, the false-positive rate of each site `p10_site`, and
its mean `p10_bar`).

## Details

The recorder effects, when the model has them, are kept at their values
in the chosen draw (or their posterior means), so the simulated data
have the same recorders with the same effects. Site effects on the
false-positive rate (`fp_occu(random_fp = TRUE)`) are kept at their
values in the draw with `truth = "posterior"`; with `truth = "mean"`,
whose site effects would be shrunk towards zero, new ones are drawn from
a normal distribution with the posterior mean of `sd_fp_site` as
standard deviation. Under the multiple-methods design, an automated
occasion at an occupied site is a detection with probability \\p\_{11}\\
and at an unoccupied site with probability \\p\_{10}\\; verified
occasions are detections with probability \\r\_{11}\\ at occupied sites
only. Under the observation-confirmation design, true and false
observations are drawn independently, as in
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).

## See also

[`fp_calibrate_fit()`](https://provetelab.org/crywolf/reference/fp_calibrate_fit.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
new <- simulate(m, nsim = 2, seed = 1)
new[[1]]$data
new[[1]]$truth$p10_bar
} # }
```
