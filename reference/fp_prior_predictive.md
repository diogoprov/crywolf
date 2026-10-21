# Prior predictive check

Draws parameter values from the priors, simulates data with the design
of `data` (its sites, covariates, occasions, recorders, and
verification), and compares a few summaries of the simulated data with
the observed ones. If the observed values sit far in the tails of what
the priors produce, the priors (or the model) are at odds with the data
before any fitting.

## Usage

``` r
fp_prior_predictive(
  data,
  occ = ~1,
  det = ~1,
  fp = ~1,
  ver = ~1,
  markov = NULL,
  random_recorder = !is.null(data$recorder),
  random_fp = FALSE,
  priors = fp_priors(),
  ndraws = 200,
  seed = NULL
)

# S3 method for class 'fp_prior_predictive'
plot(x, ...)
```

## Arguments

- data:

  An `fp_data` object.

- occ, det, fp, ver, markov, random_recorder, random_fp, priors:

  As in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

- ndraws:

  Number of prior draws, one simulated dataset each.

- seed:

  Optional integer seed.

- x:

  An `fp_prior_predictive` object.

- ...:

  Unused.

## Value

An object of class `fp_prior_predictive`: a data frame with one row per
draw and one column per summary, with the observed values as the
attribute `observed`. [`print()`](https://rdrr.io/r/base/print.html)
shows the prior predictive intervals and where the observed values fall;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws them.

## Details

Each prior draw takes the intercepts and slopes from their Normal
priors, the gap from its truncated Normal prior, and the standard
deviations of the random effects from their half-Normal priors, with new
recorder and site effects drawn from them. Covariates are taken as they
are in `data`, so the priors are judged on the scale the model will use.
The summaries are the share of sites with at least one automated
detection (naive occupancy), the share of automated occasions with a
detection, the share of verified occasions with a detection, the mean
occupancy probability, and the mean false-positive rate over sites.

## See also

[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)

## Examples

``` r
sim <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
pp <- fp_prior_predictive(sim$data, occ = ~ x, fp = ~ 1,
                          random_recorder = FALSE, ndraws = 50, seed = 1)
pp
#> Prior predictive check (50 draws)
#> 
#>          summary     q05 median   q95 observed share_below_observed
#>  naive_occupancy 0.24750  0.845 1.000    0.700                 0.30
#>   automated_rate 0.03030  0.273 0.712    0.319                 0.58
#>    verified_rate 0.03594  0.208 0.532    0.208                 0.48
#>         mean_psi 0.22774  0.487 0.826       NA                   NA
#>          p10_bar 0.00465  0.107 0.686       NA                   NA
#> 
#> Observed values near 0 or 1 in `share_below_observed` sit in the tails of the prior predictive distribution.
```
