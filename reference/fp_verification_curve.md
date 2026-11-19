# Coverage as a function of verification effort

Runs
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
over a grid of verification fractions and prints the coverage of the
occupancy slope against the fraction of sites verified. If the
constant-rate model does not approach nominal coverage as verification
grows, the false positives are probably structured, and the remedy is to
model the false-positive rate (for instance on the occupancy covariate),
not to verify more.

## Usage

``` r
fp_verification_curve(
  ...,
  frac_verified = c(0.1, 0.3, 0.5, 0.7, 0.9),
  models = list(constant = ~1),
  occ = ~x,
  det = ~1,
  reps = 50,
  base_seed = 1L,
  parallel = FALSE,
  fit_args = list(chains = 4, parallel_chains = 1, iter_warmup = 1000, iter_sampling =
    1000, adapt_delta = 0.95)
)
```

## Arguments

- ...:

  Arguments of
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  that define the scenario (held fixed), e.g.
  `p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7`.

- frac_verified:

  Fractions of sites verified to sweep.

- models, occ, det, reps, base_seed, parallel, fit_args:

  As in
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md).

## Value

The `fp_calibration` object (invisibly), after printing the curve.

## Examples

``` r
# \donttest{
fp_verification_curve(p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
                      models = list(constant = ~ 1, fp_on_x = ~ x),
                      frac_verified = c(0.1, 0.5), reps = 5,
                      fit_args = list(chains = 2, threads_per_chain = 1))
#> 
#> Coverage of the occupancy slope against verification effort
#>     model frac_verified coverage mc_se bias
#>  constant           0.1      0.0  0.00 0.46
#>  constant           0.5      0.0  0.00 0.42
#>   fp_on_x           0.1      0.6  0.22 0.25
#>   fp_on_x           0.5      0.4  0.22 0.29
#> 
#>   The first model stays below nominal coverage at every verification 
#>   level. Verifying more will not fix it; model the false-positive rate.
# }
```
