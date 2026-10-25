# Convert data for unmarked::occuFP()

Turns an `fp_data` object into the `unmarkedFrameOccuFP` that
[`unmarked::occuFP()`](https://ecoverseR.github.io/unmarked/reference/occuFP.html)
fits by maximum likelihood, the frequentist counterpart of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) (see
[`vignette("unmarked")`](https://provetelab.org/crywolf/articles/unmarked.md)).
With `occasion = "occasion"` it is the reverse of
[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md).

## Usage

``` r
fp_to_unmarked(data, obs_covs = NULL, occasion = NULL)
```

## Arguments

- data:

  An `fp_data` object, for instance built with
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
  and [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
  from
  [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md).

- obs_covs:

  Optional names of numeric columns of `data$auto` to pass as
  observation covariates.

- occasion:

  Optional name of a column of `data$auto` (and of `data$verified`, if
  it has one) that numbers the occasions. Each occasion then goes to the
  column of its rank among all the values, so that occasions missing at
  a site leave `NA` in their place. Without it, a site's occasions are
  placed in the order of their rows.

## Value

An `unmarkedFrameOccuFP` object.

## Details

The detection matrix has one row per site, in the order of `data$sites`,
and the verified occasions first:
`type = c(n_verified, n_automated, 0)`, where `n_verified` and
`n_automated` are the numbers of verified and automated occasion columns
(see `occasion`); sites with fewer occasions get `NA`.

The observation covariate `METH` is `"verified"` on the verified columns
and `"auto"` on the automated ones (the reference level), because
`occuFP()` describes true detection by both methods with one formula:
use `detformula = ~ METH` for a separate detection probability of the
verified method, as
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) has.
The covariates in `obs_covs` exist for the automated occasions only;
they are 0 on the verified columns, because `occuFP()` drops every
observation whose covariate is missing. Centre them, so that 0 is a
typical value, and read a covariate in `detformula` as acting on the
automated method.

Recorders are not passed on (`occuFP()` has no random effects), and data
of the observation-confirmation design (`confirmed` in
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)) are
not supported.

## See also

[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
[`vignette("unmarked")`](https://provetelab.org/crywolf/articles/unmarked.md)

## Examples

``` r
sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
umf <- fp_to_unmarked(sim$data, occasion = "occasion")
unmarked::occuFP(detformula = ~ METH, FPformula = ~ 1, stateformula = ~ x,
                 data = umf)
#> 
#> Call:
#> unmarked::occuFP(detformula = ~METH, FPformula = ~1, stateformula = ~x, 
#>     data = umf)
#> 
#> Occupancy (logit-scale):
#>             Estimate    SE      z  P(>|z|)
#> (Intercept)    0.061 0.151  0.404 0.686182
#> x             -0.516 0.156 -3.303 0.000956
#> 
#> Detection (logit-scale):
#>              Estimate     SE     z  P(>|z|)
#> (Intercept)     0.381 0.0681 5.589 2.28e-08
#> METHverified    0.222 0.2499 0.888 3.75e-01
#> 
#> false positive (logit-scale):
#>  Estimate    SE     z  P(>|z|)
#>     -2.84 0.161 -17.7 6.34e-70
#> 
#> AIC: 2165.575 
#> Number of sites: 200
#> 
```
