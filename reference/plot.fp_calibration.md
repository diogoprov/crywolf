# Plot the results of a calibration check

Plots the coverage (or bias) of one scored quantity against the column
of the design that varies, one colour per model. Coverage is shown with
an interval of two Monte Carlo standard errors and a dashed line at the
nominal level; bias, with a dashed line at zero.

## Usage

``` r
# S3 method for class 'fp_calibration'
plot(
  x,
  parameter = "occ_x",
  type = c("coverage", "bias"),
  by = NULL,
  verified = c("fraction", "sites", "occasions"),
  target = NULL,
  ...
)
```

## Arguments

- x:

  An `fp_calibration` object from
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  or
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md).

- parameter:

  Scored quantity to plot, by default the occupancy slope (`"occ_x"`);
  for a dynamic calibration without it, the first slope scored, or else
  the colonization intercept. See
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  for the others.

- type:

  `"coverage"` (default) or `"bias"`.

- by:

  Name of the design column for the horizontal axis. By default, the
  first design column that varies among cells; further varying columns
  become panels. With no varying column, the models are on the axis.

- verified:

  How to show the verification effort when `by` is `"frac_verified"`: as
  the fraction of sites (default), as the expected number of verified
  sites (`"sites"`, `frac_verified * nsite`), or as the expected number
  of verified occasions, e.g. clips listened to (`"occasions"`, times
  `nver`), which is closer to the cost of the design. In a dynamic
  calibration, the fraction is of site-seasons, and the counts are
  multiplied by the number of seasons.

- target:

  Optional coverage to mark with a second line, for instance the lowest
  coverage you are willing to accept (0.90).

- ...:

  Unused.

## Value

A 'ggplot2' object, which can be modified with the usual `+` syntax.

## Examples

``` r
# \donttest{
cal <- fp_verification_curve(p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
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
plot(cal)

plot(cal, type = "bias")

# }
```
