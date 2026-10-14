# Plot the results of a calibration check

Plots the coverage (or bias) of one scored quantity against the column
of the design that varies, one colour per model. Coverage is shown with
an interval of two Monte Carlo standard errors and a dashed line at the
nominal level; bias, with a dashed line at zero.

## Usage

``` r
# S3 method for class 'fp_calibration'
plot(x, parameter = "occ_x", type = c("coverage", "bias"), by = NULL, ...)
```

## Arguments

- x:

  An `fp_calibration` object from
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  or
  [`fp_verification_curve()`](https://provetelab.org/crywolf/reference/fp_verification_curve.md).

- parameter:

  Scored quantity to plot, by default the occupancy slope (`"occ_x"`).
  See
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
  for the others.

- type:

  `"coverage"` (default) or `"bias"`.

- by:

  Name of the design column for the horizontal axis. By default, the
  first design column that varies among cells; further varying columns
  become panels. With no varying column, the models are on the axis.

- ...:

  Unused.

## Value

A 'ggplot2' object, which can be modified with the usual `+` syntax.

## Examples

``` r
if (FALSE) { # \dontrun{
cal <- fp_verification_curve(p10 = 0.05, conf_effect = 1.5, conf_cor = 0.7,
                             models = list(constant = ~ 1, fp_on_x = ~ x),
                             reps = 50)
plot(cal)
plot(cal, type = "bias")
} # }
```
