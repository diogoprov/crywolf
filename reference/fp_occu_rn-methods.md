# Methods for Royle-Nichols fits

[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
for fits of
[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md).

## Usage

``` r
# S3 method for class 'fp_occu_rn'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_rn'
summary(object, prob = 0.95, ...)

# S3 method for class 'fp_occu_rn'
as_draws_df(x, ...)
```

## Arguments

- x, object:

  A fit of
  [`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md).

- digits:

  Digits to print.

- ...:

  Not used.

- prob:

  Probability mass of the equal-tailed credible intervals.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly;
[`summary()`](https://rdrr.io/r/base/summary.html) an object of class
`summary.fp_occu` (see
[`summary.fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu-methods.md));
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
the posterior draws with readable names.
