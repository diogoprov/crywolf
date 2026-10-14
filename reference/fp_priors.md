# Prior settings for fp_occu()

Normal priors on the logit scale. The defaults are weakly informative
for standardized covariates.

## Usage

``` r
fp_priors(
  intercept_sd = 1.5,
  slope_sd = 2.5,
  gap_mean = 2,
  gap_sd = 1.5,
  recorder_sd = 2
)
```

## Arguments

- intercept_sd:

  Standard deviation of the Normal(0, sd) priors on the intercepts of
  occupancy, true detection, and verified detection.

- slope_sd:

  Standard deviation of the Normal(0, sd) priors on all slopes,
  including those of the false-positive rate.

- gap_mean, gap_sd:

  Mean and standard deviation of the Normal prior, truncated at zero, on
  the gap between the true-detection and false-positive intercepts. The
  default, Normal(2, 1.5), puts the false-positive rate at the reference
  roughly an order of magnitude below true detection, while leaving room
  for the data.

- recorder_sd:

  Scale of the half-Normal(0, sd) prior on the standard deviation of the
  recorder random intercept.

## Value

A list of class `fp_priors`.

## Examples

``` r
fp_priors()
#> <fp_priors> (logit scale)
#>   intercepts      ~ Normal(0, 1.5)
#>   slopes          ~ Normal(0, 2.5)
#>   fp_gap          ~ Normal(2, 1.5), truncated at 0
#>   sd_recorder     ~ half-Normal(0, 2)
fp_priors(slope_sd = 1)
#> <fp_priors> (logit scale)
#>   intercepts      ~ Normal(0, 1.5)
#>   slopes          ~ Normal(0, 1)
#>   fp_gap          ~ Normal(2, 1.5), truncated at 0
#>   sd_recorder     ~ half-Normal(0, 2)
```
