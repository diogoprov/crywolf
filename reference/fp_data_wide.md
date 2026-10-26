# Assemble data from detection matrices

Convenience wrapper around
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) for
data stored as site-by-occasion matrices, the layout used by 'unmarked'.
Missing values (`NA`) mark occasions that were not sampled, so sites may
differ in effort.

## Usage

``` r
fp_data_wide(
  y_auto,
  y_verified = NULL,
  site_covs = NULL,
  obs_covs = NULL,
  recorder = NULL
)
```

## Arguments

- y_auto:

  Matrix of automated detections (sites in rows, occasions in columns),
  coded 0/1/`NA`.

- y_verified:

  Optional matrix of verified detections with the same number of rows;
  unverified sites are rows of `NA`.

- site_covs:

  Optional data frame of site covariates, one row per site, in the same
  order as the rows of `y_auto`.

- obs_covs:

  Optional named list of matrices, each the same size as `y_auto`,
  holding observation covariates of the automated method.

- recorder:

  Optional vector, one value per site, identifying the recorder or
  observer.

## Value

An object of class `fp_data`, as from
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)

## Examples

``` r
set.seed(1)
ya <- matrix(rbinom(40, 1, 0.3), 10, 4)
yv <- matrix(NA, 10, 2); yv[1:3, ] <- rbinom(6, 1, 0.5)
fp_data_wide(ya, yv, site_covs = data.frame(x = rnorm(10)))
#> <fp_data>
#>   sites:                  10
#>   automated occasions:    40 (naive occupancy 0.80)
#>   verified sites:         3 (6 occasions; verified detections at 3 sites)
#>   sites known occupied:   3
#>   site covariates:        x 
```
