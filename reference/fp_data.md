# Assemble data for a false-positive occupancy model

Collects the automated detections, the verified detections, and the site
table into one object that
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) can
fit. All tables are in long format, one row per occasion, which
accommodates unequal numbers of occasions per site without padding. Use
[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md)
if your data are detection matrices, as in 'unmarked'.

## Usage

``` r
fp_data(auto, sites, verified = NULL, site = "site", y = "y", recorder = NULL)
```

## Arguments

- auto:

  Data frame of the automated method (the one that can return false
  positives), one row per site and occasion. It needs a site column and
  a 0/1 detection column, and may hold observation covariates for the
  true-detection process and a recorder (or observer) identifier.

- sites:

  Data frame with one row per site: the site identifier plus any site
  covariates used in the occupancy (`occ`) or false-positive (`fp`)
  formulas of
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  Sites with no automated occasions are allowed and contribute only
  through their verified data, if any.

- verified:

  Optional data frame of the verified method (detections checked by an
  expert and therefore free of false positives), one row per site and
  occasion, with the same site and detection columns. Without it the
  false-positive rate is identified only by the priors, and
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
  warns about it.

- site, y:

  Names of the site and detection columns, shared by all tables.

- recorder:

  Optional name of a column of `auto` that identifies the recorder or
  observer, used for a random intercept on true detection.

## Value

An object of class `fp_data`: a list with the three tables (rows with a
missing detection dropped), the column names, and the vector
`known_present`, which flags sites with at least one verified detection.

## See also

[`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
sim <- simulate_fp(nsite = 50, nocc = 5, seed = 1)
sim$data
#> <fp_data>
#>   sites:                  50
#>   automated occasions:    250 (naive occupancy 0.68)
#>   verified sites:         18 (54 occasions, 10 with a confirmed detection)
#>   recorders:              5
#>   site covariates:        x, w 
#>   observation covariates: bio, aci 

# the same structure built by hand
d <- fp_data(auto = sim$data$auto, sites = sim$data$sites,
             verified = sim$data$verified, recorder = "recorder")
d
#> <fp_data>
#>   sites:                  50
#>   automated occasions:    250 (naive occupancy 0.68)
#>   verified sites:         18 (54 occasions, 10 with a confirmed detection)
#>   recorders:              5
#>   site covariates:        x, w 
#>   observation covariates: bio, aci 
```
