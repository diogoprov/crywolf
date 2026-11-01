# Simulate data for a dynamic false-positive occupancy model

Generates one dataset from known parameters for
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md):
the occupancy state of each site changes between seasons by colonization
and extinction, and within each season it is observed as in
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md),
by an automated method with false positives and, at a random subset of
site-seasons, a verified method without them. **This function is
experimental** (see
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md)).

## Usage

``` r
simulate_fp_dyn(
  nsite = 100,
  nseason = 5,
  nocc = 4,
  psi1 = 0.5,
  b_psi1 = 0,
  col = 0.2,
  b_col = 0,
  ext = 0.2,
  b_ext = 0,
  p11 = 0.6,
  p10 = 0.05,
  r11 = 0.8,
  frac_verified = 0.2,
  nver = 2,
  frac_confirmed = 0,
  seed = NULL
)
```

## Arguments

- nsite:

  Number of sites.

- nseason:

  Number of seasons.

- nocc:

  Number of automated occasions per site and season.

- psi1:

  Occupancy probability in the first season, at the mean of the site
  covariate `x`.

- b_psi1:

  Slope of initial occupancy on `x` (logit scale).

- col:

  Colonization probability at the mean of the season covariate `rain`.

- b_col:

  Slope of colonization on `rain`.

- ext:

  Extinction probability at the mean of `x`.

- b_ext:

  Slope of extinction on `x`.

- p11:

  True-detection probability of the automated method.

- p10:

  False-positive probability of the automated method.

- r11:

  Detection probability of the verified method.

- frac_verified:

  Probability that a site-season is verified.

- nver:

  Number of verified occasions per verified site-season.

- frac_confirmed:

  Probability that an automated occasion is checked under the
  observation-confirmation design, as in
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).
  0 (default) for none.

- seed:

  Optional integer seed.

## Value

A list of class `fp_sim_dyn` with elements `data` (an `fp_data` object
with a `season` column and the season covariate `rain`) and `truth` (the
generating values, their logits named as the output of
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md)
in `pars`, and the latent states `z`, a matrix with sites in rows and
seasons in columns).

## Details

The site covariate `x` and the site-season covariate `rain` are
standardized normal draws; `rain` is returned in `season_covs` of the
data. A transition into season \\t\\ uses the covariates of season
\\t\\, as in
[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md).
With `frac_confirmed > 0`, true and false observations occur
independently on every occasion, and checked occasions get the columns
`n_true` and `n_false` of `auto`.

## See also

[`fp_occu_dyn()`](https://provetelab.org/crywolf/reference/fp_occu_dyn.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
s <- simulate_fp_dyn(nsite = 50, nseason = 4, seed = 1)
s
#> <fp_sim_dyn> 50 sites over 4 seasons
#>   truth: psi1 = 0.50, col = 0.20, ext = 0.20, p11 = 0.60, p10 = 0.050, r11 = 0.80
#>   realized occupancy per season: 0.56 0.50 0.46 0.48 
#> <fp_data>
#>   sites:                  50
#>   automated occasions:    800 (naive occupancy 0.84)
#>   verified sites:         35 (94 occasions; verified detections at 16 sites)
#>   sites known occupied:   16
#>   seasons:                4 (sites with data per season: 50 to 50, median 50)
#>   site covariates:        x 
#>   season covariates:      rain 
colMeans(s$truth$z)   # realized occupancy per season
#> [1] 0.56 0.50 0.46 0.48
```
