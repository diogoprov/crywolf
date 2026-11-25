# Simulate data for a Royle-Nichols model with false positives

Generates one dataset from the model of
[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md):
Poisson abundance with a covariate and optional site effect,
per-individual detection with an optional site effect, false positives
on every automated occasion, verified occasions at a share of the sites,
and, optionally, confirmed automated occasions (observation-confirmation
design).

## Usage

``` r
simulate_fp_rn(
  nsite = 200,
  nocc = 10,
  psi = 0.6,
  b1 = 0.5,
  r = 0.3,
  p10 = 0.05,
  sigma_r = 0,
  sigma_lam = 0,
  rv = 0.5,
  frac_verified = 0.3,
  nver = 3,
  frac_confirmed = 0,
  seed = NULL
)
```

## Arguments

- nsite:

  Number of sites.

- nocc:

  Number of automated occasions per site: one number, or one per site
  (unbalanced designs).

- psi:

  Occupancy probability at the mean of the covariate `x`; the intercept
  of log abundance is \\\log(-\log(1 - \psi))\\.

- b1:

  Slope of log abundance on the standardized covariate `x`.

- r:

  Per-individual detection probability of the automated method.

- p10:

  False-positive probability of the automated method.

- sigma_r:

  Standard deviation of a site effect on the logit of `r`.

- sigma_lam:

  Standard deviation of a site effect on log abundance.

- rv:

  Per-individual detection probability of the verified method.

- frac_verified:

  Probability that a site has verified occasions.

- nver:

  Number of verified occasions per verified site.

- frac_confirmed:

  Probability that an automated occasion is confirmed; 1 confirms every
  occasion, as in a fully annotated survey.

- seed:

  Optional integer seed.

## Value

A list of class `fp_sim` with elements `data` (an `fp_data` object;
`sites` holds `x`) and `truth` (the generating values, the intercepts
`b0` of log abundance and `a0` of the logit of `r`, the mean
true-detection probability at occupied sites `p11`, the abundances `N`,
the states `z` and the site effects).

## Details

On each automated occasion a true observation occurs with probability
\\1 - (1 - r_i)^{N_i}\\ and a false one with probability `p10`,
independently; the occasion is a detection if either occurs. Confirmed
occasions get the columns `n_true` and `n_false` of `auto` (0 or 1, `NA`
elsewhere). A verified occasion is a detection with probability \\1 -
(1 - r^v)^{N_i}\\.

## See also

[`fp_occu_rn()`](https://provetelab.org/crywolf/reference/fp_occu_rn.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
s <- simulate_fp_rn(nsite = 50, nocc = 10, sigma_r = 1, seed = 1)
s
#> <fp_sim> Royle-Nichols false positives
#>   truth: psi = 0.60, b1 = 0.50, p11 = 0.46, p10 = 0.050 (realized mean 0.050), r11 = 0.50
#>   realized occupancy: 0.62
#> <fp_data>
#>   sites:                  50
#>   automated occasions:    500 (naive occupancy 0.76)
#>   verified sites:         13 (39 occasions; verified detections at 6 sites)
#>   sites known occupied:   6
#>   site covariates:        x 
```
