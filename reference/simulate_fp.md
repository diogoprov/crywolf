# Simulate data for a false-positive occupancy model

Generates one dataset from known parameters under the multiple-methods
design: an automated method that returns true detections and false
positives at every site, and a verified method, free of false positives,
at a random subset of sites. The false-positive rate can be homogeneous,
vary at random among sites (unstructured), or track a latent gradient
correlated with the occupancy covariate (structured confusion, as when a
soundalike species becomes commoner along the same gradient as the
target).

## Usage

``` r
simulate_fp(
  nsite = 300,
  nocc = 20,
  psi = 0.5,
  b1 = -0.8,
  p11 = 0.6,
  p10 = 0.05,
  sigma_p10 = 0,
  p10_conc = NULL,
  conf_effect = 0,
  conf_cor = 0,
  r11 = 0.6,
  frac_verified = 0.3,
  nver = 3,
  a_bio = 0,
  a_aci = 0,
  sd_recorder = 0,
  nrec = 5,
  seed = NULL
)
```

## Arguments

- nsite:

  Number of sites.

- nocc:

  Number of automated occasions per site.

- psi:

  Occupancy probability at the mean of the covariate.

- b1:

  Occupancy slope on the standardized covariate `x`.

- p11:

  True-detection probability of the automated method at the mean of the
  detection covariates.

- p10:

  False-positive probability of the automated method (its base value;
  the realized mean over sites is returned as `truth$p10_marginal`).

- sigma_p10:

  Among-site standard deviation of the logit of the false-positive rate
  (unstructured, logit-normal heterogeneity). 0 for none.

- p10_conc:

  Optional concentration of a Beta distribution for the site-level
  false-positive rate, with mean `p10` (unstructured, aggregated
  heterogeneity; small values push rates towards 0 or 1). Overrides
  `sigma_p10`.

- conf_effect:

  Slope of the logit of the false-positive rate on the confusing
  gradient `w` (structured confusion). 0 for none.

- conf_cor:

  Correlation between `w` and the occupancy covariate `x`.

- r11:

  Detection probability of the verified method.

- frac_verified:

  Probability that a site is verified.

- nver:

  Number of verified occasions per verified site.

- a_bio, a_aci:

  Effects of two standardized observation covariates (`bio`, `aci`) on
  true detection. 0 by default.

- sd_recorder:

  Standard deviation of a recorder random intercept on true detection. 0
  by default.

- nrec:

  Number of recorders.

- seed:

  Optional integer seed.

## Value

A list of class `fp_sim` with elements `data` (an `fp_data` object whose
`sites` table holds `x` and `w`, and whose `auto` table holds `bio`,
`aci`, and `recorder`) and `truth` (the generating values, the realized
mean false-positive rate `p10_marginal`, the realized correlation
between `w` and `x`, and the latent states `z`).

## Details

The occupancy covariate `x` is standardized, so `b1` is the true slope
exactly. The confusing gradient is \\w = \rho x + \sqrt{1 - \rho^2}
\epsilon\\, with \\\epsilon\\ standard normal, so it has unit variance
and correlation \\\rho\\ with `x`. In practice `w` is not observed; it
is returned in `sites` so that a model given the true structure can
serve as a positive control. The function draws random numbers in the
same order as `simulate_fp()` of the fpcalib research compendium
(v0.1.0), with `nsite`, `nocc`, `frac_verified`, `nver`, `sd_recorder`,
and `nrec` corresponding to its `nsite`, `K`, `verified_coverage`,
`n_verified`, `sig_u`, and `nobser`.

## See also

[`fp_occu()`](https://diogoprov.github.io/crywolf/reference/fp_occu.md),
[`fp_calibrate()`](https://diogoprov.github.io/crywolf/reference/fp_calibrate.md)

## Examples

``` r
# homogeneous false positives
s0 <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
s0
#> <fp_sim> homogeneous false positives
#>   truth: psi = 0.50, b1 = -0.80, p11 = 0.60, p10 = 0.050 (realized mean 0.050), r11 = 0.60
#>   realized occupancy: 0.49
#> <fp_data>
#>   sites:                  100
#>   automated occasions:    1000 (naive occupancy 0.70)
#>   verified sites:         32 (96 occasions, 12 with a confirmed detection)
#>   recorders:              5
#>   site covariates:        x, w 
#>   observation covariates: bio, aci 

# structured confusion
s1 <- simulate_fp(nsite = 100, nocc = 10, conf_effect = 1.5,
                  conf_cor = 0.7, seed = 1)
s1$truth$p10_marginal
#> [1] 0.1081844
```
