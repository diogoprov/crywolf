# Simulate data for a multispecies false-positive occupancy model

Generates one dataset from known parameters for
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md):
several species on the same sites, occasions, recorders and verified
occasions, each with its own occupancy, true detection and false
positives drawn from distributions of the community, as in the model.
**In development**, as
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).

## Usage

``` r
simulate_fp_ms(
  nspecies = 5,
  nsite = 100,
  nocc = 10,
  psi = 0.4,
  sd_psi = 1,
  b1 = -0.5,
  sd_b1 = 0.5,
  p11 = 0.5,
  sd_p11 = 0.5,
  p10 = 0.05,
  sd_gap = 0.3,
  a_bio = 0,
  r11 = 0.6,
  frac_verified = 0.3,
  nver = 3,
  sd_recorder = 0,
  nrec = 5,
  frac_confirmed = 0,
  ngroup = NULL,
  sd_fp_group = 0,
  sd_fp_site = 0,
  psi_rare = NULL,
  seed = NULL
)
```

## Arguments

- nspecies:

  Number of species.

- nsite, nocc, frac_verified, nver, nrec, frac_confirmed:

  As in
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).

- psi, b1:

  Community mean occupancy probability at the mean of `x`, and community
  mean slope on `x`.

- sd_psi, sd_b1:

  Standard deviations among species of the logit of occupancy at the
  mean of `x` and of the slope.

- p11, p10:

  True-detection and false-positive probabilities of the community (at
  the mean detection intercept and the median gap); `p10` must be below
  `p11`.

- sd_p11:

  Standard deviation among species of the true-detection intercept
  (logit scale).

- sd_gap:

  Standard deviation among species of the log of the gap.

- a_bio:

  Shared slope of true detection on `bio`.

- r11:

  Detection probability of the verified method, below 1.

- sd_recorder:

  Standard deviation of the recorder effects on true detection, shared
  by the species. 0 by default.

- ngroup:

  Optional number of groups of sites (column `group` of the sites,
  assigned at random in equal numbers), as `fp_occu_ms(fp_group = )`
  expects.

- sd_fp_group, sd_fp_site:

  Standard deviations of random effects on the logit of the
  false-positive rate for each group and species (needs `ngroup`) and
  for each site and species, as in
  [`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).
  0 by default.

- psi_rare:

  Optional occupancy probability of the first species at the mean of
  `x`, in place of its draw from the community, to make it rare (its
  slope is still drawn).

- seed:

  Optional integer seed.

## Value

A list of class `fp_sim_ms` with `data` (an `fp_data` object with the
species column `species`; species are named `sp1`, `sp2`, ...) and
`truth`: `pars`, the generating values named as the output of
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
with `occ = ~ x` and `det = ~ bio` (community values and the values of
each species, `occ_x[sp1]`, ...); `z`, the latent states, a matrix of
sites by species; and `p10_bar`, the false-positive rate of each
species.

## Details

The sites, the occupancy covariate `x`, the occasions with their
covariates `bio` and `aci`, the recorders, and which occasions are
verified (or checked, with `frac_confirmed > 0`) come from
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
and are shared by the species. Each species \\s\\ then has occupancy
intercept and slope on `x` \\\beta_s = \mu\_\beta + \sigma\_\beta z_s\\,
true-detection intercept \\\alpha_s = \alpha + \sigma\_\alpha z_s\\, and
gap \\g_s = \mu_g \exp(\sigma_g z_s)\\ between its true-detection and
false-positive intercepts, with \\z_s\\ standard normal, so that its
false-positive rate is \\\mathrm{logit}^{-1} (\alpha_s - g_s)\\. The
slope of true detection on `bio`, the verified detection probability and
the recorder effects are shared. Given these values, the detections of
each species are drawn as in
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md).

## See also

[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)

## Examples

``` r
s <- simulate_fp_ms(nspecies = 4, nsite = 50, nocc = 8, seed = 1)
s
#> <fp_sim_ms> 4 species at 50 sites
#>  species occupancy    p10
#>      sp1      0.28 0.0459
#>      sp2      0.66 0.0450
#>      sp3      0.10 0.0893
#>      sp4      0.58 0.0200
#> <fp_data>
#>   sites:                  50
#>   automated occasions:    1600 (naive occupancy 0.62)
#>   species:                4
#>   verified sites:         13 (156 occasions; verified detections at 8 sites)
#>   site x species known occupied: 13
#>   recorders:              5
#>   site covariates:        x, w 
#>   observation covariates: bio, aci 
#>   by species:
#>  species occasions detections sites_detected naive_occ verified verified_det
#>      sp1       400         70             24      0.48       39            4
#>      sp2       400        101             35      0.70       39           14
#>      sp3       400         55             32      0.64       39            2
#>      sp4       400        104             33      0.66       39            4
```
