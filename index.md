# crywolf

Bayesian false-positive occupancy models for data from automated
detectors, with simulation-based checks of whether the correction is
calibrated for your design.

Acoustic classifiers and other automated detectors return false
positives. **crywolf** fits the multiple-methods (site-confirmation)
false-positive occupancy model of Miller et al. (2011) in Stan: the
classifier output may be wrong, and a verified subset of detections is
certain. The false-positive rate can be constant or depend on site
covariates, which matters when false positives track the same
environmental gradient as the species.

## Installation

``` r

install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev",
                                       getOption("repos")))
cmdstanr::install_cmdstan()
# install.packages("remotes")
remotes::install_github("diogoprov/crywolf")
```

## Example

``` r

library(crywolf)

# simulated data in which false positives track a gradient correlated
# with the occupancy covariate
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 2026)

m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1)   # constant false-positive rate
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x)   # rate depends on x

coef(m0, "occ"); coef(m1, "occ")                # true slope: -0.8
loo::loo_compare(loo(m0), loo(m1))

# is a constant rate calibrated for this design, and does verifying more help?
fp_verification_curve(nsite = 300, nocc = 20, conf_effect = 1.5,
                      conf_cor = 0.7,
                      models = list(constant = ~ 1, fp_on_x = ~ x),
                      reps = 50)
```

## Main functions

| Function | Purpose |
|----|----|
| [`fp_data()`](https://diogoprov.github.io/crywolf/reference/fp_data.md), [`fp_data_wide()`](https://diogoprov.github.io/crywolf/reference/fp_data_wide.md) | assemble automated and verified detections and site covariates |
| [`fp_occu()`](https://diogoprov.github.io/crywolf/reference/fp_occu.md) | fit the model; methods `print`, `summary`, `coef`, `predict`, `loo`, `as_draws_df` |
| [`fp_priors()`](https://diogoprov.github.io/crywolf/reference/fp_priors.md) | prior settings |
| [`simulate_fp()`](https://diogoprov.github.io/crywolf/reference/simulate_fp.md) | simulate data under homogeneous, unstructured, or structured false positives |
| [`fp_calibrate()`](https://diogoprov.github.io/crywolf/reference/fp_calibrate.md) | coverage and bias of the credible intervals for a design, by simulation |
| [`fp_verification_curve()`](https://diogoprov.github.io/crywolf/reference/fp_verification_curve.md) | coverage against the fraction of sites verified |

## Citation

If you use crywolf, please cite the package (`citation("crywolf")`) and
Miller et al. (2011).

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
<https://doi.org/10.1890/10-1396.1>

## License

GPL (\>= 3). The wolf silhouette in the logo is *Canis lupus
monstrabilis* by Tracy A. Heath, from PhyloPic, in the public domain
(CC0 1.0).
