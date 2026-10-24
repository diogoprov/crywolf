# Simulation-based calibration

Checks that the sampler recovers the posterior correctly for a given
design and set of priors (Talts et al. 2018; Modrák et al. 2025). Each
replicate draws parameter values from the priors, simulates data with
the design of `data` (as
[`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)
does), fits the model, and records the rank of the true value among the
posterior draws. If the computation is correct, the ranks are uniform.
Besides the parameters, it ranks the log-likelihood of the simulated
data at the true values among its values at the posterior draws
(`log_lik`), a test quantity that depends on the data and is sensitive
to problems the parameter ranks can miss (Modrák et al. 2025).

## Usage

``` r
fp_sbc(
  data,
  occ = ~1,
  det = ~1,
  fp = ~1,
  ver = ~1,
  markov = NULL,
  random_recorder = !is.null(data$recorder),
  random_fp = FALSE,
  random_det = FALSE,
  priors = fp_priors(),
  reps = 100,
  ndraws = 99,
  base_seed = 1L,
  parallel = FALSE,
  fit_args = list(chains = 4, parallel_chains = 1, iter_warmup = 1000, iter_sampling =
    1000, adapt_delta = 0.95)
)

# S3 method for class 'fp_sbc'
plot(x, bins = 10, ...)
```

## Arguments

- data:

  An `fp_data` object that gives the design.

- occ, det, fp, ver, markov, random_recorder, random_fp, random_det,
  priors:

  As in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  The priors are both the ones the data are drawn from and the ones the
  model is fitted with, as the method requires.

- reps:

  Number of replicates.

- ndraws:

  Number of posterior draws, evenly thinned, among which the rank is
  computed; thinning reduces autocorrelation. It must not exceed the
  number of draws of each fit, and `ndraws + 1` should be a multiple of
  the number of bins used by
  [`base::summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) (10 by
  default).

- base_seed, parallel, fit_args:

  As in
  [`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md).

- x:

  An `fp_sbc` object.

- bins:

  Number of histogram bins.

- ...:

  Unused.

## Value

A data frame of class `fp_sbc` with one row per replicate and parameter:
`rep`, `parameter`, `truth`, `rank` (0 to `ndraws`), `divergences`, and
`ok`. [`base::summary()`](https://rdrr.io/r/base/summary.html) tests the
uniformity of the ranks of each parameter, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws their
histograms with the band expected under uniformity.

## Details

This differs from
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md):
there, data come from fixed parameter values and the question is whether
a model's intervals cover them, which tests the model against the
data-generating process (including misspecification). Here, data come
from the priors of the model itself, so the question is only whether the
computation is right: a sampler problem, an error in the model code, or
a prior that leads to extreme data shows up as ranks that are not
uniform. A U-shaped histogram means posteriors that are too narrow; a
hump, too wide; a slope, a bias.

Prior draws can produce extreme data (for instance, no detection at
all); such replicates are fitted as they are, and their divergent
transitions are recorded. The `log_lik` quantity is not computed with
`random_fp = TRUE` or `random_det = TRUE`, whose site effects are not
kept by the simulation.
[`base::summary()`](https://rdrr.io/r/base/summary.html) uses a
chi-squared test of the binned ranks, which is coarse; look at the
histograms from [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
as well.

## References

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., & Gelman, A.
(2018). Validating Bayesian inference algorithms with simulation-based
calibration. arXiv:1804.06788.

Modrák, M., Moon, A. H., Kim, S., Bürkner, P., Huurre, N., Faltejsková,
K., Gelman, A., & Vehtari, A. (2025). Simulation-based calibration
checking for Bayesian computation: the choice of test quantities shapes
sensitivity. *Bayesian Analysis*, 20, 461–488.
[doi:10.1214/23-BA1404](https://doi.org/10.1214/23-BA1404)

## See also

[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md),
[`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
sbc <- fp_sbc(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
              reps = 100, parallel = TRUE)
summary(sbc)
plot(sbc)
} # }
```
