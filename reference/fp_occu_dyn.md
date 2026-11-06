# Dynamic (multi-season) false-positive occupancy model

Fits a false-positive occupancy model over several seasons (primary
periods), in which the occupancy state of each site changes between
seasons by colonization and extinction. Within each season the
observation model is that of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md): the
multiple-methods and observation-confirmation designs, verified
detections, recorder and site effects, false-positive covariates and
runs of detections.

## Usage

``` r
fp_occu_dyn(
  data,
  psi1 = ~1,
  col = ~1,
  ext = ~1,
  det = ~1,
  fp = ~1,
  random_recorder = !is.null(data$recorder),
  ver = ~1,
  markov = NULL,
  random_fp = FALSE,
  random_det = FALSE,
  priors = fp_priors(),
  chains = 4,
  parallel_chains = chains,
  threads_per_chain = "auto",
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = NULL,
  refresh = NULL,
  ...
)

# S3 method for class 'fp_occu_dyn'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_dyn'
summary(object, prob = 0.95, ...)

# S3 method for class 'fp_occu_dyn'
as_draws_df(x, ...)
```

## Arguments

- data:

  An `fp_data` object with a season column; see the `season` and
  `season_covs` arguments of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).

- psi1:

  One-sided formula for the occupancy probability in the first season.
  It may use site covariates and covariates of the first season.

- col, ext:

  One-sided formulas for the colonization and extinction probabilities.
  They may use site covariates and season covariates
  (`fp_data(season_covs = )`). All three formulas keep their intercept,
  with the priors of
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)
  for intercepts and slopes.

- det, fp, random_recorder, ver, markov, random_fp, random_det, priors:

  As in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md);
  formulas may also use season covariates. `markov` runs of detections
  do not cross seasons. Site effects (`random_fp`, `random_det`) belong
  to sites and are shared by their seasons.

- chains, parallel_chains, threads_per_chain, iter_warmup,
  iter_sampling, adapt_delta, seed, refresh, ...:

  Passed to the sampler, as in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

- x, object:

  An `fp_occu_dyn` object.

- digits:

  Number of digits to print.

- prob:

  Width of the equal-tailed credible interval.

## Value

An object of class `fp_occu_dyn`, a list with `fit`, `draws` (readable
names: see
[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md)),
`stan_data`, `formulas` (`psi1`, `col`, `ext`, `det`, `fp`, `ver`),
`markov`, `design`, `data`, `random_recorder`, `random_fp`,
`random_det`, `priors`, `sampler`, `diagnostics` and `call`. Its methods
are [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html)
and [`simulate()`](https://rdrr.io/r/stats/simulate.html);
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md)
gives the occupancy of each site-season and
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) checks
the fit.

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.fp_occu`, as for
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
whose table has also a `season` column: the season of the derived
quantities `psi_mean[t]`, `psi_fs[t]` and `turnover[t]`, and empty for
the other rows (see
[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md)).

## Details

The occupancy state of site \\i\\ in the first season is Bernoulli with
probability \\\psi\_{1,i}\\. In each later season \\t\\, an occupied
site stays occupied with probability \\1 - \epsilon\_{it}\\ (extinction
\\\epsilon\\) and an unoccupied one becomes occupied with probability
\\\gamma\_{it}\\ (colonization), on the logit scale linear in the
covariates of `col` and `ext`. The covariates of a transition are those
of the season it enters (row \\t\\ of the site-season table for the
transition from \\t - 1\\ to \\t\\). The states are summed out by the
forward algorithm, so the log-likelihood (`log_lik`) is per site, over
all its seasons.

A site without data in a season contributes to that season only through
the transitions. With seasons that are dates,
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) fills
in the missing dates, so colonization and extinction are probabilities
per step (per night, for nightly seasons). A single false positive looks
like a colonization followed by an extinction: without verified or
confirmed data in the seasons, the two processes are hard to separate
from false positives. In simulations with 100 sites (see
[`vignette("dynamic")`](https://provetelab.org/crywolf/articles/dynamic.md)),
the model was calibrated with \\p\_{11} \ge 0.30\\ and at least 10% of
site-seasons verified. With \\p\_{11} = 0.15\\ and no verification,
colonization and extinction were biased upwards by 0.09 to 0.17 and fits
often did not converge; with 10% verified, the bias stayed above 0.05
with fewer than 10 seasons. With site effects on detection of SD 1,
`random_det` needed 2000 + 2000 iterations for R-hat below 1.01 in
nearly all fits.
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
does not simulate dynamic designs.

The fit also gives, per season, the expected and the realized occupancy
over sites (`psi_mean`, `psi_fs`) and the turnover, and, with constant
transitions, the equilibrium occupancy (`psi_eq`); see
[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md).
[`summary()`](https://rdrr.io/r/base/summary.html) lists them with their
season; [`print()`](https://rdrr.io/r/base/print.html) shows only the
parameters.

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`simulate_fp_dyn()`](https://provetelab.org/crywolf/reference/simulate_fp_dyn.md),
[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md);
[`vignette("dynamic")`](https://provetelab.org/crywolf/articles/dynamic.md)
fits it to the `curiango` data, with the nights of each site as seasons.

## Examples

``` r
if (FALSE) { # \dontrun{
s <- simulate_fp_dyn(nsite = 150, nseason = 4, b_col = 1, b_ext = -1, seed = 1)
m <- fp_occu_dyn(s$data, psi1 = ~ x, col = ~ rain, ext = ~ x,
                 random_recorder = FALSE)
m
head(fp_sites(m))   # one row per site-season
} # }
```
