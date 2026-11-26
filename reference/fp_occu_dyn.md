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
# \donttest{
s <- simulate_fp_dyn(nsite = 150, nseason = 4, b_col = 1, b_ext = -1, seed = 1)
m <- fp_occu_dyn(s$data, psi1 = ~ x, col = ~ rain, ext = ~ x,
                 random_recorder = FALSE,
  chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 30.3 seconds.
#> Chain 2 finished in 31.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 30.8 seconds.
#> Total execution time: 31.5 seconds.
#> 
m
#> Dynamic Bayesian false-positive occupancy model (multiple-methods design)
#> 
#>   initial occupancy: ~x 
#>   colonization:      ~rain 
#>   extinction:        ~x 
#>   true detection:    ~1   
#>   false positive:    ~1  
#>   data:              150 sites over 4 seasons, 2400 automated and 216 verified occasions
#> 
#>              variable     mean      sd   lower   upper  rhat ess_bulk ess_tail
#>      psi1_(Intercept)  0.05292 0.18184 -0.3000  0.3954 1.005     2790     1388
#>                psi1_x  0.00795 0.17215 -0.3282  0.3321 1.000     3190     1441
#>       det_(Intercept)  0.53337 0.06903  0.3975  0.6686 1.002     2766     1408
#>  verified_(Intercept)  1.27509 0.24949  0.8080  1.7743 1.000     2525     1264
#>            p_verified  0.77867 0.04232  0.6917  0.8550 1.000     2525     1264
#>               p11_ref  0.63012 0.01607  0.5981  0.6612 1.002     2766     1408
#>        fp_(Intercept) -2.79939 0.15384 -3.1292 -2.4955 0.999     2803     1384
#>                fp_gap  3.33276 0.15460  3.0234  3.6531 1.001     2512     1258
#>               p10_ref  0.05792 0.00833  0.0419  0.0762 0.999     2803     1384
#>               p10_bar  0.05792 0.00833  0.0419  0.0762 0.999     2803     1384
#>       col_(Intercept) -1.55576 0.25887 -2.1048 -1.1059 0.999     1689     1121
#>              col_rain  1.40508 0.30023  0.8588  2.0664 1.003     1608     1179
#>       ext_(Intercept) -1.35923 0.20684 -1.7684 -0.9779 0.999     2542     1554
#>                 ext_x -0.94826 0.23847 -1.4219 -0.4963 1.001     2665     1317
#> 
#>   Occupancy per season (psi_mean, psi_fs), turnover in summary().
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.005, min bulk/tail ESS 1608/1121
head(fp_sites(m))   # one row per site-season
#>   site          x season        rain       psi psi_lower psi_upper     occupied
#> 1    1 -0.7169114      1  0.45018710 0.5116521 0.4030334 0.6202032 0.0152747771
#> 2    1 -0.7169114      2  0.89367370 0.5464184 0.4706894 0.6259788 0.0006151916
#> 3    1 -0.7169114      3  0.96169663 0.5651681 0.4759732 0.6565700 0.0080835775
#> 4    1 -0.7169114      4 -0.34106698 0.4275434 0.3304659 0.5291732 0.0036400833
#> 5    2  0.1790354      1 -0.01855983 0.5134748 0.4256221 0.5977721 0.7445787710
#> 6    2  0.1790354      2 -1.04729815 0.4458094 0.3681717 0.5202678 0.7817167652
#>   occupied_lower occupied_upper detections verified        p10
#> 1   0.0074172531    0.028612973          0    FALSE 0.05792068
#> 2   0.0001668958    0.001527347          0    FALSE 0.05792068
#> 3   0.0038096811    0.014956343          0    FALSE 0.05792068
#> 4   0.0016852553    0.006600718          0    FALSE 0.05792068
#> 5   0.5377077629    0.906142370          1    FALSE 0.05792068
#> 6   0.6065914030    0.915079458          1    FALSE 0.05792068
# }
```
