# Fit the model to draws of the detection histories

Fits [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
or
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
when the data hold several species, to each of several versions of the
same data, which differ only in which occasions count as detections, and
pools the posterior draws in equal shares. The versions are typically
detection histories drawn from the posterior of a calibration of the
classifier, as
[`sheepdog::history_draws()`](https://provetelab.org/sheepdog/reference/history_draws.html)
produces, so that the uncertainty about which detections are true
reaches the occupancy estimates (Monte Carlo propagation).

## Usage

``` r
fp_occu_mc(
  data_list,
  ...,
  seed = NULL,
  parallel = FALSE,
  max_rhat = 1.05,
  keep_fits = FALSE
)

# S3 method for class 'fp_occu_mc'
summary(object, prob = 0.95, ...)

# S3 method for class 'summary.fp_occu_mc'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_mc'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_mc'
as_draws_df(x, ...)
```

## Arguments

- data_list:

  A list of `fp_data` objects with the same sites, one per draw of the
  detection histories, as
  [`fp_draws_data()`](https://provetelab.org/crywolf/reference/fp_draws_data.md)
  returns. When they hold several species (`fp_data(species = )`), every
  version must have the same species and each is fitted with
  [`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md).

- ...:

  Arguments passed to
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) or
  [`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
  for every fit (formulas, priors, sampler settings).

- seed:

  Optional integer: fit `i` uses seed `seed + i - 1`. The fits are
  reproduced draw for draw only with one thread per chain
  (`threads_per_chain = 1`, the setting when `parallel = TRUE`).

- parallel:

  Fit the versions in parallel with 'furrr', on the workers set by
  [`future::plan()`](https://future.futureverse.org/reference/plan.html).
  Unless set in `...`, `parallel_chains` and `threads_per_chain` are
  then 1, so that each worker uses one core.

- max_rhat:

  Fits whose largest R-hat exceeds this value are left out of the pooled
  draws, with a warning that names them. `Inf` keeps every fit.

- keep_fits:

  Keep the
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) fit
  of every version in the result (`fits`)? They hold the backend objects
  and can be large. Their draws are in memory, but with
  `parallel = TRUE` the CmdStan output files of a worker may be gone
  once it ends.

- object, x:

  An `fp_occu_mc` object.

- prob:

  Width of the equal-tailed intervals.

- digits:

  Digits to print.

## Value

An object of class `fp_occu_mc`: `draws` (the pooled draws, in the
readable names of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), the
chains of fit `i` following those of fit `i - 1`), `summaries` (the
posterior summary of each fit), `diagnostics` (one row per fit that ran:
number of draws, divergences, largest R-hat, smallest bulk and tail
effective sample sizes, and `kept`), `fits` (with `keep_fits = TRUE`,
the fits kept), `n_fits` (the fits in the pool), and `call`.

[`summary()`](https://rdrr.io/r/base/summary.html) returns, for each
quantity, the pooled posterior mean, standard deviation and interval,
the standard deviations within and between fits (`sd_within`,
`sd_between`), the share of the variance between fits (`share_between`),
and the largest R-hat and smallest bulk effective sample size over the
fits.

## Details

**What the pooled posterior is.** Each fit conditions on one version of
the detection histories; pooling the fits in equal shares averages the
posterior over the versions. When the versions are draws from the
posterior of a calibration, this approximates the posterior that
accounts for the uncertainty of the calibration, but the calibration
itself is not updated by the occupancy model (a "cut" model):
information about the state of a site does not flow back to how likely
its detections are to be true.

**Between and within.**
[`summary()`](https://rdrr.io/r/base/summary.html) splits the posterior
variance of each quantity into the part within fits (uncertainty given
one version of the histories) and the part between fits (uncertainty
about the histories), by the law of total variance: the squares of
`sd_within` and `sd_between` add up to the variance of the pooled draws.
A large share between fits means that the result depends on which
detections are taken as true. The between part includes the Monte Carlo
error of each fit's mean, which is small next to it unless the fits are
short.

**Failed fits.** A version whose fit stops with an error is left out
with a warning; the versions should weigh equally, so the function also
warns when the fits have different numbers of draws.

**Convergence.** R-hat over the pooled draws would mix differences
between fits with lack of convergence, so it is not reported for the
pool. Each fit is checked on its own: a fit whose largest R-hat exceeds
`max_rhat` (1.05 by default) has not converged and is left out, with a
warning; report how many were left out. This is looser than the 1.01
that [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
warns about, because mixing error within one fit mostly widens the
within-fit part of the pooled variance; fits between the two values are
kept, and their warnings show in `diagnostics`. Divergent transitions do
not exclude a fit. Leaving fits out conditions the propagation on the
versions that converged: if non-convergence depends on the version (for
instance, histories with very few detections), the pooled posterior is
biased and the between-fit variance understated. Rather than dropping
them, refit those versions with longer chains or a higher `adapt_delta`.
The summary reports the largest R-hat and the smallest effective sample
size over the fits kept, and `diagnostics` lists every fit, with `kept`
marking those in the pool.

**Several species.** The versions are fitted with
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
and the pooled draws hold its parameters (species effects, community
means and standard deviations among species). The parametrization
(`centered`, `centered_occ`) is the same in every fit: when left at
`"auto"`, it is chosen once from the number of species and sites. With
`keep_fits = TRUE`,
[`fp_richness()`](https://provetelab.org/crywolf/reference/fp_richness.md)
and the other methods of
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md)
apply to each fit in `fits`; they are not defined for the pool.

**Number of versions.** The cost is that of one fit times the number of
versions. The between-fit variance is estimated from the versions, so it
needs several of them; check that the estimates of interest are stable
as versions are added.

## Experimental

This function is experimental (new in 2.0.0, part of the integration
with 'sheepdog'): its arguments, its class and methods, and its results
may change in a minor version. See
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md).

## See also

[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`fp_occu_ms()`](https://provetelab.org/crywolf/reference/fp_occu_ms.md),
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md),
[`fp_draws_data()`](https://provetelab.org/crywolf/reference/fp_draws_data.md)

## Examples

``` r
# \donttest{
# three versions of a simulated dataset, with detections switched off at
# random as if they had been judged false
sim <- simulate_fp(nsite = 100, nocc = 8, seed = 1)
versions <- purrr::map(1:3, function(i) {
  d <- sim$data
  set.seed(i)
  hit <- d$auto$y == 1
  d$auto$y[hit] <- stats::rbinom(sum(hit), 1, 0.9)
  d
})
mc <- fp_occu_mc(versions, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                 chains = 2, threads_per_chain = 1, seed = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 2 finished in 9.8 seconds.
#> Chain 1 finished in 10.8 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 10.3 seconds.
#> Total execution time: 10.9 seconds.
#> 
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 2 finished in 9.6 seconds.
#> Chain 1 finished in 9.7 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 9.6 seconds.
#> Total execution time: 9.9 seconds.
#> 
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 9.7 seconds.
#> Chain 2 finished in 11.1 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 10.4 seconds.
#> Total execution time: 11.2 seconds.
#> 
summary(mc)
#> Posterior pooled over 3 fits (95% equal-tailed intervals)
#> 
#>              variable    mean     sd   lower  upper sd_within sd_between
#>       occ_(Intercept) -0.0817 0.2279 -0.5475  0.361    0.2271    0.01886
#>                 occ_x -1.0383 0.2826 -1.6298 -0.516    0.2825    0.00806
#>       det_(Intercept)  0.3471 0.1209  0.1142  0.586    0.1103    0.04959
#>  verified_(Intercept)  0.4101 0.3078 -0.1718  1.030    0.3068    0.02504
#>            p_verified  0.5989 0.0722  0.4572  0.737    0.0720    0.00581
#>               p11_ref  0.5856 0.0292  0.5285  0.643    0.0267    0.01199
#>        fp_(Intercept) -2.5953 0.2414 -3.1008 -2.162    0.2159    0.10801
#>                fp_gap  2.9425 0.2419  2.4821  3.453    0.2329    0.06530
#>               p10_ref  0.0710 0.0156  0.0431  0.103    0.0140    0.00676
#>               p10_bar  0.0710 0.0156  0.0431  0.103    0.0140    0.00676
#>  share_between max_rhat min_ess_bulk
#>       0.006850        1         1674
#>       0.000813        1         1738
#>       0.168119        1         1756
#>       0.006616        1         1642
#>       0.006463        1         1642
#>       0.168151        1         1756
#>       0.200109        1         2109
#>       0.072854        1         1686
#>       0.188876        1         2109
#>       0.188876        1         2109
#> 
#>   Per fit: 0 divergent transitions in all, largest R-hat 1.003, smallest bulk ESS 1642

# several species, from the small sheepdog export shipped with crywolf;
# it has no verified data, so the model without false positives
if (requireNamespace("jsonlite", quietly = TRUE) &&
    requireNamespace("arrow", quietly = TRUE)) {
  dir <- system.file("extdata", "sheepdog-export", package = "crywolf")
  x <- fp_read_sheepdog(dir, species = c("ela_bic", "lep_fus", "phy_cuv"))
  versions <- fp_draws_data(x, n_draws = 3)
  mc_ms <- fp_occu_mc(versions, occ = ~ 1, fp = NULL, random_recorder = FALSE,
                      chains = 2, threads_per_chain = 1, seed = 1)
  summary(mc_ms)
}
#> Warning: 670 reviewed clip(s) drawn with sampling_mode = "top" were dropped: they were chosen for their high scores and would make the false-positive rate look lower than it is. Use keep_top = TRUE to keep them.
#> Warning: 1 species have no automated detection (e.g. lep_fus); their occupancy is informed mostly by the other species.
#> Warning: 3 species have no verified detection (e.g. ela_bic, lep_fus, phy_cuv); their false-positive rates are informed only by the other species and the priors.
#> Parametrization of every version (3 species at 4 sites): centered = FALSE, centered_occ = FALSE; see ?fp_occu_ms.
#> Warning: No verified data: the false-positive rates are identified only by the priors. Treat the estimates with caution.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 3.3 seconds.
#> Chain 2 finished in 3.3 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 3.3 seconds.
#> Total execution time: 3.4 seconds.
#> 
#> Warning: 7 of 2000 (0.0%) transitions ended with a divergence.
#> See https://mc-stan.org/misc/warnings for details.
#> Warning: 7 divergent transition(s); consider a higher `adapt_delta`.
#> Warning: No verified data: the false-positive rates are identified only by the priors. Treat the estimates with caution.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 7.0 seconds.
#> Chain 2 finished in 7.8 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 7.4 seconds.
#> Total execution time: 7.9 seconds.
#> 
#> Warning: No verified data: the false-positive rates are identified only by the priors. Treat the estimates with caution.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 7.9 seconds.
#> Chain 2 finished in 8.9 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 8.4 seconds.
#> Total execution time: 9.0 seconds.
#> 
#> Posterior pooled over 3 fits (95% equal-tailed intervals)
#> 
#>                  variable   mean     sd    lower upper sd_within sd_between
#>           occ_(Intercept)  1.477 1.1901 -0.83616 3.868    1.0365    0.58509
#>           det_(Intercept) -0.290 1.0121 -2.35841 1.192    0.7649    0.66296
#>                   p11_ref  0.447 0.2100  0.08640 0.767    0.1508    0.14614
#>        sd_occ_(Intercept)  1.131 0.7498  0.04868 2.739    0.6419    0.38780
#>  occ_(Intercept)[ela_bic]  2.322 1.3216  0.12848 5.300    1.3200    0.06966
#>  occ_(Intercept)[lep_fus]  0.663 2.1808 -3.87667 4.455    1.4998    1.58334
#>  occ_(Intercept)[phy_cuv]  2.327 1.3346  0.14939 5.383    1.3335    0.05989
#>        sd_det_(Intercept)  1.480 0.9110  0.02616 2.998    0.4602    0.78620
#>  det_(Intercept)[ela_bic]  0.647 0.1135  0.42890 0.871    0.1133    0.00658
#>  det_(Intercept)[lep_fus] -3.083 2.7949 -6.74873 1.432    0.8711    2.65550
#>  det_(Intercept)[phy_cuv]  0.665 0.1082  0.45385 0.875    0.1076    0.01165
#>          p11_ref[ela_bic]  0.656 0.0255  0.60561 0.705    0.0255    0.00148
#>          p11_ref[lep_fus]  0.223 0.3090  0.00117 0.807    0.0689    0.30123
#>          p11_ref[phy_cuv]  0.660 0.0242  0.61155 0.706    0.0241    0.00263
#>  share_between max_rhat min_ess_bulk
#>        0.24165     1.01         1211
#>        0.42898     1.00          456
#>        0.48446     1.00          456
#>        0.26742     1.00         1160
#>        0.00278     1.00         1835
#>        0.52707     1.00         1604
#>        0.00201     1.00         1736
#>        0.74476     1.00          493
#>        0.00337     1.00         2040
#>        0.90284     1.00          850
#>        0.01158     1.00         2089
#>        0.00337     1.00         2040
#>        0.95029     1.00          850
#>        0.01176     1.00         2089
#> 
#>   Per fit: 7 divergent transitions in all, largest R-hat 1.006, smallest bulk ESS 456
# }
```
