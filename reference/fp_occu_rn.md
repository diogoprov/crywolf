# Fit a Royle-Nichols occupancy model with false positives

Fits the abundance-induced heterogeneity model of Royle and Nichols
(2003) with false positives in Stan. The number of individuals at site
\\i\\ is \\N_i \sim \mathrm{Poisson}(\lambda_i)\\, each detected on an
occasion with probability \\r\\, so the automated method detects the
species with probability \\1 - (1 - r)^{N_i}\\, and the site is occupied
when \\N_i \> 0\\. A site effect on \\r\\ (`random_det`) then measures
variation in detection among sites beyond that due to abundance, which
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
cannot separate from it.

## Usage

``` r
fp_occu_rn(
  data,
  occ = ~1,
  det = ~1,
  fp = ~1,
  random_recorder = !is.null(data$recorder),
  ver = ~1,
  random_fp = FALSE,
  random_det = FALSE,
  random_abund = FALSE,
  n_max = 20,
  priors = fp_priors(),
  chains = 4,
  parallel_chains = chains,
  threads_per_chain = "auto",
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = NULL,
  refresh = NULL,
  prior_only = FALSE,
  ...
)
```

## Arguments

- data:

  An `fp_data` object from
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
  [`fp_data_wide()`](https://provetelab.org/crywolf/reference/fp_data_wide.md),
  or
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md).

- occ:

  One-sided formula for log abundance (\\\log \lambda\\); site
  covariates only. The occupancy probability is \\\psi = 1 -
  e^{-\lambda}\\.

- det:

  One-sided formula for the logit of the per-individual detection \\r\\
  of the automated method (site or occasion covariates).

- fp:

  One-sided formula for the logit of the false-positive rate, or `NULL`
  for the Royle-Nichols model without false positives. Site covariates
  only.

- random_recorder:

  Logical. Add a recorder random intercept to true detection? Defaults
  to `TRUE` when `data` has a recorder column. Data from
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  always have one; with its default `sd_recorder = 0` set this to
  `FALSE`, because a random effect whose true variance is zero causes
  divergent transitions.

- ver:

  One-sided formula for the logit of the per-individual detection of the
  verified method.

- random_fp:

  Logical. Add a site random intercept to the logit of the
  false-positive rate.

- random_det:

  Logical. Add a site random intercept to the logit of the
  per-individual detection \\r\\: its standard deviation (`sd_det_site`)
  is the variation in detection among sites that abundance does not
  explain.

- random_abund:

  Logical. Add a site random intercept to log abundance
  (`sd_abund_site`), for variation in abundance beyond the covariates of
  `occ`.

- n_max:

  Largest abundance summed over. The fit warns when the Poisson mass
  above it is not negligible at some site.

- priors:

  Prior settings from
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  Single numbers (`seed` may be `NULL`), checked before sampling and
  passed to the `$sample()` method of 'cmdstanr', which runs Stan's
  No-U-Turn sampler (Hoffman & Gelman, 2014). `chains` is the number of
  Markov chains (4 by default); `iter_warmup` and `iter_sampling` are
  the warm-up and kept iterations of each chain (1000 each); warm-up
  draws adapt the sampler and are discarded. `adapt_delta` is the target
  acceptance rate during adaptation (0.95, above the Stan default of
  0.8, because the false-positive models often have curved posteriors);
  raise it towards 0.99 when divergent transitions are reported. `seed`
  is the random seed; with `NULL`, CmdStan draws one, and the seed used
  is stored in `$sampler$seed` either way.

- parallel_chains:

  Number of chains run in parallel.

- threads_per_chain:

  Number of threads each chain uses to compute the likelihood, which is
  summed over sites in parallel with Stan's `reduce_sum()`. One of:

  - `"auto"`: the physical cores divided among the parallel chains of
    this fit (not among other processes, such as calibration workers).
    The fastest choice for most data, but a fit is not reproduced
    exactly by its `seed`: with more than one thread the partial sums
    are added in an order that varies from run to run, so the draws
    drift apart after a while. The posterior is the same.

  - `1`: each chain runs on one core, and the same `seed` gives the same
    draws. Use it when exact reproduction matters (published analyses,
    tests, saved results).

  - a larger number: that many threads per chain, with the same caveat
    as `"auto"`. Keep `parallel_chains * threads_per_chain` at or below
    the number of cores.

  Threads speed up fits with many sites or occasions (on the `curiango`
  data, 4 threads made each chain 3.4 times faster). The default is
  `"auto"` since crywolf 2.0.0 (it was 1 before). Under `R CMD check`
  (when `_R_CHECK_LIMIT_CORES_` is set, as on CRAN) a fit uses at most
  two cores, counting `parallel_chains * threads_per_chain`.

- refresh:

  How often each chain reports its progress, in iterations. The default,
  `NULL`, reports every 10% of the iterations in interactive sessions
  and nothing otherwise; `0` silences it. Warm-up iterations are usually
  slower than sampling ones, so the first half of the run takes more
  than half of the time.

- prior_only:

  Logical. If `TRUE`, the likelihood is left out and the sampler draws
  from the priors alone, through the same Stan program. This checks what
  the priors imply and that the sampler recovers them;
  [`fp_prior_predictive()`](https://provetelab.org/crywolf/reference/fp_prior_predictive.md)
  gives the prior predictive distribution without CmdStan.

- ...:

  Further arguments passed to the `$sample()` method of 'cmdstanr', for
  example `max_treedepth`, `show_messages = FALSE` and
  `show_exceptions = FALSE` (with `refresh = 0`, a silent run),
  `save_warmup = TRUE` (so that `plot(fit, type = "trace")` shows the
  warm-up), or `init`. `init` takes starting values (a number, a list
  per chain, or a function), or a previous fit, `init = old$fit`, so
  that a run starts where an earlier one ended; the draws of the earlier
  run are not pooled with the new ones. Only Stan's No-U-Turn sampler is
  offered.

## Value

An object of class `fp_occu_rn`, with the elements of an
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) fit
(`fit`, `draws`, `stan_data`, `formulas`, `design`, `data`, `priors`,
`sampler`, `diagnostics`, `call`) and `n_max`. The draws hold the
abundance coefficients (`abund_*`), the detection coefficients on the
logit of \\r\\ (`det_*`), the false-positive intercept and slopes, the
verified detection, the derived quantities `lambda_ref` and `psi_ref`
(at the reference of the covariates), `lambda_bar` and `psi_bar` (means
over sites), `r_ref`, `p10_ref` and `r_verified`, and the standard
deviations of the site effects.

## Details

**Model.** With \\N_i \sim \mathrm{Poisson}(\lambda_i)\\, \$\$\log
\lambda_i = X_i \beta \\ \[+ v_i\],\$\$ an automated occasion \\k\\ of
site \\i\\ is a detection with probability \$\$1 - (1 - r\_{ik})^{N_i}
(1 - p\_{10,i}),\$\$ \$\$\mathrm{logit}(r\_{ik}) = Z\_{ik} \alpha +
u\_{\mathrm{rec}\[ik\]} \\ \[+ d_i\],\$\$ \$\$\mathrm{logit}(p\_{10,i})
= \gamma_0 + X^{fp}\_i \gamma \\ \[+ e_i\],\$\$ so false positives occur
on every occasion, also at occupied sites, and an unoccupied site (\\N_i
= 0\\) is detected with probability \\p\_{10,i}\\. A verified occasion
is a detection with probability \\1 - (1 - r^{v})^{N_i}\\. The site
effects \\v_i\\, \\d_i\\ and \\e_i\\ are normal with standard deviations
`sd_abund_site`, `sd_det_site` and `sd_fp_site` and half-normal priors
(`abund_site_sd`, `det_site_sd` and `fp_site_sd` of
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md)).
Unlike
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), the
false-positive intercept \\\gamma_0\\ is not tied to the detection
intercept by a gap: detection at an occupied site always exceeds
\\p\_{10}\\, so the labels cannot switch, and \\p\_{10}\\ may exceed the
per-individual detection \\r\\. Its prior is the marginal prior the gap
implies in
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md):
Normal with mean `-gap_mean` and standard deviation
`sqrt(intercept_sd^2 + gap_sd^2)`. \\N_i\\ is summed out over 0 to
`n_max`.

**Observation-confirmation design.** Confirmed occasions (argument
`confirmed` of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md))
enter as in
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md), with
probabilities of at least one true observation \\s_1 = 1 - (1 -
r)^{N_i}\\ and at least one false observation \\s_0 = p\_{10}\\. They
can be the same occasions as the automated detections, so a fully
annotated survey enters as one row per occasion with both columns
filled.

**Outputs.** [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
[`loo()`](https://mc-stan.org/loo/reference/loo.html) (leave one site
out), [`simulate()`](https://rdrr.io/r/stats/simulate.html) (new
abundances and detections with the design of the data) and
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md) and
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md), which
use the posterior of each site's abundance given its data; see
[fp_occu_rn-methods](https://provetelab.org/crywolf/reference/fp_occu_rn-methods.md).

**Not supported.** False-positive covariates measured on each occasion,
runs of detections (`markov` of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)), and
the other outputs of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
([`predict()`](https://rdrr.io/r/stats/predict.html),
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fp_effects()`](https://provetelab.org/crywolf/reference/fp_effects.md),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)).

**Identifiability.** Simulations at 200 sites and at 21 sites with 17
occasions (in `dev/rn-fp/` of the package repository) recovered
`sd_det_site` near 1 and gave few false signals when it was 0; values
near 0.5 were not distinguishable from 0, and `sd_abund_site` was weakly
identified. With few sites the prior of `sd_det_site` carries much of
the posterior near 0: compare fits under several values of
`det_site_sd`.

## Experimental

This function is experimental (new after 2.3.0): its arguments and
results may change in a minor version. See
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md).

## References

Royle, J. A. and Nichols, J. D. (2003). Estimating abundance from
repeated presence-absence data or point counts. *Ecology*, 84, 777–790.
[doi:10.1890/0012-9658(2003)084\[0777:EAFRPA\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282003%29084%5B0777%3AEAFRPA%5D2.0.CO%3B2)

Chambert, T., Miller, D. A. W. and Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339.
[doi:10.1890/14-1507.1](https://doi.org/10.1890/14-1507.1)

## See also

[`simulate_fp_rn()`](https://provetelab.org/crywolf/reference/simulate_fp_rn.md),
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)

## Examples

``` r
# \donttest{
s <- simulate_fp_rn(nsite = 60, nocc = 10, sigma_r = 1, seed = 1)
m <- fp_occu_rn(s$data, random_det = TRUE, chains = 2,
                iter_warmup = 300, iter_sampling = 300, seed = 1)
#> Running MCMC with 2 parallel chains, with 2 thread(s) per chain...
#> 
#> Chain 1 finished in 4.0 seconds.
#> Chain 2 finished in 4.8 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 4.4 seconds.
#> Total execution time: 4.8 seconds.
#> 
#> Warning: Smallest bulk ESS is 187 (< 400).
m
#> Bayesian Royle-Nichols occupancy model with false positives
#> 
#>   log abundance:    ~1  
#>   per-individual detection:  ~1  + (1 | site) 
#>   false positive:   ~1  
#>   data:             60 sites, 600 automated (0 confirmed) and 42 verified occasions; N summed over 0..20
#> 
#>              variable    mean     sd   lower  upper  rhat ess_bulk ess_tail
#>     abund_(Intercept) -0.0765 0.3054 -0.6281  0.563 1.000      264      172
#>       det_(Intercept) -1.0104 0.3620 -1.7585 -0.352 0.998      247      368
#>  verified_(Intercept)  0.1929 0.5032 -0.7819  1.189 0.999      434      325
#>            r_verified  0.5456 0.1178  0.3139  0.767 0.999      434      325
#>            lambda_ref  0.9716 0.3145  0.5336  1.756 1.000      264      172
#>               psi_ref  0.6045 0.1075  0.4135  0.827 1.000      264      172
#>            lambda_bar  0.9716 0.3145  0.5336  1.756 1.000      264      172
#>               psi_bar  0.6045 0.1075  0.4135  0.827 1.000      264      172
#>                 r_ref  0.2727 0.0686  0.1470  0.413 0.999      247      368
#>        fp_(Intercept) -2.6733 0.5843 -4.0899 -1.965 1.003      267      128
#>               p10_ref  0.0712 0.0268  0.0165  0.123 1.004      267      128
#>           sd_det_site  0.7320 0.3105  0.1833  1.431 1.007      187      441
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.007, min bulk/tail ESS 187/128
# }
```
