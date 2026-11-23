# Multispecies false-positive occupancy model

Fits the false-positive occupancy model of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) to
several species at once, with the coefficients of each species drawn
from distributions of the community, so that rare species borrow
strength from the others. Besides
[`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html),
it works with
[`fp_sites()`](https://provetelab.org/crywolf/reference/fp_sites.md) (by
site and species),
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md),
[simulate()](https://provetelab.org/crywolf/reference/simulate.fp_occu_ms.md)
and
[`fp_richness()`](https://provetelab.org/crywolf/reference/fp_richness.md);
the other methods of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) are
not available.

## Usage

``` r
fp_occu_ms(
  data,
  occ = ~1,
  det = ~1,
  fp = ~1,
  random_recorder = !is.null(data$recorder),
  ver = ~1,
  priors = fp_priors(),
  by_species = character(),
  random_fp = FALSE,
  fp_group = NULL,
  centered = "auto",
  centered_occ = "auto",
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

# S3 method for class 'fp_occu_ms'
print(x, digits = 3, ...)

# S3 method for class 'fp_occu_ms'
summary(object, prob = 0.95, ...)

# S3 method for class 'fp_occu_ms'
as_draws_df(x, ...)
```

## Arguments

- data:

  An `fp_data` object with a species column; see the `species` argument
  of [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).

- occ, det, fp, random_recorder, ver, priors:

  As in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  The standard deviations among species have the half-Normal priors of
  `fp_priors(species_sd = )`; the community means have the priors of the
  intercepts, slopes and gap of
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

- by_species:

  Terms that vary by species besides the intercepts of `occ`, `det` and
  the gap: any of `"det"` (the slopes of `det`), `"fp"` (the slopes of
  `fp`) and `"ver"` (the intercept of the verified method, so that the
  reviewer may miss some species more often than others). The default,
  [`character()`](https://rdrr.io/r/base/character.html), shares them
  among species. A term without slopes (for example `"det"` with
  `det = ~ 1`) adds nothing.

- random_fp:

  Add a random effect on the logit of the false-positive rate for each
  site and species? See Details.

- fp_group:

  Optional name of a column of the sites table that groups the sites;
  adds a random effect on the logit of the false-positive rate for each
  group and species. See Details.

- centered:

  Use the centred parametrization of the detection intercepts, gaps,
  recorder effects and the terms of `by_species`? `TRUE`, `FALSE` or
  `"auto"` (the default). The model is the same; only the sampling
  changes. The centred form suits data with much information on each
  species: many sites, or thousands of occasions per species in passive
  acoustic monitoring, where the non-centred form takes very small steps
  and reaches the maximum tree depth. With few sites the centred form
  samples worse.

- centered_occ:

  Use the centred parametrization of the occupancy coefficients of the
  species? `TRUE`, `FALSE` or `"auto"` (the default). Their posteriors
  narrow with the number of sites, not of occasions. `"auto"` makes
  `centered` and `centered_occ` `TRUE` with at least 20 species and 100
  sites, and `FALSE` otherwise, with a message. The rule comes from
  simulations of 5, 20 and 50 species at 30, 100 and 300 sites with 10
  occasions (see
  [`vignette("multispecies", package = "crywolf")`](https://provetelab.org/crywolf/articles/multispecies.md)):
  with 20 or more species and 100 or more sites, the fits converged with
  both centred, and with 30 sites or 5 species they converged more often
  without. Other designs, for example 10 species, were not simulated.

- chains, parallel_chains, threads_per_chain, iter_warmup,
  iter_sampling, adapt_delta, seed, refresh, ...:

  Passed to the sampler, as in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  Runs of detections (`markov`) and site effects on true detection
  (`random_det`) are not available.

- x, object:

  An `fp_occu_ms` object.

- digits:

  Number of digits to print.

- prob:

  Width of the equal-tailed credible interval.

## Value

An object of class `fp_occu_ms`, a list with `fit`, `draws` (readable
names: see Details), `stan_data`, `formulas`, `design`, `data`,
`by_species`, `random_fp`, `fp_group`, `random_recorder`, `priors`,
`sampler`, `diagnostics` and `call`.

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `summary.fp_occu`, as for
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

## Details

Each species \\s\\ has its own occupancy coefficients \\\beta_s =
\mu\_\beta + \sigma\_\beta z_s\\ (all the terms of `occ`, with \\z_s\\
standard normal), its own true-detection intercept \\\alpha_s = \alpha +
\sigma\_\alpha z_s\\, and its own gap between the true-detection and
false-positive intercepts, \\g_s = \mu_g \exp(\sigma_g z_s)\\, which
keeps the false-positive rate of every species below its true-detection
rate at the reference. By default the slopes of `det`, `fp` and `ver`,
the intercept of the verified method and the recorder effect are shared
by the species; `by_species` gives each species its own slopes of `det`
or `fp`, or its own intercept of the verified method, drawn around the
community value as the intercepts are (\\\alpha\_{sk} = \alpha_k +
\sigma\_{\alpha k} z\_{sk}\\). The false positives of a species are
independent of those of the others: the classifier errs in labelling a
species, not in confusing two species.

Two random effects let the false-positive rate vary beyond its
covariates, on its logit, both with the half-Normal prior of
`fp_priors(fp_site_sd = )`: `fp_group` names a column of the sites that
groups them (for example the physical site, when each "site" of the data
is a site in one month), and adds an effect for each group and species,
\\v\_{gs} \sim N(0, \sigma_v)\\; `random_fp = TRUE` adds an effect for
each site and species, \\t\_{is} \sim N(0, \sigma_t)\\. With sites that
are site-months, `fp_group` lets the false-positive rate of a species
differ among physical sites, and `random_fp` among months of the same
site as well.

The species are those of the data (`fp_data(species = )`), and a species
is never assumed to be absent from the list: there is no data
augmentation for species never detected.

**Output names.** The community means have the names of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
(`occ_*`, `det_*`, `fp_*`, `fp_gap`, the median gap, and `p11_ref`,
`p10_ref` at the community means); the standard deviations among species
are `sd_occ_*`, `sd_det_(Intercept)` and `sd_log_fp_gap`; and the values
of each species carry its name in brackets: `occ_*[species]`,
`det_(Intercept)[species]`, `fp_gap[species]`, `p11_ref[species]` and
`p10_ref[species]`. With `by_species`, the slopes add `sd_det_*`,
`det_*[species]`, `sd_fp_*` and `fp_*[species]`, and the verified method
adds `sd_verified_(Intercept)`, `verified_(Intercept)[species]` and
`p_verified[species]`. The standard deviations of the random effects on
the false-positive rate are `sd_fp_group` (`fp_group`) and `sd_fp_site`
(`random_fp`).

**Time and threads.** Each gradient evaluates the likelihood of every
pair of site and species, so its cost grows with the number of species.
The likelihood of each chain is summed over the pairs of site and
species in parallel, on `threads_per_chain` threads. On simulated data
of 20 species with the design of the `curiango` data (1753 site-nights
of 56 recordings), 8 threads per chain made each gradient 5 times faster
than one thread, and 16 threads 7 times faster. Four chains of 1000
warm-up and 1000 sampling iterations, side by side with 8 threads each,
took about 9.5 hours on 32 cores. The default,
`threads_per_chain = "auto"`, divides the cores among the chains in this
way. See
[`vignette("multispecies", package = "crywolf")`](https://provetelab.org/crywolf/articles/multispecies.md).

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`simulate_fp_ms()`](https://provetelab.org/crywolf/reference/simulate_fp_ms.md),
[`fp_richness()`](https://provetelab.org/crywolf/reference/fp_richness.md)

## Examples

``` r
# \donttest{
s <- simulate_fp_ms(nspecies = 4, nsite = 60, seed = 1)
m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, random_recorder = FALSE,
                chains = 2, threads_per_chain = 1)
#> Parametrization (4 species at 60 sites): centered = FALSE, centered_occ = FALSE; see ?fp_occu_ms.
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 100.7 seconds.
#> Chain 2 finished in 110.1 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 105.4 seconds.
#> Total execution time: 110.2 seconds.
#> 
#> Warning: Largest R-hat is 1.013 (> 1.01); run longer chains.
m
#> Multispecies Bayesian false-positive occupancy model (multiple-methods design)
#> 
#>   occupancy:         ~x by species
#>   true detection:    ~bio (intercept by species)  
#>   false positive:    ~1 (gap by species)   
#>   data:              4 species at 60 sites, 2400 automated and 276 verified occasions
#> 
#>              variable     mean     sd    lower  upper rhat ess_bulk ess_tail
#>       occ_(Intercept) -0.02792 0.3459 -0.69571  0.715 1.00      990      967
#>                 occ_x  0.00629 0.2721 -0.56296  0.574 1.00     1044      790
#>       det_(Intercept) -0.10870 0.3840 -0.91817  0.656 1.00      739      955
#>               det_bio  0.00975 0.0623 -0.10939  0.134 1.00     2966     1414
#>  verified_(Intercept)  0.36560 0.1786  0.01120  0.720 1.01     3422     1350
#>            p_verified  0.58971 0.0429  0.50280  0.673 1.01     3422     1350
#>               p11_ref  0.47379 0.0908  0.28533  0.658 1.00      739      955
#>        fp_(Intercept) -2.94096 0.7629 -4.50989 -1.466 1.00      832     1055
#>                fp_gap  2.83226 0.6600  1.58060  4.279 1.00      885      959
#>               p10_ref  0.06242 0.0443  0.01088  0.187 1.00      832     1055
#>    sd_occ_(Intercept)  0.54855 0.3535  0.03238  1.414 1.01      461      341
#>              sd_occ_x  0.32705 0.3045  0.00911  1.193 1.00      716      868
#>    sd_det_(Intercept)  0.76610 0.3368  0.32544  1.628 1.00     1065     1239
#>         sd_log_fp_gap  0.56797 0.3111  0.19801  1.340 1.00      922      956
#> 
#>   Values of each species in summary().
#> 
#>   Sampler: 0 divergent transitions, max R-hat 1.013, min bulk/tail ESS 461/341
# }
```
