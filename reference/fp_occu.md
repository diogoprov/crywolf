# Fit a Bayesian false-positive occupancy model

Fits the multiple-methods false-positive occupancy model of Miller et
al. (2011) in Stan. Detections from an automated method may be false
positives, whereas detections from a verified subset are certain. The
occupancy probability depends on site covariates, true detection on
observation (or site) covariates and an optional recorder random
intercept, and the false-positive rate is either constant or a function
of site covariates.

## Usage

``` r
fp_occu(
  data,
  occ = ~1,
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

  One-sided formula for the logit of occupancy, using columns of
  `data$sites`.

- det:

  One-sided formula for the logit of the true-detection probability of
  the automated method. It may use columns of `data$auto` and of
  `data$sites`. It must keep the intercept.

- fp:

  One-sided formula for the logit of the false-positive rate, using
  columns of `data$sites` or, for covariates measured on each occasion
  (background noise, for instance), of `data$auto`; the rate then varies
  among the occasions of a site. `~ 1` (the default) fits a constant
  rate; `~ w` lets the rate change along `w`. Use `NULL` to fit a
  standard occupancy model without false positives. The intercept is
  fixed by the model (see Details) and cannot be removed.

- random_recorder:

  Logical. Add a recorder random intercept to true detection? Defaults
  to `TRUE` when `data` has a recorder column. Data from
  [`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
  always have one; with its default `sd_recorder = 0` set this to
  `FALSE`, because a random effect whose true variance is zero causes
  divergent transitions.

- ver:

  One-sided formula for the logit of the detection probability of the
  verified method, using columns of `data$verified` or `data$sites` (for
  example, who listened, or the length of the clip). The default, `~ 1`,
  keeps it constant.

- markov:

  Optional name of a numeric, date, or date-time column of `data$auto`
  that places the occasions of each site in time (an occasion number or
  a recording time). With it, true detection depends on whether the
  previous occasion of the same site was a detection, through the
  covariate `y_prev` added to `det` (reported as `det_y_prev`): a
  first-order Markov term for runs of detections, such as a bird calling
  for several consecutive minutes. See Details.

- random_fp:

  Logical. Add a site random intercept to the logit of the
  false-positive rate, for rates that vary among sites in ways the
  covariates in `fp` do not capture. `FALSE` by default; see the section
  "Site random effect on the false-positive rate" before turning it on.

- random_det:

  Logical. Add a site random intercept to the logit of the
  true-detection probability of the automated method, for detection that
  varies among occupied sites in ways the covariates in `det` do not
  capture (distance to the recorder, abundance, local noise). `FALSE` by
  default; it can be combined with `random_fp`. See the section "Site
  random effect on true detection".

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

An object of class `fp_occu`, with methods for
[`base::print()`](https://rdrr.io/r/base/print.html),
[`base::summary()`](https://rdrr.io/r/base/summary.html),
[`stats::coef()`](https://rdrr.io/r/stats/coef.html),
[`stats::predict()`](https://rdrr.io/r/stats/predict.html),
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html), and
[`posterior::as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html).
The 'cmdstanr' fit is in element `fit`.

## Details

**Checks before sampling.** Covariates must be complete and finite (a
missing or infinite value stops the fit; impute or drop the rows before
fitting), and the columns of each design matrix must be linearly
independent. A constant or collinear covariate, or an unused factor
level, triggers a warning of class `crywolf_design_warning` that names
the redundant column. The fit still runs, because the priors keep the
posterior proper, but the posterior of that coefficient only repeats its
prior; remove the column, or drop unused levels with
[`droplevels()`](https://rdrr.io/r/base/droplevels.html). Factor
covariates use the contrasts of R's
[`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html):
treatment contrasts for unordered factors and orthogonal polynomials for
ordered ones (columns such as `habitat.L`); set
[`contrasts()`](https://rdrr.io/r/stats/contrasts.html) on the factor to
change them.

**Diagnostics and warnings.** After sampling, the fit warns when there
are divergent transitions, when the largest R-hat exceeds 1.01, when the
smallest bulk effective sample size is below 400, when transitions hit
the maximum tree depth, and when the E-BFMI of a chain is below 0.3.
These warnings have class `crywolf_sampler_warning`, so they can be
muffled selectively, for example with
`withCallingHandlers(fp_occu(...), crywolf_sampler_warning = function(w) invokeRestart("muffleWarning"))`;
the diagnostics are kept in `$diagnostics` either way.

**Identifiability.** The intercept of the false-positive rate equals the
intercept of true detection minus a positive gap, \\\gamma_0 =
\alpha_1 - \Delta\\, so that at the reference values of the covariates
(all zero) a true detection is more likely than a false positive. This
removes the label switching between the occupied and unoccupied states.
Centre the covariates, so that the reference is an average site;
[`simulate_fp()`](https://provetelab.org/crywolf/reference/simulate_fp.md)
already returns them standardized.

**Verified data.** The verified detections anchor the false-positive
rate. With none, the model still runs, but the rate is informed only by
the priors and the gap, and the fit should not be trusted.

**Confusion between species.** The false-positive rate does not depend
on which other species occur at the site. The model assumes that the
classifier mislabels noise, not the calls of another species. When the
false positives of a species come from the calls of a species that is
often present, they follow the occupancy of that species, and the model
reads them as presence. In a small simulation (100 sites, 10 automated
occasions, 30% of sites verified, three datasets per level), occupancy
was not biased when a confused detection had a probability of 0.15 per
occasion. At 0.30, close to the rate of true detection, the occupancy
intercept was overestimated by about 1.7 on the logit scale in two of
the three datasets. Where a classifier is known to confuse two species,
verify detections of each at sites where the other occurs.

**Observation-confirmation design.** When `data` has confirmed occasions
(argument `confirmed` of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)), the
model follows the observation-confirmation design of Chambert et al.
(2015): an expert checks every observation of some automated occasions
and finds only false ones, only true ones, both, or none. False
observations can then occur at occupied sites too, and `det` and `fp`
model \\s_1\\ and \\s_0\\, the probabilities of at least one true and at
least one false observation on an occasion; a site-level detection at an
occupied site has probability \\p\_{11} = 1 - (1 - s_1)(1 - s_0)\\, and
`p11_ref` reports it at the reference. The confirmed occasions identify
\\s_0\\ directly, so verified data are optional; they can be combined
with this design. The gap still applies, so that \\s_1 \> s_0\\ at the
reference.

**Site random effect on the false-positive rate.** With
`random_fp = TRUE` the logit of the false-positive rate of site \\i\\
gains a term \\e_i \sim \mathrm{Normal}(0, \sigma\_{fp})\\, reported as
`sd_fp_site`, with a half-Normal prior (`fp_priors(fp_site_sd = )`). In
simulations with 300 sites, 20 automated occasions, and 30% of sites
verified, it restored the coverage of the mean false-positive rate
(`p10_bar`) when the rate varied among sites at random (from 0.00 and
0.55 with a constant rate to 0.95) and cost nothing when the rate was
homogeneous. Use it with care:

- It adds one parameter per site, and fits took 6 to 9 times longer.

- It mixes less well: in those simulations the largest R-hat exceeded
  1.01 in a quarter to three quarters of the fits. Check the
  diagnostics, and run longer chains (`iter_warmup`, `iter_sampling`)
  when they are poor.

- Each site has a single latent state, so a site effect on false
  positives is told apart from occupancy mainly by the verified data;
  with few occasions or little verification it may be weakly identified.

- It does not fix false positives that follow an environmental gradient.
  Those bias the occupancy effects, and the remedy is to put the
  occupancy covariates, and any measured correlate of the confusion, in
  `fp` (see below). Prefer covariates, and add the random effect for the
  variation that remains.

Predictions with `predict(process = "fp")` are for a site with \\e_i =
0\\, a typical site.

**Site random effect on true detection.** With `random_det = TRUE` the
logit of the true-detection probability of the automated method at site
\\i\\ gains a term \\d_i \sim \mathrm{Normal}(0, \sigma\_{det})\\,
reported as `sd_det_site`, with a half-Normal prior
(`fp_priors(det_site_sd = )`). It is for detection that varies among
occupied sites for reasons the covariates in `det` do not capture, such
as the distance from the recorder to the calling animals, their
abundance, or local noise. The detection probability of the verified
method has no site effect. `p11_ref` is the true-detection probability
of a typical site (\\d_i = 0\\).

The two site effects answer different questions, and one cannot stand in
for the other: in simulations where only detection varied among sites, a
model with `random_fp = TRUE` alone attributed that variation to false
positives and estimated a large `sd_fp_site` although the true rate was
constant. Both can be fitted together. They are then identified by
different data (true detection only at occupied sites, false positives
only at unoccupied ones, with the verified data telling the two apart),
but they mix slowly, and `fp_occu()` warns about it: check R-hat and the
effective sample sizes, and run longer chains when they are poor.
[`vignette("site_effects")`](https://provetelab.org/crywolf/articles/site_effects.md)
shows both effects in simulations and on the `curiango` data.

**Detections that come in runs.** Occasions close in time are not
independent when the species, once present, is active for a while: a
calling bird is recorded in consecutive minutes, which makes detections
more clustered among sites than independent occasions allow (a posterior
predictive check flags it). `markov` adds the detection on the previous
occasion as a covariate of true detection, so that \\p\_{11}\\ can be
higher right after a detection. It acts on the occupied state only:
false positives stay independent among occasions. `y_prev` is the
detection on the occasion one time step earlier, where the step is the
shortest interval between consecutive occasions of a site in the data
(one occasion, or 15 minutes for recordings every 15 minutes). The first
occasion of each site, and an occasion that follows a gap (for instance,
a missing occasion), has `y_prev` = 0: its probability of detection is
the one after no detection, so a run is not carried across a gap. Two
occasions of a site at the same time are an error.
[`fp_ppc()`](https://provetelab.org/crywolf/reference/fp_ppc.md) checks
one step ahead, replicating each occasion given the observed previous
one;
[`simulate.fp_occu()`](https://provetelab.org/crywolf/reference/simulate.fp_occu.md)
simulates the runs sequentially.

**When to model the false-positive rate.** In the simulations behind
this package, a constant rate gave calibrated intervals for the
occupancy slope when false positives were homogeneous; random variation
of the rate among sites made the estimate of the rate unreliable but
largely spared the slope. When the rate tracked an environmental
gradient correlated with the occupancy covariate and the confusion was
strong, the slope was biased and verifying more sites did not fix it.
Placing the occupancy covariate in `fp` recovered part of the
calibration, and placing the gradient itself recovered all of it. Use
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)
to check your own design.

**Parameter names.** Coefficients are labelled by process and column of
the design matrix: `occ_*` (occupancy), `det_*` (true detection), `fp_*`
(false positive, where `fp_(Intercept)` is \\\gamma_0\\), plus `fp_gap`
(\\\Delta\\), `sd_recorder`, `p_verified` (detection probability of the
verified method), `p11_ref` and `p10_ref` (true- and false-positive
probabilities at the reference), and `p10_bar` (false-positive
probability averaged over sites or automated occasions, depending on the
FP covariates). See
[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md)
for all names, scales, and availability conditions.

## References

Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
positive detections in species occurrence data under different study
designs. *Ecology*, 96, 332–339.
[doi:10.1890/14-1507.1](https://doi.org/10.1890/14-1507.1)

Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
two types of observational error occur: non-detection and species
misidentification. *Ecology*, 92, 1422–1428.
[doi:10.1890/10-1396.1](https://doi.org/10.1890/10-1396.1)

Hoffman, M. D., & Gelman, A. (2014). The No-U-Turn sampler: adaptively
setting path lengths in Hamiltonian Monte Carlo. *Journal of Machine
Learning Research*, 15, 1593–1623.
<https://jmlr.org/papers/v15/hoffman14a.html>

## See also

[crywolf-output](https://provetelab.org/crywolf/reference/crywolf-output.md),
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)

## Examples

``` r
# \donttest{
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 42)

# constant false-positive rate (the simulated recorders have no effect,
# so the recorder term is left out)
m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
              chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 29.4 seconds.
#> Chain 2 finished in 32.4 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 30.9 seconds.
#> Total execution time: 32.5 seconds.
#> 
# false-positive rate that tracks the occupancy covariate
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE,
              chains = 2, threads_per_chain = 1)
#> Running MCMC with 2 parallel chains, with 1 thread(s) per chain...
#> 
#> Chain 1 finished in 38.1 seconds.
#> Chain 2 finished in 38.1 seconds.
#> 
#> Both chains finished successfully.
#> Mean chain execution time: 38.1 seconds.
#> Total execution time: 38.3 seconds.
#> 

coef(m0)["occ_x"]; coef(m1)["occ_x"]   # truth: -0.8
#>      occ_x 
#> -0.5030135 
#>      occ_x 
#> -0.6773906 
loo::loo_compare(loo(m0), loo(m1))
#>   model elpd_diff se_diff p_worse diag_diff diag_elpd
#>  model2       0.0     0.0      NA                    
#>  model1     -24.3    11.5    0.98                    
# }
```
