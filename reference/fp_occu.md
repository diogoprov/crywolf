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
  priors = fp_priors(),
  chains = 4,
  parallel_chains = chains,
  threads_per_chain = 1,
  iter_warmup = 1000,
  iter_sampling = 1000,
  adapt_delta = 0.95,
  seed = NULL,
  refresh = NULL,
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

- priors:

  Prior settings from
  [`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md).

- chains, iter_warmup, iter_sampling, adapt_delta, seed:

  Passed to the `$sample()` method of 'cmdstanr'.

- parallel_chains:

  Number of chains run in parallel.

- threads_per_chain:

  Number of threads each chain uses to compute the likelihood, which is
  summed over sites in parallel with Stan's `reduce_sum()`. The default,
  1, runs each chain on one core. Larger values speed up fits with many
  sites or occasions (on the `curiango` data, 4 threads made each chain
  3.4 times faster); keep `parallel_chains * threads_per_chain` at or
  below the number of cores. `"auto"` uses the physical cores divided
  among the parallel chains of this fit (not among other processes, such
  as calibration workers). With more than one thread the partial sums
  are added in an order that varies from run to run, so a fit is not
  reproduced exactly by its `seed`; the posterior is the same. Keep 1
  when exact reproduction matters.

- refresh:

  How often each chain reports its progress, in iterations. The default,
  `NULL`, reports every 10% of the iterations in interactive sessions
  and nothing otherwise; `0` silences it. Warm-up iterations are usually
  slower than sampling ones, so the first half of the run takes more
  than half of the time.

- ...:

  Further arguments passed to `$sample()`.

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
probability averaged over the sites).

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

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_priors()`](https://provetelab.org/crywolf/reference/fp_priors.md),
[`fp_calibrate()`](https://provetelab.org/crywolf/reference/fp_calibrate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                   conf_cor = 0.7, seed = 42)

# constant false-positive rate (the simulated recorders have no effect,
# so the recorder term is left out)
m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
# false-positive rate that tracks the occupancy covariate
m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)

coef(m0)["occ_x"]; coef(m1)["occ_x"]   # truth: -0.8
loo::loo_compare(loo(m0), loo(m1))
} # }
```
