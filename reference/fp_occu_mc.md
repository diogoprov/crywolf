# Fit the model to draws of the detection histories

Fits [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)
to each of several versions of the same data, which differ only in which
occasions count as detections, and pools the posterior draws in equal
shares. The versions are typically detection histories drawn from the
posterior of a calibration of the classifier, as
[`sheepdog::history_draws()`](https://provetelab.org/sheepdog/reference/history_draws.html)
produces, so that the uncertainty about which detections are true
reaches the occupancy estimates (Monte Carlo propagation).

## Usage

``` r
fp_occu_mc(data_list, ..., seed = NULL, parallel = FALSE, keep_fits = FALSE)

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
  detection histories.

- ...:

  Arguments passed to
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) for
  every fit (formulas, priors, sampler settings).

- seed:

  Optional integer: fit `i` uses seed `seed + i - 1`.

- parallel:

  Fit the versions in parallel with 'furrr', on the workers set by
  [`future::plan()`](https://future.futureverse.org/reference/plan.html).
  Unless set in `...`, `parallel_chains` and `threads_per_chain` are
  then 1, so that each worker uses one core.

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
posterior summary of each fit), `diagnostics` (one row per fit: number
of draws, divergences, largest R-hat, smallest bulk and tail effective
sample sizes), `fits` (with `keep_fits = TRUE`), `n_fits` (the fits that
succeeded), and `call`.

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
pool. Each fit is checked on its own: the summary reports the largest
R-hat and the smallest effective sample size over the fits, and
`diagnostics` lists them by fit.

**Number of versions.** The cost is that of one fit times the number of
versions. The between-fit variance is estimated from the versions, so it
needs several of them; check that the estimates of interest are stable
as versions are added.

## See also

[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md)

## Examples

``` r
if (FALSE) { # \dontrun{
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
                 threads_per_chain = 1, seed = 1)
summary(mc)
} # }
```
