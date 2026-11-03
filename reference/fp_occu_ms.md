# Multispecies false-positive occupancy model

Fits the false-positive occupancy model of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) to
several species at once, with the coefficients of each species drawn
from distributions of the community, so that rare species borrow
strength from the others. **This function is in development** (see
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md)):
its arguments, results and output names may change in a minor version,
and it lacks the methods of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md) other
than [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html) and
[`as_draws_df()`](https://mc-stan.org/posterior/reference/draws_df.html).

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

- chains, parallel_chains, threads_per_chain, iter_warmup,
  iter_sampling, adapt_delta, seed, refresh, ...:

  Passed to the sampler, as in
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).
  Runs of detections (`markov`) and site effects (`random_fp`,
  `random_det`) are not available.

- x, object:

  An `fp_occu_ms` object.

- digits:

  Number of digits to print.

- prob:

  Width of the equal-tailed credible interval.

## Value

An object of class `fp_occu_ms`, a list with `fit`, `draws` (readable
names: see Details), `stan_data`, `formulas`, `design`, `data`,
`random_recorder`, `priors`, `sampler`, `diagnostics` and `call`.

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
rate at the reference. The slopes of `det`, `fp` and `ver`, the
intercept of the verified method and the recorder effect are shared by
the species. The false positives of a species are independent of those
of the others: the classifier errs in labelling a species, not in
confusing two species.

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
`p10_ref[species]`.

## See also

[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# name every species to keep (`schema.json` lists them)
x <- fp_read_sheepdog("occupancy_export",
                      species = c("ela_bic", "lep_fus", "phy_cuv"))
d <- fp_draws_data(x, n_draws = 1)[[1]]
m <- fp_occu_ms(d, fp = ~ 1)
m
} # }
```
