# Detection histories over the draws of a classifier calibration

Builds one version of the data per posterior draw of the calibration
exported by sheepdog (`calibration_draws.parquet`, the probability that
a prediction with a given score is a true detection, at each site and
species), for
[`fp_occu_mc()`](https://provetelab.org/crywolf/reference/fp_occu_mc.md).
In each version an occasion is a detection when at least one of its
predictions is judged true under that draw.

## Usage

``` r
fp_draws_data(
  x,
  rule = c("sample", "threshold"),
  p_min = 0.95,
  n_draws = NULL,
  seed = NULL,
  occasions = x$occasions,
  sites = x$sites,
  site = "site"
)
```

## Arguments

- x:

  An export read by
  [`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md),
  with calibration draws. Read it with `species` set to the target
  species.

- rule:

  How a prediction is judged true under a draw with probability `p` of
  being true: `"sample"` (default) draws it true with probability `p`,
  so that no threshold is involved; `"threshold"` takes it as true when
  `p >= p_min`, as
  [`sheepdog::history_draws()`](https://provetelab.org/sheepdog/reference/history_draws.html)
  does.

- p_min:

  Threshold of `rule = "threshold"`.

- n_draws:

  Optional number of calibration draws to use (the first ones, in the
  order of their identifiers); all by default.

- seed:

  Optional seed of `rule = "sample"`. The random numbers of the session
  continue afterwards as if the call had drawn none.

- occasions:

  The occasions, one row per `site` and `occasion_id` of the export,
  with the column `site` names (the sites of the model, such as
  site-nights) and any detection covariates. Defaults to `x$occasions`.
  Every version carries a copy of it, so keep only the columns the model
  needs when there are many occasions and draws.

- sites:

  The sites table of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
  with the column `site` names. Defaults to `x$sites`.

- site:

  Column of `occasions` and `sites` that identifies the sites of the
  model.

## Value

A list of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
objects, one per calibration draw, with the attributes `draws` (the draw
identifiers) and `rule`.

## Details

**The rules.** `"sample"` judges each prediction on its own, so an
occasion with many low-scoring predictions is likely to hold at least
one judged true; it assumes that the predictions of an occasion are true
or false independently, given the calibration. `"threshold"` keeps only
the predictions the calibration is nearly sure of, so it gives fewer
detections: on the `curiango` data (custom classifier, minute
occasions), about 8% of the minutes with `"sample"` against about 1%
with `"threshold"`.

When the calibration has keys beyond site and species
(`schema$calibration$keys`, such as `"month"`), a prediction takes the
value of its occasion: each such key is a column of the occasions of the
export and of `calibration_draws.parquet`. In an export of schema
version 1.0, without a `month` column, the month is taken from the
recording date of the occasion.

A prediction takes the calibration value of the grid point nearest its
score, as in sheepdog: grid index
`round((score - grid_min) / grid_step) + 1`, within the grid. A
prediction without calibration draws (no draws for its site, species and
keys, or a key missing for its occasion) cannot be scored: its occasion
is set to `NA` (not surveyed) in every version, with a warning, rather
than counted as a non-detection.

**Histories drawn by sheepdog.** An export of schema 1.2 may hold the
histories sheepdog drew with `sheepdog::occasion_draws()`
(`occasion_draws.parquet`), by the same rules. They are read instead of
drawn again when they were made for the target species with the same
`rule`, the same `seed` (`"sample"`, which then needs `seed`) or `p_min`
(`"threshold"`), at least the draws asked for, and every prediction has
its occasion in `occasions`; the values are the same.

**Calibrations with continuous covariates.** When the calibration has
continuous covariates (`schema$calibration$covariates`, schema 1.3, such
as a standardised temperature of the night), the probability that a
prediction is true depends on its occasion, not only on its score, and
sheepdog draws the histories when it exports; the export holds no
calibration draws. The histories are then read as they are: `rule`,
`p_min` and `seed` left at their defaults take the values of the export,
and other values, or more draws than were drawn, are an error (export
again). Occasions sheepdog could not score (without a value of a
covariate) are treated as not surveyed, with a warning, and occasions
not in `occasions` are left out.

The recordings annotated by ear (`x$verified`) are the verified
occasions of every version: they do not depend on the calibration. The
reviewed clips (`x$validation`) are not used here; they are the data of
the calibration.

Histories built this way hold few false positives, so a false-positive
process fitted to them (the default of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)) is
informed mostly by its priors; compare with `fp = NULL`.

## Experimental

This function is experimental (new in 2.0.0, part of the integration
with 'sheepdog'): its arguments and results may change in a minor
version. See
[crywolf-api](https://provetelab.org/crywolf/reference/crywolf-api.md).

## See also

[`fp_occu_mc()`](https://provetelab.org/crywolf/reference/fp_occu_mc.md),
[`fp_read_sheepdog()`](https://provetelab.org/crywolf/reference/fp_read_sheepdog.md)

## Examples

``` r
if (FALSE) { # \dontrun{
x <- fp_read_sheepdog("occupancy_export", species = "nyc_alb")
versions <- fp_draws_data(x, n_draws = 20, seed = 1)
mc <- fp_occu_mc(versions, occ = ~ 1, fp = NULL, threads_per_chain = 1)
summary(mc)
} # }
```
