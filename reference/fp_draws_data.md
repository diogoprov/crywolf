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

When the calibration depends on the month (`"month"` among
`schema$calibration$keys`), a prediction takes the value of the month of
its recording, from the recording date of its occasion in the export.

A prediction takes the calibration value of the grid point nearest its
score, as in sheepdog: grid index
`round((score - grid_min) / grid_step) + 1`, within the grid. A
prediction without calibration draws (no draws for its site and species,
or month) cannot be scored: its occasion is set to `NA` (not surveyed)
in every version, with a warning, rather than counted as a
non-detection.

The recordings annotated by ear (`x$verified`) are the verified
occasions of every version: they do not depend on the calibration. The
reviewed clips (`x$validation`) are not used here; they are the data of
the calibration.

Histories built this way hold few false positives, so a false-positive
process fitted to them (the default of
[`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md)) is
informed mostly by its priors; compare with `fp = NULL`.

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
