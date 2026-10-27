# Read an export of the sheepdog package

Reads the folder written by
[`sheepdog::export_occupancy_data()`](https://provetelab.org/sheepdog/reference/export_occupancy_data.html),
which follows the `sheepdog-occupancy` schema (version 1.x): every
occasion of a passive acoustic survey, the classifier predictions, the
reviewed clips with their sampling design, and the recordings annotated
by ear. The tables come back in the forms that
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
and [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)
expect.

## Usage

``` r
fp_read_sheepdog(path, species = NULL, keep_top = FALSE)
```

## Arguments

- path:

  Folder of the export, holding `schema.json`.

- species:

  Optional species to keep, matched against `species_code` (the
  classifier label, such as `"Nyctidromus albicollis_Common Pauraque"`)
  and `species` (the short code that sheepdog derives from it). Without
  it, an export of several species is kept whole, and the verdicts and
  annotated recordings of all species are mixed: name the target
  species.

- keep_top:

  Keep reviewed clips drawn with `sampling_mode = "top"`? See Details.

## Value

A list of class `fp_sheepdog` with the data frames `occasions` (`site`,
`occasion_id`, `source_audio`, `recording_date`, `recording_hour`,
`night`, `start`, `duration`, then any occasion covariates),
`detections` (one row per prediction: `site`, `occasion_id`,
`species_code`, `species`, `start`, `confidence`), `scores` (highest
score and number of predictions per occasion and species), `validation`
(one row per reviewed clip, with `start`, `end`, `confidence`,
`sampling_mode`, `inclusion_prob`, `verdict` as `TRUE`, `FALSE` or
`NA`), `verified` (one row per annotated recording, with `y`, 1 or 0,
and `inclusion_prob`), `sites`, the parsed `schema`, and `files` (paths
of the tables found). Tables absent from the export are `NULL`.

## Details

**Which reviewed clips.** The false-positive models of crywolf treat the
checked detections as chosen independently of the classifier's score.
Clips of the `"random_detections"` design are; clips of the `"top"`
design, the highest-scoring detections, are not: they hold fewer false
positives than a random sample and would make the false-positive rate
look lower than it is. By default they are dropped with a warning
(`keep_top = FALSE`); they remain useful to calibrate the classifier
within sheepdog.

**Clips and recordings.** A positive clip proves that the species called
in its occasion; a negative clip says only that this prediction was
false. These go to `validation`, for the observation-confirmation design
through `fp_aggregate(validation = )`. A recording annotated by ear
(`"random_recordings"`) gives the state of the whole occasion; these go
to `verified`, the verified method of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md).

**Censored scores.** Only predictions with a confidence of at least
`schema$min_confidence` were exported, so thresholds below it count
every exported prediction and no more.

Occasions are identified by `site` and `occasion_id`; pass
`by = c("site", "occasion_id")` to
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).
Identifiers, dates and the keys of the calibration (such as `month`) are
read as text. The occasion covariates of the export
(`schema$occasion_covariates`, from schema version 1.1) are extra
columns of `occasions`. The calibration draws
(`calibration_draws.parquet`) and the detection histories sheepdog drew
from them (`occasion_draws.parquet`, schema 1.2), when present, are not
read here; their path is in `files`.

## See also

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md),
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# the tables hold only the target species, so fp_aggregate() needs no
# `species`
x <- fp_read_sheepdog("occupancy_export", species = "compau")
auto <- fp_aggregate(x$occasions, x$detections, threshold = 0.5,
                     by = c("site", "occasion_id"), validation = x$validation)
d <- fp_data(auto, x$sites, x$verified, confirmed = c("n_true", "n_false"))
} # }
```
