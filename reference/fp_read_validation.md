# Read the verdicts of a clip review

Reads `validation_by_clip.csv`, the clip-by-clip result of the review
step of the 'sheepdog' package, and returns the verdicts in a form that
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
can attach to the occasions. Files written by the sampler up to version
1.2.0 (`validacao_por_clipe.csv`, with column names in Portuguese) are
read as well. This legacy format is read throughout crywolf 2.x and will
no longer be read from crywolf 3.0.0; reading it gives no warning before
then.

## Usage

``` r
fp_read_validation(path)
```

## Arguments

- path:

  Path of `validation_by_clip.csv` (or of the legacy
  `validacao_por_clipe.csv`).

## Value

A data frame with one row per clip: `site`, `species_code`, `file` (the
audio file the clip was cut from), `start` and `end` (in seconds within
it), `confidence`, `verdict` (`TRUE` for a true detection, `FALSE` for a
false one, `NA` for clips not reviewed), and `occasion_id` when the file
has it.

## Details

The format is detected from the columns: `site` marks the format of
'sheepdog' 2.0.0 and later, `local` the legacy one. The columns used are
`site`, `species_code`, `source_audio`, `clip_start`, `clip_end`,
`confidence`, and `verdict` (legacy: `local`, `species_code`,
`audio_origem`, `inicio_clipe`, `fim_clipe`, `confidence`, `veredito`),
and `occasion_id` (`ocasiao_id`) when present. Verdicts `positive` and
`negative` (legacy: `positivo`, `negativo`) become `TRUE` and `FALSE`;
`not_reviewed` and `clip_missing` (`nao_validado`, `clipe_ausente`)
become `NA`, and any other value becomes `NA` with a warning.

The sampler can pick clips among the detections of a species by their
confidence or at random. The models in crywolf treat the checked
occasions as chosen independently of the classifier's score: clips
picked for their high scores hold fewer false positives than a random
sample, and would make the false-positive rate look lower than it is.
Use clips drawn at random among the detections for this purpose. Clips
of whole recordings (`sampling_mode = "random_recordings"`, with no
confidence) are a second, manual method of detection: give them to
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md) as
verified occasions, not to the `validation` argument of
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).

## See also

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
