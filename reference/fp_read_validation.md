# Read the verdicts of the birdnet-validation-sampler

Reads `validacao_por_clipe.csv`, the clip-by-clip result of the review
step of the birdnet-validation-sampler pipeline, and returns the
verdicts in a form that
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
can attach to the occasions.

## Usage

``` r
fp_read_validation(path)
```

## Arguments

- path:

  Path of `validacao_por_clipe.csv`.

## Value

A data frame with one row per clip: `site` (`local` in the sampler),
`species_code`, `file` (the audio file the clip was cut from), `start`
and `end` (in seconds within it), `confidence`, and `verdict` (`TRUE`
for a true detection, `FALSE` for a false one, `NA` for clips not
reviewed).

## Details

The sampler picks clips among the detections of a species, by default
the ones with the highest confidence. The models in crywolf treat the
checked occasions as chosen independently of the classifier's score:
clips picked for their high scores hold fewer false positives than a
random sample, and would make the false-positive rate look lower than it
is. Use clips drawn at random among the detections for this purpose.

## See also

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)
