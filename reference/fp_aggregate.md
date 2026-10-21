# Detections and confirmations per occasion

Turns classifier detections into the automated occasions of
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md): an
occasion is a detection when the classifier found the species in it with
a confidence at or above `threshold`. Occasions come from a table you
provide, one row per occasion (for instance, one per recording), because
the classifier lists only the detections, not the recordings where it
found nothing.

## Usage

``` r
fp_aggregate(
  occasions,
  detections,
  threshold,
  by = "file",
  species = NULL,
  validation = NULL,
  standardize = NULL
)
```

## Arguments

- occasions:

  Data frame with one row per occasion, with the columns named in `by`,
  the site, and any covariates.

- detections:

  Data frame of detections, such as the output of
  [`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md),
  with the columns in `by` and `confidence`.

- threshold:

  Confidence at or above which a detection counts.

- by:

  Columns that identify an occasion in both tables.

- species:

  Optional species to keep, matched against `species_code` and
  `common_name` of `detections` (and `species_code` of `validation`).

- validation:

  Optional verdicts on detections, such as the output of
  [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
  with the columns in `by`, `confidence`, and `verdict`. They become
  counts of confirmed true and false observations for the
  observation-confirmation design (columns `n_true` and `n_false`, `NA`
  on occasions without a reviewed clip).

- standardize:

  Optional names of numeric columns of `occasions` to standardize (mean
  0, standard deviation 1), as the model expects for detection
  covariates. Their means and standard deviations are kept in the
  attribute `"scaling"` of the result. Site covariates live in the sites
  table of
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md);
  standardize them there.

## Value

`occasions` with three or five more columns: `y` (1 for a detection),
`score` (highest confidence in the occasion, 0 without any),
`n_detections` (detections at or above the threshold), and, with
`validation`, `n_true` and `n_false`. Pass it as `auto` to
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md), with
`confirmed = c("n_true", "n_false")` when it holds confirmations.

## Details

The observation-confirmation design assumes that every detection of a
checked occasion was checked. The function warns when an occasion has
more detections at or above the threshold than reviewed clips, and drops
verdicts on clips below the threshold, which are not detections at this
threshold.

## See also

[`fp_read_birdnet()`](https://provetelab.org/crywolf/reference/fp_read_birdnet.md),
[`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
[`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md)

## Examples

``` r
occasions <- data.frame(site = c("A", "A", "B"),
                        file = c("a1.wav", "a2.wav", "b1.wav"))
detections <- data.frame(file = c("a1.wav", "a1.wav", "b1.wav"),
                         species_code = "compau",
                         confidence = c(0.9, 0.3, 0.2))
fp_aggregate(occasions, detections, threshold = 0.5)
#>   site   file n_detections y score
#> 1    A a1.wav            1 1   0.9
#> 2    A a2.wav            0 0   0.0
#> 3    B b1.wav            0 0   0.2
```
