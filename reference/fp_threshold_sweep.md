# Fit the model at several classifier thresholds

The threshold at which a confidence score counts as a detection trades
false positives against missed detections, and it can change the
estimated environmental effects (Piatti et al. 2026). This function
builds the automated occasions at each threshold with
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md),
fits the same model to each, and collects the estimates, so that their
sensitivity to the threshold can be seen at once.

## Usage

``` r
fp_threshold_sweep(
  occasions,
  detections,
  thresholds = c("MaxF1", "MaxPrecision"),
  sites,
  verified = NULL,
  site = "site",
  by = "file",
  species = NULL,
  recorder = NULL,
  occ = ~1,
  det = ~1,
  fp = ~1,
  random_recorder = !is.null(recorder),
  prob = 0.95,
  validation = NULL,
  ...
)

# S3 method for class 'fp_threshold_sweep'
plot(x, parameters = NULL, ...)
```

## Arguments

- occasions, detections, by, species:

  Passed to
  [`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).

- thresholds:

  Numeric confidence thresholds in `[0, 1]`, or names `"MaxF1"` and/or
  `"MaxPrecision"` to derive thresholds from `validation`. Defaults to
  both named methods. Each distinct threshold is fitted once.

- sites, verified, site, recorder:

  Passed to
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md);
  `site` must also be a column of `occasions`.

- occ, det, fp, random_recorder:

  Passed to
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

- prob:

  Width of the credible intervals in the summary.

- validation:

  Optional reviewed clips, as returned by
  [`fp_read_validation()`](https://provetelab.org/crywolf/reference/fp_read_validation.md),
  with numeric `confidence` in `[0, 1]` and logical (or binary 0/1)
  `verdict`. Required for named thresholds. Rows with missing scores or
  verdicts are dropped with a warning. When `species` is supplied,
  validation is filtered by `species_code` or `common_name`.

- ...:

  Further arguments for
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md),
  such as `chains`, `threads_per_chain`, or `seed`.

- x:

  An `fp_threshold_sweep` object.

- parameters:

  Parameters to plot; by default the occupancy slopes and `p10_bar`.

## Value

A data frame of class `fp_threshold_sweep` with one row per threshold
and parameter: the threshold, the share of occasions that are
detections, the parameter, its posterior mean and interval, and the
number of divergent transitions and largest R-hat of the fit.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
estimates against the threshold. For named thresholds, attributes
`threshold_selection` and `validation_metrics` contain the selected
methods and the metrics at every candidate threshold, respectively
(columns `threshold`, `tp`, `fp`, `fn`, `precision`, `recall`, `f1`).

## Details

Named methods evaluate every distinct validation score, counting scores
at or above the threshold as positive. Precision is TP / (TP + FP),
recall is TP / (TP + FN), and F1 is 2 TP / (2 TP + FP + FN). The
candidates are the observed scores, so each predicts at least one
positive. At least one reviewed positive is required. The selected
threshold is a score of the validation table, applied with `>=` to
`detections`; if the validation file stores rounded scores, the clips at
the cutoff can fall on the other side. Ties are resolved by greater
recall, then the smaller threshold. Recall and F1 refer to positives in
the supplied validation set: clips missed by the classifier or excluded
by an earlier score cutoff are not represented. Use a representative
validation sample for the target species. Validation selects thresholds
only; it is not passed as confirmation data to
[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md).

## References

Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo, A. V., Sugai,
L. S. M., & Provete, D. B. (2026). Methodological choices influence
ecological inference in passive acoustic monitoring of a Neotropical
nightjar. *Journal of Avian Biology*, in press.

## See also

[`fp_aggregate()`](https://provetelab.org/crywolf/reference/fp_aggregate.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# occasions: one row per one-minute recording, with site_night and covariates
sweep <- fp_threshold_sweep(occasions, detections, thresholds = c(0.1, 0.4, 0.7),
                            sites = nights, verified = verified,
                            site = "site_night", occ = ~ temp + rh, fp = ~ 1,
                            random_recorder = FALSE)
plot(sweep)
# Derive both thresholds from reviewed classifier clips:
reviews <- fp_read_validation("validation_by_clip.csv")
selected <- fp_threshold_sweep(occasions, detections, sites = nights,
                               site = "site_night", validation = reviews,
                               random_recorder = FALSE)
attr(selected, "threshold_selection")
} # }
```
