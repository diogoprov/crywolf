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
  thresholds,
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

  Confidence thresholds to try, for instance those that maximise F1 or
  precision on a validation set (this function does not derive them).

- sites, verified, site, recorder:

  Passed to
  [`fp_data()`](https://provetelab.org/crywolf/reference/fp_data.md);
  `site` must also be a column of `occasions`.

- occ, det, fp, random_recorder:

  Passed to
  [`fp_occu()`](https://provetelab.org/crywolf/reference/fp_occu.md).

- prob:

  Width of the credible intervals in the summary.

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
estimates against the threshold.

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
} # }
```
