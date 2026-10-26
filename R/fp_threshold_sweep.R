#' Fit the model at several classifier thresholds
#'
#' The threshold at which a confidence score counts as a detection trades
#' false positives against missed detections, and it can change the
#' estimated environmental effects (Piatti et al. 2026). This function
#' builds the automated occasions at each threshold with [fp_aggregate()],
#' fits the same model to each, and collects the estimates, so that their
#' sensitivity to the threshold can be seen at once.
#'
#' @param occasions,detections,by,species Passed to [fp_aggregate()].
#' @param thresholds Numeric confidence thresholds in `[0, 1]`, or names
#'   `"MaxF1"` and/or `"MaxPrecision"` to derive thresholds from `validation`.
#'   Defaults to both named methods. Each distinct threshold is fitted once.
#' @param validation Optional reviewed clips, as returned by
#'   [fp_read_validation()], with numeric `confidence` in `[0, 1]` and logical
#'   (or binary 0/1) `verdict`. Required for named thresholds. Rows with missing
#'   scores or verdicts are dropped with a warning. When `species` is supplied,
#'   validation is filtered by `species_code` or `common_name`.
#'
#' @details Named methods evaluate every distinct validation score, counting
#'   scores at or above the threshold as positive. Precision is TP / (TP + FP),
#'   recall is TP / (TP + FN), and F1 is 2 TP / (2 TP + FP + FN). The
#'   candidates are the observed scores, so each predicts at least one
#'   positive. At least one reviewed positive is required. The selected
#'   threshold is a score of the validation table, applied with `>=` to
#'   `detections`; if the validation file stores rounded scores, the clips
#'   at the cutoff can fall on the other side. Ties are resolved by greater recall, then the smaller threshold.
#'   Recall and F1 refer to positives in the supplied validation set: clips
#'   missed by the classifier or excluded by an earlier score cutoff are not
#'   represented. Use a representative validation sample for the target species.
#'   Validation selects thresholds only; it is not passed as confirmation data
#'   to [fp_aggregate()].
#' @param sites,verified,site,recorder Passed to [fp_data()]; `site` must
#'   also be a column of `occasions`.
#' @param occ,det,fp,random_recorder Passed to [fp_occu()].
#' @param prob Width of the credible intervals in the summary.
#' @param ... Further arguments for [fp_occu()], such as `chains`,
#'   `threads_per_chain`, or `seed`.
#'
#' @return A data frame of class `fp_threshold_sweep` with one row per
#'   threshold and parameter: the threshold, the share of occasions that are
#'   detections, the parameter, its posterior mean and interval, and the
#'   number of divergent transitions and largest R-hat of the fit. [plot()]
#'   draws the estimates against the threshold. For named thresholds, attributes
#'   `threshold_selection` and `validation_metrics` contain the selected methods
#'   and the metrics at every candidate threshold, respectively (columns
#'   `threshold`, `tp`, `fp`, `fn`, `precision`, `recall`, `f1`).
#' @references Piatti, L., de Sousa, D. L. H., Oliveira, B. S., de Melo,
#'   A. V., Sugai, L. S. M., & Provete, D. B. (2026). Methodological choices
#'   influence ecological inference in passive acoustic monitoring of a
#'   Neotropical nightjar. *Journal of Avian Biology*, in press.
#' @seealso [fp_aggregate()]
#' @export
#' @examples
#' \dontrun{
#' # occasions: one row per one-minute recording, with site_night and covariates
#' sweep <- fp_threshold_sweep(occasions, detections, thresholds = c(0.1, 0.4, 0.7),
#'                             sites = nights, verified = verified,
#'                             site = "site_night", occ = ~ temp + rh, fp = ~ 1,
#'                             random_recorder = FALSE)
#' plot(sweep)
#' # Derive both thresholds from reviewed classifier clips:
#' reviews <- fp_read_validation("validation_by_clip.csv")
#' selected <- fp_threshold_sweep(occasions, detections, sites = nights,
#'                                site = "site_night", validation = reviews,
#'                                random_recorder = FALSE)
#' attr(selected, "threshold_selection")
#' }
fp_threshold_sweep <- function(occasions, detections,
                               thresholds = c("MaxF1", "MaxPrecision"), sites,
                               verified = NULL, site = "site", by = "file",
                               species = NULL, recorder = NULL,
                               occ = ~ 1, det = ~ 1, fp = ~ 1,
                               random_recorder = !is.null(recorder),
                               prob = 0.95, validation = NULL, ...) {
  selection <- metrics <- NULL
  if (is.character(thresholds)) {
    if (!length(thresholds)) {
      stop("`thresholds` must be nonempty numeric values in [0, 1] or named methods.",
           call. = FALSE)
    }
    methods <- match.arg(thresholds, c("MaxF1", "MaxPrecision"), several.ok = TRUE)
    if (is.null(validation)) {
      stop("Named thresholds (", paste(unique(methods), collapse = ", "), ") are ",
           "derived from reviewed clips: give them in `validation`, or give ",
           "numeric `thresholds`.", call. = FALSE)
    }
    metrics <- .threshold_metrics(validation, species)
    selection <- purrr::map(unique(methods), function(method) {
      metric <- if (method == "MaxF1") metrics$f1 else metrics$precision
      best <- order(-metric, -metrics$recall, metrics$threshold)[1]
      data.frame(method = method, metrics[best, ], row.names = NULL)
    }) |> purrr::list_rbind()
    thresholds <- selection$threshold
  } else {
    if (!is.numeric(thresholds) || !length(thresholds) || anyNA(thresholds) ||
        any(!is.finite(thresholds)) || any(thresholds < 0 | thresholds > 1)) {
      stop("`thresholds` must be nonempty numeric values in [0, 1] or named methods.",
           call. = FALSE)
    }
    if (!is.null(validation)) {
      warning("`validation` is used only to derive named thresholds; it is ",
              "ignored with numeric `thresholds`.", call. = FALSE)
    }
  }
  out <- purrr::map(sort(unique(thresholds)), function(th) {
    auto <- fp_aggregate(occasions, detections, threshold = th, by = by,
                         species = species)
    d <- fp_data(auto, sites, verified, site = site, recorder = recorder)
    m <- suppressWarnings(fp_occu(d, occ = occ, det = det, fp = fp,
                                  random_recorder = random_recorder, ...))
    s <- summary(m, prob = prob)$table
    keep <- grepl("^(occ_|det_|fp_|verified_)", s$variable) | s$variable == "p10_bar"
    s <- s[keep & s$variable != "fp_gap", ]
    data.frame(threshold = th, detection_rate = mean(auto$y),
               parameter = s$variable, mean = s$mean, lower = s$lower,
               upper = s$upper, divergences = m$diagnostics$divergences,
               max_rhat = m$diagnostics$max_rhat)
  }) |>
    purrr::list_rbind() |>
    structure(class = c("fp_threshold_sweep", "data.frame"), prob = prob)
  if (!is.null(selection)) {
    attr(out, "threshold_selection") <- selection
    attr(out, "validation_metrics") <- metrics
  }
  out
}

# Metrics use complete score groups so tied scores always cross together.
.threshold_metrics <- function(validation, species = NULL) {
  .check_df(validation, "validation", c("confidence", "verdict"))
  sp_cols <- intersect(c("species_code", "common_name"), names(validation))
  if (!is.null(species)) {
    if (!length(sp_cols)) {
      stop("`validation` needs species_code or common_name when `species` is supplied.",
           call. = FALSE)
    }
    hit <- purrr::map(sp_cols, function(cl) validation[[cl]] %in% species) |>
      purrr::reduce(`|`)
    validation <- validation[hit, , drop = FALSE]
  } else if (length(sp_cols) && length(unique(validation[[sp_cols[1]]])) > 1) {
    stop("`validation` holds several species; name the target with `species`.",
         call. = FALSE)
  }
  score <- validation$confidence
  truth <- validation$verdict
  if (!is.numeric(score) || any(!is.na(score) & (!is.finite(score) | score < 0 | score > 1))) {
    stop("`validation$confidence` must be numeric scores in [0, 1].", call. = FALSE)
  }
  if (!(is.logical(truth) || is.numeric(truth)) || any(!is.na(truth) & !truth %in% c(0, 1))) {
    stop("`validation$verdict` must be logical or binary 0/1.", call. = FALSE)
  }
  keep <- !is.na(score) & !is.na(truth)
  if (any(!keep)) warning(sum(!keep), " validation row(s) with missing confidence or verdict were dropped.",
                          call. = FALSE)
  score <- score[keep]
  truth <- truth[keep]
  if (!length(score) || !any(truth == 1)) {
    stop("`validation` must contain at least one reviewed positive with a score.",
         call. = FALSE)
  }
  ord <- order(score, decreasing = TRUE)
  score <- score[ord]
  truth <- truth[ord]
  ends <- which(!duplicated(score, fromLast = TRUE))
  tp <- cumsum(truth)[ends]
  fp <- ends - tp
  fn <- sum(truth) - tp
  out <- data.frame(threshold = score[ends], tp = tp, fp = fp, fn = fn,
                    precision = tp / ends, recall = tp / sum(truth),
                    f1 = 2 * tp / (2 * tp + fp + fn))
  out[order(out$threshold), , drop = FALSE]
}

#' @rdname fp_threshold_sweep
#' @param x An `fp_threshold_sweep` object.
#' @param parameters Parameters to plot; by default the occupancy slopes and
#'   `p10_bar`.
#' @export
plot.fp_threshold_sweep <- function(x, parameters = NULL, ...) {
  .require_pkg("ggplot2")
  if (is.null(parameters)) {
    parameters <- c(grep("^occ_", unique(x$parameter), value = TRUE), "p10_bar")
    parameters <- setdiff(parameters, "occ_(Intercept)")
  }
  x <- x[x$parameter %in% parameters, ]
  ggplot2::ggplot(x, ggplot2::aes(threshold, mean)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_pointrange(ggplot2::aes(ymin = lower, ymax = upper), colour = "#241042") +
    ggplot2::geom_line(colour = "#241042", alpha = 0.5) +
    ggplot2::facet_wrap(~ parameter, scales = "free_y") +
    ggplot2::labs(x = "Confidence threshold",
                  y = sprintf("Posterior mean and %g%% interval", 100 * attr(x, "prob"))) +
    ggplot2::theme_minimal()
}

utils::globalVariables("threshold")
