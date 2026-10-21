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
#' @param thresholds Confidence thresholds to try, for instance those that
#'   maximise F1 or precision on a validation set (this function does not
#'   derive them).
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
#'   draws the estimates against the threshold.
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
#' }
fp_threshold_sweep <- function(occasions, detections, thresholds, sites,
                               verified = NULL, site = "site", by = "file",
                               species = NULL, recorder = NULL,
                               occ = ~ 1, det = ~ 1, fp = ~ 1,
                               random_recorder = !is.null(recorder),
                               prob = 0.95, ...) {
  a <- (1 - prob) / 2
  purrr::map(sort(thresholds), function(th) {
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
