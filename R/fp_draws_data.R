#' Detection histories over the draws of a classifier calibration
#'
#' Builds one version of the data per posterior draw of the calibration
#' exported by sheepdog (`calibration_draws.parquet`, the probability that
#' a prediction with a given score is a true detection, at each site and
#' species), for [fp_occu_mc()]. In each version an occasion is a detection
#' when at least one of its predictions is judged true under that draw.
#'
#' @param x An export read by [fp_read_sheepdog()], with calibration draws.
#'   Read it with `species` set to the target species.
#' @param rule How a prediction is judged true under a draw with
#'   probability `p` of being true: `"sample"` (default) draws it true with
#'   probability `p`, so that no threshold is involved; `"threshold"` takes
#'   it as true when `p >= p_min`, as `sheepdog::history_draws()` does.
#' @param p_min Threshold of `rule = "threshold"`.
#' @param n_draws Optional number of calibration draws to use (the first
#'   ones, in the order of their identifiers); all by default.
#' @param seed Optional seed of `rule = "sample"`. The random numbers of
#'   the session continue afterwards as if the call had drawn none.
#' @param occasions The occasions, one row per `site` and `occasion_id` of
#'   the export, with the column `site` names (the sites of the model, such
#'   as site-nights) and any detection covariates. Defaults to
#'   `x$occasions`. Every version carries a copy of it, so keep only the
#'   columns the model needs when there are many occasions and draws.
#' @param sites The sites table of [fp_data()], with the column `site`
#'   names. Defaults to `x$sites`.
#' @param site Column of `occasions` and `sites` that identifies the sites
#'   of the model.
#'
#' @details
#' **The rules.** `"sample"` judges each prediction on its own, so an
#' occasion with many low-scoring predictions is likely to hold at least one
#' judged true; it assumes that the predictions of an occasion are true or
#' false independently, given the calibration. `"threshold"` keeps only the
#' predictions the calibration is nearly sure of, so it gives fewer
#' detections: on the `curiango` data (custom classifier, minute
#' occasions), about 8% of the minutes with `"sample"` against about 1%
#' with `"threshold"`.
#'
#' When the calibration depends on the month (`"month"` among
#' `schema$calibration$keys`), a prediction takes the value of the month of
#' its recording, from the recording date of its occasion in the export.
#'
#' A prediction takes the calibration value of the grid point nearest its
#' score, as in sheepdog: grid index
#' `round((score - grid_min) / grid_step) + 1`, within the grid. A prediction without calibration draws (no draws
#' for its site and species, or month) cannot be scored: its occasion is
#' set to `NA` (not surveyed) in every version, with a warning, rather than
#' counted as a non-detection.
#'
#' The recordings annotated by ear (`x$verified`) are the verified
#' occasions of every version: they do not depend on the calibration. The
#' reviewed clips (`x$validation`) are not used here; they are the data of
#' the calibration.
#'
#' Histories built this way hold few false positives, so a false-positive
#' process fitted to them (the default of [fp_occu()]) is informed mostly by
#' its priors; compare with `fp = NULL`.
#'
#' @return A list of [fp_data()] objects, one per calibration draw, with
#'   the attributes `draws` (the draw identifiers) and `rule`.
#' @seealso [fp_occu_mc()], [fp_read_sheepdog()]
#' @section Experimental:
#' This function is experimental (new in 2.0.0, part of the integration with
#' 'sheepdog'): its arguments and results may change in a minor version.
#' See [crywolf-api].
#'
#' @export
#' @examples
#' \dontrun{
#' x <- fp_read_sheepdog("occupancy_export", species = "nyc_alb")
#' versions <- fp_draws_data(x, n_draws = 20, seed = 1)
#' mc <- fp_occu_mc(versions, occ = ~ 1, fp = NULL, threads_per_chain = 1)
#' summary(mc)
#' }
fp_draws_data <- function(x, rule = c("sample", "threshold"), p_min = 0.95,
                          n_draws = NULL, seed = NULL, occasions = x$occasions,
                          sites = x$sites, site = "site") {
  rule <- match.arg(rule)
  if (!inherits(x, "fp_sheepdog")) {
    stop("`x` must be an export read by fp_read_sheepdog().", call. = FALSE)
  }
  if (!"calibration_draws" %in% names(x$files) || is.null(x$schema$calibration)) {
    stop("The export has no calibration draws; write it with ",
         "sheepdog::export_occupancy_data(fit = ).", call. = FALSE)
  }
  .require_pkg("arrow")
  .check_df(occasions, "occasions", c("site", "occasion_id", site))
  cal <- as.data.frame(arrow::read_parquet(x$files[["calibration_draws"]]))
  # the origin of the grid with the exact bits sheepdog used: the value
  # read back from schema.json can differ in the last digit and send a
  # score half-way between two points to the other one
  grid_min <- min(cal$confidence[cal$grid_index == 1])
  det <- x$detections
  if (length(unique(det$species)) > 1) {
    stop("The export holds several species; read it with ",
         "fp_read_sheepdog(species = ) set to the target species.", call. = FALSE)
  }
  cal <- cal[cal$species %in% unique(det$species), , drop = FALSE]
  ids <- sort(unique(cal$draw))
  if (!is.null(n_draws)) ids <- utils::head(ids, n_draws)

  # a calibration by month (keys holding "month") gives P(TP) for each
  # month of recording, "YYYYMM" from the recording date of the occasion
  by_month <- "month" %in% x$schema$calibration$keys
  if (by_month != ("month" %in% names(cal))) {
    stop("schema.json and calibration_draws.parquet disagree on whether the ",
         "calibration depends on the month.", call. = FALSE)
  }
  if (by_month) {
    oe <- x$occasions
    mo <- substr(oe$recording_date, 1, 6)
    det_month <- mo[match(paste(det$site, det$occasion_id, sep = "\r"),
                          paste(oe$site, oe$occasion_id, sep = "\r"))]
  } else {
    det_month <- rep("", nrow(det))
    cal$month <- rep("", nrow(cal))
  }
  # P(TP) as a matrix: one row per site x species (x month) x grid point,
  # one column per draw
  key <- function(s, sp, m, g) paste(s, sp, m, g, sep = "\r")
  cal <- cal[cal$draw %in% ids, , drop = FALSE]
  rows <- unique(key(cal$site, cal$species, cal$month, cal$grid_index))
  if (anyDuplicated(cbind(key(cal$site, cal$species, cal$month, cal$grid_index),
                          cal$draw))) {
    stop("calibration_draws.parquet has more than one value for some key and ",
         "draw; is a key column missing from schema.json?", call. = FALSE)
  }
  P <- matrix(NA_real_, length(rows), length(ids))
  P[cbind(match(key(cal$site, cal$species, cal$month, cal$grid_index), rows),
          match(cal$draw, ids))] <- cal$p_tp
  step <- x$schema$calibration$grid_step
  g <- round((det$confidence - grid_min) / step) + 1
  g <- pmin(pmax(g, 1), max(cal$grid_index))
  occ_key <- paste(occasions$site, occasions$occasion_id, sep = "\r")
  if (anyDuplicated(occ_key)) {
    stop("`occasions` must have one row per site and occasion_id.", call. = FALSE)
  }
  o <- match(paste(det$site, det$occasion_id, sep = "\r"), occ_key)
  k <- match(key(det$site, det$species, det_month, g), rows)
  if (any(is.na(o))) {
    warning(sum(is.na(o)), " prediction(s) are not in `occasions` and were ",
            "left out.", call. = FALSE)
  }
  # an occasion with a prediction that has no calibration cannot be scored:
  # it is treated as not surveyed (y = NA), not as a non-detection
  unscored <- unique(o[is.na(k) & !is.na(o)])
  if (length(unscored)) {
    warning(sum(is.na(k) & !is.na(o)), " prediction(s) of site x species",
            if (by_month) " x month", " without calibration draws: their ",
            length(unscored), " occasion(s) are treated as not surveyed.",
            call. = FALSE)
  }
  keep <- !is.na(k) & !is.na(o) & !(o %in% unscored)
  if (!any(keep)) {
    stop("No prediction could be matched to the calibration draws.", call. = FALSE)
  }
  k <- k[keep]
  o <- o[keep]

  verified <- x$verified
  if (!is.null(verified) && site != "site") {
    v <- match(paste(verified$site, verified$occasion_id, sep = "\r"), occ_key)
    verified[[site]] <- occasions[[site]][v]
    if (anyNA(v)) {
      warning(sum(is.na(v)), " annotated recording(s) are not in `occasions` and ",
              "were left out.", call. = FALSE)
    }
    verified <- verified[!is.na(v), , drop = FALSE]
  }
  if (!is.null(verified)) verified <- verified[c(site, setdiff(names(verified), site))]
  if (!is.null(seed) && rule == "sample") {
    # the user's random numbers continue as if this call had drawn none
    old <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    }
    .set_seed(seed)
    on.exit(if (is.null(old)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old, envir = globalenv())
    }, add = TRUE)
  }

  out <- purrr::imap(ids, function(id, j) {
    p <- P[k, j]
    true <- if (rule == "sample") stats::runif(length(p)) < p else p >= p_min
    auto <- occasions
    auto$y <- as.integer(tabulate(o[true], nbins = nrow(occasions)) > 0)
    auto$y[unscored] <- NA
    # the design checks of fp_data() speak once, for the first version
    if (j == 1) return(fp_data(auto, sites, verified, site = site))
    withCallingHandlers(fp_data(auto, sites, verified, site = site),
                        crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  })
  structure(out, draws = ids, rule = rule)
}
