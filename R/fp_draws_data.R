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
#' predictions the calibration is nearly sure of. On the `curiango` data
#' (custom classifier, one calibration draw), `"sample"` gave a detection
#' in 18% of the recordings and `"threshold"` in 0.4%, against 21% of the
#' recordings annotated by ear holding the species.
#'
#' A prediction takes the calibration value of the grid point nearest its
#' score, as in sheepdog: grid index `round((score - grid_min) / grid_step)
#' + 1`, within the grid. Predictions of a site and species without
#' calibration draws are left out, with a warning (no value is imputed).
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

  # P(TP) as a matrix: one row per site x species x grid point, one column
  # per draw
  key <- function(s, sp, g) paste(s, sp, g, sep = "\r")
  cal <- cal[cal$draw %in% ids, , drop = FALSE]
  rows <- unique(key(cal$site, cal$species, cal$grid_index))
  P <- matrix(NA_real_, length(rows), length(ids))
  P[cbind(match(key(cal$site, cal$species, cal$grid_index), rows),
          match(cal$draw, ids))] <- cal$p_tp
  step <- x$schema$calibration$grid_step
  g <- round((det$confidence - grid_min) / step) + 1
  g <- pmin(pmax(g, 1), max(cal$grid_index))
  k <- match(key(det$site, det$species, g), rows)
  if (anyNA(k)) {
    warning(sum(is.na(k)), " prediction(s) of site x species without calibration ",
            "draws were left out.", call. = FALSE)
    det <- det[!is.na(k), , drop = FALSE]
    k <- k[!is.na(k)]
  }
  occ_key <- paste(occasions$site, occasions$occasion_id, sep = "\r")
  if (anyDuplicated(occ_key)) {
    stop("`occasions` must have one row per site and occasion_id.", call. = FALSE)
  }
  o <- match(paste(det$site, det$occasion_id, sep = "\r"), occ_key)
  if (anyNA(o)) {
    warning(sum(is.na(o)), " prediction(s) are not in `occasions` and were ",
            "left out.", call. = FALSE)
    k <- k[!is.na(o)]
    o <- o[!is.na(o)]
  }

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
    # the design checks of fp_data() speak once, for the first version
    if (j == 1) return(fp_data(auto, sites, verified, site = site))
    withCallingHandlers(fp_data(auto, sites, verified, site = site),
                        crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  })
  structure(out, draws = ids, rule = rule)
}
