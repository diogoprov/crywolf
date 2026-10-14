#' Assemble data for a false-positive occupancy model
#'
#' Collects the automated detections, the verified detections, and the site
#' table into one object that [fp_occu()] can fit. All tables are in long
#' format, one row per occasion, which accommodates unequal numbers of
#' occasions per site without padding. Use [fp_data_wide()] if your data are
#' detection matrices, as in 'unmarked'.
#'
#' @param auto Data frame of the automated method (the one that can return
#'   false positives), one row per site and occasion. It needs a site column
#'   and a 0/1 detection column, and may hold observation covariates for the
#'   true-detection process and a recorder (or observer) identifier.
#' @param sites Data frame with one row per site: the site identifier plus any
#'   site covariates used in the occupancy (`occ`) or false-positive (`fp`)
#'   formulas of [fp_occu()]. Sites with no automated occasions are allowed and
#'   contribute only through their verified data, if any.
#' @param verified Optional data frame of the verified method (detections
#'   checked by an expert and therefore free of false positives), one row per
#'   site and occasion, with the same site and detection columns. Without it
#'   the false-positive rate is identified only by the priors, and
#'   [fp_occu()] warns about it.
#' @param site,y Names of the site and detection columns, shared by all
#'   tables.
#' @param recorder Optional name of a column of `auto` that identifies the
#'   recorder or observer, used for a random intercept on true detection.
#'
#' @return An object of class `fp_data`: a list with the three tables (rows
#'   with a missing detection dropped), the column names, and the vector
#'   `known_present`, which flags sites with at least one verified detection.
#' @seealso [fp_data_wide()], [fp_occu()], [simulate_fp()]
#' @export
#' @examples
#' sim <- simulate_fp(nsite = 50, nocc = 5, seed = 1)
#' sim$data
#'
#' # the same structure built by hand
#' d <- fp_data(auto = sim$data$auto, sites = sim$data$sites,
#'              verified = sim$data$verified, recorder = "recorder")
#' d
fp_data <- function(auto, sites, verified = NULL, site = "site", y = "y",
                    recorder = NULL) {
  .check_df(auto, "auto", c(site, y, recorder))
  .check_df(sites, "sites", site)
  if (anyDuplicated(sites[[site]])) {
    stop("`sites` must have one row per site; duplicated identifiers found in `",
         site, "`.", call. = FALSE)
  }
  auto <- auto[!is.na(auto[[y]]), , drop = FALSE]
  .check_binary(auto[[y]], "auto")
  .check_sites_known(auto[[site]], sites[[site]], "auto", site)

  if (!is.null(verified)) {
    .check_df(verified, "verified", c(site, y))
    verified <- verified[!is.na(verified[[y]]), , drop = FALSE]
    .check_binary(verified[[y]], "verified")
    .check_sites_known(verified[[site]], sites[[site]], "verified", site)
  }
  if (!is.null(recorder) && anyNA(auto[[recorder]])) {
    stop("Column `", recorder, "` of `auto` has missing values.", call. = FALSE)
  }

  known_present <- rep(0L, nrow(sites))
  if (!is.null(verified) && nrow(verified) > 0) {
    hit <- unique(verified[[site]][verified[[y]] == 1])
    known_present[match(hit, sites[[site]])] <- 1L
  }

  structure(
    list(auto = auto, verified = verified, sites = sites,
         site = site, y = y, recorder = recorder,
         known_present = known_present),
    class = "fp_data"
  )
}

#' Assemble data from detection matrices
#'
#' Convenience wrapper around [fp_data()] for data stored as site-by-occasion
#' matrices, the layout used by 'unmarked'. Missing values (`NA`) mark
#' occasions that were not sampled, so sites may differ in effort.
#'
#' @param y_auto Matrix of automated detections (sites in rows, occasions in
#'   columns), coded 0/1/`NA`.
#' @param y_verified Optional matrix of verified detections with the same
#'   number of rows; unverified sites are rows of `NA`.
#' @param site_covs Optional data frame of site covariates, one row per site,
#'   in the same order as the rows of `y_auto`.
#' @param obs_covs Optional named list of matrices, each the same size as
#'   `y_auto`, holding observation covariates of the automated method.
#' @param recorder Optional vector, one value per site, identifying the
#'   recorder or observer.
#'
#' @return An object of class `fp_data`.
#' @seealso [fp_data()]
#' @export
#' @examples
#' set.seed(1)
#' ya <- matrix(rbinom(40, 1, 0.3), 10, 4)
#' yv <- matrix(NA, 10, 2); yv[1:3, ] <- rbinom(6, 1, 0.5)
#' fp_data_wide(ya, yv, site_covs = data.frame(x = rnorm(10)))
fp_data_wide <- function(y_auto, y_verified = NULL, site_covs = NULL,
                         obs_covs = NULL, recorder = NULL) {
  y_auto <- as.matrix(y_auto)
  nsite <- nrow(y_auto)
  sites <- data.frame(site = seq_len(nsite))
  if (!is.null(site_covs)) {
    if (nrow(site_covs) != nsite) {
      stop("`site_covs` must have one row per row of `y_auto`.", call. = FALSE)
    }
    sites <- cbind(sites, site_covs)
  }
  to_long <- function(m) {
    data.frame(site = rep(seq_len(nrow(m)), times = ncol(m)),
               occasion = rep(seq_len(ncol(m)), each = nrow(m)),
               y = as.vector(m))
  }
  auto <- to_long(y_auto)
  if (!is.null(obs_covs)) {
    if (is.null(names(obs_covs)) || any(!nzchar(names(obs_covs)))) {
      stop("`obs_covs` must be a named list of matrices.", call. = FALSE)
    }
    purrr::iwalk(obs_covs, function(m, nm) {
      if (!identical(dim(as.matrix(m)), dim(y_auto))) {
        stop("Observation covariate `", nm, "` must have the same dimensions ",
             "as `y_auto`.", call. = FALSE)
      }
    })
    auto <- cbind(auto, purrr::map(obs_covs, function(m) as.vector(as.matrix(m))) |>
                        as.data.frame())
  }
  rec_col <- NULL
  if (!is.null(recorder)) {
    if (length(recorder) != nsite) {
      stop("`recorder` must have one value per site.", call. = FALSE)
    }
    auto$recorder <- recorder[auto$site]
    rec_col <- "recorder"
  }
  verified <- if (is.null(y_verified)) NULL else to_long(as.matrix(y_verified))
  fp_data(auto = auto, sites = sites, verified = verified, recorder = rec_col)
}

#' @export
print.fp_data <- function(x, ...) {
  ns <- nrow(x$sites)
  nver_sites <- if (is.null(x$verified)) 0L else length(unique(x$verified[[x$site]]))
  naive_auto <- mean(tapply(x$auto[[x$y]], x$auto[[x$site]], max))
  cat("<fp_data>\n")
  cat(sprintf("  sites:                  %d\n", ns))
  cat(sprintf("  automated occasions:    %d (naive occupancy %.2f)\n",
              nrow(x$auto), naive_auto))
  cat(sprintf("  verified sites:         %d (%d occasions, %d with a confirmed detection)\n",
              nver_sites, if (is.null(x$verified)) 0L else nrow(x$verified),
              sum(x$known_present)))
  if (!is.null(x$recorder)) {
    cat(sprintf("  recorders:              %d\n", length(unique(x$auto[[x$recorder]]))))
  }
  covs <- setdiff(names(x$sites), x$site)
  if (length(covs)) cat("  site covariates:       ", paste(covs, collapse = ", "), "\n")
  ocovs <- setdiff(names(x$auto), c(x$site, x$y, x$recorder, "occasion"))
  if (length(ocovs)) cat("  observation covariates:", paste(ocovs, collapse = ", "), "\n")
  invisible(x)
}

# validation helpers ----------------------------------------------------------

.check_df <- function(df, what, cols) {
  if (!is.data.frame(df)) stop("`", what, "` must be a data frame.", call. = FALSE)
  miss <- setdiff(cols, names(df))
  if (length(miss)) {
    stop("`", what, "` lacks column(s): ", paste(miss, collapse = ", "), ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.check_binary <- function(v, what) {
  if (!all(v %in% c(0, 1))) {
    stop("Detections in `", what, "` must be coded 0/1.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_sites_known <- function(v, known, what, site) {
  bad <- setdiff(unique(v), known)
  if (length(bad)) {
    stop(length(bad), " site(s) in `", what, "` are missing from `sites` (e.g. ",
         paste(utils::head(bad, 3), collapse = ", "), ").", call. = FALSE)
  }
  invisible(TRUE)
}
