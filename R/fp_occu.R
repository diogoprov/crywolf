#' Fit a Bayesian false-positive occupancy model
#'
#' Fits the multiple-methods false-positive occupancy model of Miller et al.
#' (2011) in Stan. Detections from an automated method may be false positives,
#' whereas detections from a verified subset are certain. The occupancy
#' probability depends on site covariates, true detection on observation
#' (or site) covariates and an optional recorder random intercept, and the
#' false-positive rate is either constant or a function of site covariates.
#'
#' @param data An `fp_data` object from [fp_data()], [fp_data_wide()], or
#'   [simulate_fp()].
#' @param occ One-sided formula for the logit of occupancy, using columns of
#'   `data$sites`.
#' @param det One-sided formula for the logit of the true-detection probability
#'   of the automated method. It may use columns of `data$auto` and of
#'   `data$sites`. It must keep the intercept.
#' @param fp One-sided formula for the logit of the false-positive rate, using
#'   columns of `data$sites`. `~ 1` (the default) fits a constant rate; `~ w`
#'   lets the rate change along `w`. Use `NULL` to fit a standard occupancy
#'   model without false positives. The intercept is fixed by the model (see
#'   Details) and cannot be removed.
#' @param random_recorder Logical. Add a recorder random intercept to true
#'   detection? Defaults to `TRUE` when `data` has a recorder column. Data
#'   from [simulate_fp()] always have one; with its default `sd_recorder = 0`
#'   set this to `FALSE`, because a random effect whose true variance is zero
#'   causes divergent transitions.
#' @param priors Prior settings from [fp_priors()].
#' @param chains,iter_warmup,iter_sampling,adapt_delta,seed Passed to the
#'   `$sample()` method of 'cmdstanr'.
#' @param refresh How often each chain reports its progress, in iterations.
#'   The default, `NULL`, reports every 10% of the iterations in interactive
#'   sessions and nothing otherwise; `0` silences it. Warm-up iterations are
#'   usually slower than sampling ones, so the first half of the run takes
#'   more than half of the time.
#' @param parallel_chains Number of chains run in parallel.
#' @param ... Further arguments passed to `$sample()`.
#'
#' @details
#' **Identifiability.** The intercept of the false-positive rate equals the
#' intercept of true detection minus a positive gap, \eqn{\gamma_0 = \alpha_1 -
#' \Delta}, so that at the reference values of the covariates (all zero) a
#' true detection is more likely than a false positive. This removes the label
#' switching between the occupied and unoccupied states. Centre the
#' covariates, so that the reference is an average site; [simulate_fp()]
#' already returns them standardized.
#'
#' **Verified data.** The verified detections anchor the false-positive rate.
#' With none, the model still runs, but the rate is informed only by the priors
#' and the gap, and the fit should not be trusted.
#'
#' **Observation-confirmation design.** When `data` has confirmed occasions
#' (argument `confirmed` of [fp_data()]), the model follows the
#' observation-confirmation design of Chambert et al. (2015): an expert checks
#' every observation of some automated occasions and finds only false ones,
#' only true ones, both, or none. False observations can then occur at
#' occupied sites too, and `det` and `fp` model \eqn{s_1} and \eqn{s_0}, the
#' probabilities of at least one true and at least one false observation on an
#' occasion; a site-level detection at an occupied site has probability
#' \eqn{p_{11} = 1 - (1 - s_1)(1 - s_0)}, and `p11_ref` reports it at the
#' reference. The confirmed occasions identify \eqn{s_0} directly, so verified
#' data are optional; they can be combined with this design. The gap still
#' applies, so that \eqn{s_1 > s_0} at the reference.
#'
#' **When to model the false-positive rate.** In the simulations behind this
#' package, a constant rate gave calibrated intervals for the occupancy slope
#' when false positives were homogeneous; random variation of the rate among
#' sites made the estimate of the rate unreliable but largely spared the slope.
#' When the rate tracked an environmental gradient correlated with the
#' occupancy covariate and the confusion was strong, the slope was biased and
#' verifying more sites did not fix it. Placing the occupancy covariate in `fp`
#' recovered part of the calibration, and placing the gradient itself
#' recovered all of it. Use [fp_calibrate()] to check your own design.
#'
#' **Parameter names.** Coefficients are labelled by process and column of the
#' design matrix: `occ_*` (occupancy), `det_*` (true detection), `fp_*` (false
#' positive, where `fp_(Intercept)` is \eqn{\gamma_0}), plus `fp_gap`
#' (\eqn{\Delta}), `sd_recorder`, `p_verified` (detection probability of the
#' verified method), `p11_ref` and `p10_ref` (true- and false-positive
#' probabilities at the reference), and `p10_bar` (false-positive probability
#' averaged over the sites).
#'
#' @return An object of class `fp_occu`, with methods for [base::print()],
#'   [base::summary()], [stats::coef()], [stats::predict()], [loo::loo()], and
#'   [posterior::as_draws_df()]. The 'cmdstanr' fit is in element `fit`.
#' @references
#' Chambert, T., Miller, D. A. W., & Nichols, J. D. (2015). Modeling false
#' positive detections in species occurrence data under different study
#' designs. *Ecology*, 96, 332–339. \doi{10.1890/14-1507.1}
#'
#' Miller, D. A., Nichols, J. D., McClintock, B. T., Grant, E. H. C.,
#' Bailey, L. L., & Weir, L. A. (2011). Improving occupancy estimation when
#' two types of observational error occur: non-detection and species
#' misidentification. *Ecology*, 92, 1422–1428. \doi{10.1890/10-1396.1}
#' @seealso [fp_data()], [fp_priors()], [fp_calibrate()]
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
#'                    conf_cor = 0.7, seed = 42)
#'
#' # constant false-positive rate (the simulated recorders have no effect,
#' # so the recorder term is left out)
#' m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' # false-positive rate that tracks the occupancy covariate
#' m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#'
#' coef(m0)["occ_x"]; coef(m1)["occ_x"]   # truth: -0.8
#' loo::loo_compare(loo(m0), loo(m1))
#' }
fp_occu <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
                    random_recorder = !is.null(data$recorder),
                    priors = fp_priors(),
                    chains = 4, parallel_chains = chains,
                    iter_warmup = 1000, iter_sampling = 1000,
                    adapt_delta = 0.95, seed = NULL, refresh = NULL, ...) {
  if (!inherits(data, "fp_data")) {
    stop("`data` must be an `fp_data` object; see `fp_data()`.", call. = FALSE)
  }
  no_ver <- is.null(data$verified) || nrow(data$verified) == 0
  no_conf <- is.null(data$confirmed) || all(is.na(data$auto[[data$confirmed[1]]]))
  if (no_ver && no_conf) {
    warning("No verified data: the false-positive rate is identified only by ",
            "the priors. Treat the estimates with caution.", call. = FALSE)
  }
  sd_list <- .build_stan_data(data, occ, det, fp, random_recorder, priors)
  if (is.null(refresh)) {
    refresh <- if (interactive()) max(1, round((iter_warmup + iter_sampling) / 10)) else 0
  }
  mod <- fp_compile()
  fit <- mod$sample(data = sd_list$stan_data, chains = chains,
                    parallel_chains = parallel_chains,
                    iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                    adapt_delta = adapt_delta, seed = seed, refresh = refresh,
                    ...)
  draws <- .pretty_draws(fit, sd_list)
  diag <- .diagnostics(fit, draws)

  out <- structure(
    list(fit = fit, draws = draws, stan_data = sd_list$stan_data,
         formulas = list(occ = occ, det = det, fp = fp),
         design = sd_list$design, data = data,
         random_recorder = isTRUE(sd_list$stan_data$estimate_re == 1),
         priors = priors, diagnostics = diag, call = match.call()),
    class = "fp_occu"
  )
  .warn_diagnostics(diag)
  out
}

#' Prior settings for fp_occu()
#'
#' Normal priors on the logit scale. The defaults are weakly informative for
#' standardized covariates.
#'
#' @param intercept_sd Standard deviation of the Normal(0, sd) priors on the
#'   intercepts of occupancy, true detection, and verified detection.
#' @param slope_sd Standard deviation of the Normal(0, sd) priors on all
#'   slopes, including those of the false-positive rate.
#' @param gap_mean,gap_sd Mean and standard deviation of the Normal prior,
#'   truncated at zero, on the gap between the true-detection and
#'   false-positive intercepts. The default, Normal(2, 1.5), puts the
#'   false-positive rate at the reference roughly an order of magnitude below
#'   true detection, while leaving room for the data.
#' @param recorder_sd Scale of the half-Normal(0, sd) prior on the standard
#'   deviation of the recorder random intercept.
#'
#' @return A list of class `fp_priors`.
#' @export
#' @examples
#' fp_priors()
#' fp_priors(slope_sd = 1)
fp_priors <- function(intercept_sd = 1.5, slope_sd = 2.5,
                      gap_mean = 2, gap_sd = 1.5, recorder_sd = 2) {
  stopifnot(intercept_sd > 0, slope_sd > 0, gap_sd > 0, recorder_sd > 0)
  structure(list(intercept_sd = intercept_sd, slope_sd = slope_sd,
                 gap_mean = gap_mean, gap_sd = gap_sd, recorder_sd = recorder_sd),
            class = "fp_priors")
}

#' @export
print.fp_priors <- function(x, ...) {
  cat("<fp_priors> (logit scale)\n")
  cat(sprintf("  intercepts      ~ Normal(0, %g)\n", x$intercept_sd))
  cat(sprintf("  slopes          ~ Normal(0, %g)\n", x$slope_sd))
  cat(sprintf("  fp_gap          ~ Normal(%g, %g), truncated at 0\n", x$gap_mean, x$gap_sd))
  cat(sprintf("  sd_recorder     ~ half-Normal(0, %g)\n", x$recorder_sd))
  invisible(x)
}

# internal: data list for Stan ------------------------------------------------

.design <- function(form, df, what) {
  if (!inherits(form, "formula") || length(form) != 2) {
    stop("`", what, "` must be a one-sided formula, e.g. ~ x.", call. = FALSE)
  }
  vars <- all.vars(form)
  miss <- setdiff(vars, names(df))
  if (length(miss)) {
    stop("Variable(s) in `", what, "` not found: ", paste(miss, collapse = ", "),
         ".", call. = FALSE)
  }
  mf <- stats::model.frame(form, df, na.action = stats::na.pass)
  if (anyNA(mf)) {
    stop("Missing values in the covariates of `", what, "`.", call. = FALSE)
  }
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  list(X = X, terms = tt, xlevels = stats::.getXlevels(tt, mf),
       contrasts = attr(X, "contrasts"),
       intercept = attr(tt, "intercept") == 1)
}

.build_stan_data <- function(data, occ, det, fp, random_recorder, priors) {
  s <- data$site
  sites <- data$sites
  auto <- data$auto
  nsite <- nrow(sites)

  d_occ <- .design(occ, sites, "occ")
  if (!d_occ$intercept) stop("`occ` must keep the intercept.", call. = FALSE)

  site2 <- match(auto[[s]], sites[[s]])
  # the likelihood sums each site's occasions over contiguous rows
  ord <- order(site2)
  auto <- auto[ord, , drop = FALSE]
  site2 <- site2[ord]
  # det may use site covariates: bring them into the automated table
  site_only <- setdiff(names(sites), names(auto))
  auto_full <- cbind(auto, sites[site2, site_only, drop = FALSE])
  d_det <- .design(det, auto_full, "det")
  if (!d_det$intercept) stop("`det` must keep the intercept.", call. = FALSE)

  estimate_fp <- !is.null(fp)
  if (estimate_fp) {
    d_fp <- .design(fp, sites, "fp")
    if (!d_fp$intercept) {
      stop("`fp` must keep the intercept; it is tied to the detection ",
           "intercept by the gap (see ?fp_occu).", call. = FALSE)
    }
    Xfp <- d_fp$X[, -1, drop = FALSE]
  } else {
    d_fp <- NULL
    Xfp <- matrix(0, nsite, 0)
  }

  use_re <- isTRUE(random_recorder)
  if (use_re && is.null(data$recorder)) {
    stop("`random_recorder = TRUE` but `data` has no recorder column.", call. = FALSE)
  }
  if (use_re) {
    rec_f <- factor(auto[[data$recorder]])
    rec2 <- as.integer(rec_f)
    rec_levels <- levels(rec_f)
  } else {
    rec2 <- rep(1L, nrow(auto))
    rec_levels <- character(0)
  }

  ver <- data$verified
  if (is.null(ver)) {
    y1 <- integer(0); site1 <- integer(0)
  } else {
    y1 <- as.integer(ver[[data$y]]); site1 <- match(ver[[s]], sites[[s]])
  }

  # observation confirmation: checked occasions leave the automated rows and
  # enter with their confirmation state
  conf <- data$confirmed
  is3 <- if (is.null(conf)) rep(FALSE, nrow(auto)) else !is.na(auto[[conf[1]]])
  obs_confirm <- any(is3)
  if (obs_confirm && !estimate_fp) {
    stop("The observation-confirmation design needs a false-positive ",
         "process; do not use `fp = NULL`.", call. = FALSE)
  }
  v3 <- if (obs_confirm) {
    as.integer(2L * (auto[[conf[1]]][is3] > 0) + (auto[[conf[2]]][is3] > 0))
  } else integer(0)
  site3 <- site2[is3]
  rec3 <- rec2[is3]
  Z3 <- d_det$X[is3, , drop = FALSE]
  auto <- auto[!is3, , drop = FALSE]
  site2 <- site2[!is3]
  rec2 <- rec2[!is3]
  X_det <- d_det$X[!is3, , drop = FALSE]

  # collapse occasions with identical covariates within a site into
  # binomial rows (exact, and much faster when detection is site-level)
  y_raw <- as.integer(auto[[data$y]])
  key <- do.call(paste, c(list(site2, rec2), as.data.frame(X_det),
                          list(sep = "\r")))
  grp <- match(key, unique(key))
  first <- !duplicated(grp)
  n2 <- tabulate(grp)
  y2 <- as.integer(rowsum(y_raw, grp, reorder = FALSE))
  site2 <- site2[first]
  rec2 <- rec2[first]
  Z <- X_det[first, , drop = FALSE]
  len2 <- tabulate(site2, nbins = nsite)
  start2 <- cumsum(c(1L, len2))[seq_len(nsite)]

  stan_data <- list(
    nsite = nsite, Kbeta = ncol(d_occ$X), X = unname(d_occ$X),
    N2 = length(n2), n2 = as.integer(n2), y2 = y2, site2 = as.integer(site2),
    start2 = as.integer(start2), len2 = as.integer(len2),
    Kalpha = ncol(Z), Z = unname(Z),
    nrec = max(1L, length(rec_levels)), rec2 = as.integer(rec2),
    estimate_re = as.integer(use_re),
    Kfp = ncol(Xfp), Xfp = unname(Xfp), estimate_fp = as.integer(estimate_fp),
    N1 = length(y1), y1 = y1, site1 = as.integer(site1),
    known_present = as.integer(data$known_present),
    obs_confirm = as.integer(obs_confirm), N3 = length(v3), v3 = v3,
    site3 = as.integer(site3), Z3 = unname(Z3), rec3 = as.integer(rec3),
    prior_sd_int = priors$intercept_sd, prior_sd_slope = priors$slope_sd,
    gap_mean = priors$gap_mean, gap_sd = priors$gap_sd,
    recorder_sd = priors$recorder_sd
  )
  design <- list(occ = d_occ, det = d_det, fp = d_fp,
                 occ_names = colnames(d_occ$X), det_names = colnames(d_det$X),
                 fp_names = if (estimate_fp) colnames(d_fp$X) else character(0),
                 recorders = rec_levels)
  list(stan_data = stan_data, design = design)
}

# internal: readable parameter names ------------------------------------------

.name_map <- function(design, estimate_fp, use_re, has_verified = TRUE) {
  m <- c(
    stats::setNames(paste0("occ_", design$occ_names),
                    sprintf("beta[%d]", seq_along(design$occ_names))),
    stats::setNames(paste0("det_", design$det_names),
                    sprintf("alpha[%d]", seq_along(design$det_names))),
    if (has_verified) c("delta0" = "verified_(Intercept)", "r11" = "p_verified"),
    c("p11_ref" = "p11_ref")
  )
  if (estimate_fp) {
    slopes <- design$fp_names[-1]
    m <- c(m, c("gamma0" = "fp_(Intercept)"))
    if (length(slopes)) {
      m <- c(m, stats::setNames(paste0("fp_", slopes),
                                sprintf("gamma_fp[%d]", seq_along(slopes))))
    }
    m <- c(m, c("fp_gap" = "fp_gap", "p10_ref" = "p10_ref", "p10_bar" = "p10_bar"))
  }
  if (use_re) m <- c(m, c("sigma_u" = "sd_recorder"))
  m
}

.pretty_draws <- function(fit, sd_list) {
  sd <- sd_list$stan_data
  map <- .name_map(sd_list$design, sd$estimate_fp == 1, sd$estimate_re == 1,
                   has_verified = sd$N1 > 0)
  dr <- posterior::as_draws_df(fit$draws(variables = unique(sub("\\[.*", "", names(map)))))
  dr <- posterior::subset_draws(dr, variable = names(map))
  posterior::variables(dr) <- unname(map[posterior::variables(dr)])
  dr
}

.diagnostics <- function(fit, draws) {
  ds <- tryCatch(fit$diagnostic_summary(quiet = TRUE), error = function(e) NULL)
  summ <- posterior::summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")
  list(
    divergences = if (is.null(ds)) NA_integer_ else sum(ds$num_divergent),
    max_treedepth = if (is.null(ds)) NA_integer_ else sum(ds$num_max_treedepth),
    max_rhat = max(summ$rhat, na.rm = TRUE),
    min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(summ$ess_tail, na.rm = TRUE)
  )
}

.warn_diagnostics <- function(diag) {
  msgs <- c(
    if (isTRUE(diag$divergences > 0))
      sprintf("%d divergent transition(s); consider a higher `adapt_delta`.",
              diag$divergences),
    if (isTRUE(diag$max_rhat > 1.01))
      sprintf("Largest R-hat is %.3f (> 1.01); run longer chains.", diag$max_rhat),
    if (isTRUE(diag$min_ess_bulk < 400))
      sprintf("Smallest bulk ESS is %.0f (< 400).", diag$min_ess_bulk)
  )
  if (length(msgs)) warning(paste(msgs, collapse = "\n"), call. = FALSE)
  invisible(NULL)
}
