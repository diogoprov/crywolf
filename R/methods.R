#' Methods for fp_occu fits
#'
#' Summaries, point estimates, predictions, and model comparison for objects
#' returned by [fp_occu()].
#'
#' @param x,object An `fp_occu` object.
#' @param prob Width of the equal-tailed credible interval.
#' @param digits Number of digits to print.
#' @param ... Unused, or passed on (see each method).
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, seed = 1)
#' m <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' m                                   # posterior summary and diagnostics
#' summary(m, prob = 0.9)$table
#' coef(m, "occ")
#' predict(m, newdata = data.frame(x = c(-1, 0, 1)), process = "occ")
#' loo(m)
#' draws <- as_draws_df(m)
#' }
#' @name fp_occu-methods
NULL

#' @rdname fp_occu-methods
#' @export
print.fp_occu <- function(x, digits = 3, ...) {
  sd <- x$stan_data
  cat("Bayesian false-positive occupancy model (",
      if (isTRUE(sd$obs_confirm == 1)) "observation-confirmation" else "multiple-methods",
      " design)\n\n", sep = "")
  f <- x$formulas
  cat("  occupancy:      ", .fmt_formula(f$occ), "\n")
  cat("  true detection: ", .fmt_formula(f$det),
      if (x$random_recorder) "+ (1 | recorder)" else "", "\n")
  cat("  false positive: ",
      if (is.null(f$fp)) "none (standard occupancy model)" else .fmt_formula(f$fp),
      "\n")
  if (isTRUE(sd$N3 > 0)) {
    cat(sprintf("  data:            %d sites, %d automated (%d confirmed) and %d verified occasions\n\n",
                sd$nsite, sum(sd$n2) + sd$N3, sd$N3, sd$N1))
  } else {
    cat(sprintf("  data:            %d sites, %d automated and %d verified occasions\n\n",
                sd$nsite, sum(sd$n2), sd$N1))
  }
  s <- summary(x)
  print(as.data.frame(s$table), digits = digits, row.names = FALSE)
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu-methods
#' @return `summary()` returns an object of class `summary.fp_occu` whose
#'   element `table` is a tibble with the posterior mean, standard deviation,
#'   credible interval, R-hat, and effective sample sizes of every parameter.
#' @export
summary.fp_occu <- function(object, prob = 0.95, ...) {
  a <- (1 - prob) / 2
  tab <- posterior::summarise_draws(
    object$draws,
    mean = mean, sd = stats::sd,
    lower = function(v) unname(stats::quantile(v, a)),
    upper = function(v) unname(stats::quantile(v, 1 - a)),
    rhat = posterior::rhat, ess_bulk = posterior::ess_bulk,
    ess_tail = posterior::ess_tail
  )
  structure(list(table = tibble::as_tibble(tab), prob = prob,
                 diagnostics = object$diagnostics, call = object$call),
            class = "summary.fp_occu")
}

#' @rdname fp_occu-methods
#' @export
print.summary.fp_occu <- function(x, digits = 3, ...) {
  cat(sprintf("Posterior summary (%g%% equal-tailed intervals)\n\n", 100 * x$prob))
  print(as.data.frame(x$table), digits = digits, row.names = FALSE)
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu-methods
#' @param type Which coefficients to return: `"all"` (default), `"occ"`,
#'   `"det"`, or `"fp"`.
#' @return `coef()` returns a named vector of posterior means.
#' @export
coef.fp_occu <- function(object, type = c("all", "occ", "det", "fp"), ...) {
  type <- match.arg(type)
  m <- vapply(posterior::variables(object$draws),
              function(v) mean(object$draws[[v]]), numeric(1))
  if (type == "all") return(m)
  m[startsWith(names(m), paste0(type, "_")) & names(m) != "fp_gap"]
}

#' @rdname fp_occu-methods
#' @param newdata Data frame with the covariates of the requested process. If
#'   missing, the sites used in the fit.
#' @param process `"occ"` for the occupancy probability or `"fp"` for the
#'   false-positive rate, both at the site level.
#' @param summary Logical. Return a summary per row (`TRUE`) or the matrix of
#'   posterior draws (rows are draws, columns are rows of `newdata`)?
#' @return `predict()` returns a data frame with the posterior mean and
#'   credible interval of the probability for each row of `newdata`, or a
#'   matrix of draws.
#' @export
predict.fp_occu <- function(object, newdata, process = c("occ", "fp"),
                            prob = 0.95, summary = TRUE, ...) {
  process <- match.arg(process)
  if (process == "fp" && is.null(object$formulas$fp)) {
    stop("This model has no false-positive process.", call. = FALSE)
  }
  if (missing(newdata)) newdata <- object$data$sites
  des <- object$design[[process]]
  tt <- stats::delete.response(des$terms)
  mf <- stats::model.frame(tt, newdata, xlev = des$xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = des$contrasts)

  vars <- paste0(process, "_", colnames(X))
  B <- as.matrix(posterior::as_draws_matrix(
    posterior::subset_draws(object$draws, variable = vars)))[, vars, drop = FALSE]
  eta <- B %*% t(X)
  p <- stats::plogis(eta)
  if (!summary) return(p)
  a <- (1 - prob) / 2
  data.frame(
    newdata,
    estimate = colMeans(p),
    lower = apply(p, 2, stats::quantile, probs = a),
    upper = apply(p, 2, stats::quantile, probs = 1 - a),
    row.names = NULL
  )
}

#' @rdname fp_occu-methods
#' @return `loo()` returns a `psis_loo` object from 'loo'. The pointwise
#'   log-likelihood is computed per site, with the latent state marginalized,
#'   so the comparison is leave-one-site-out.
#' @export
loo.fp_occu <- function(x, ...) {
  x$fit$loo(variables = "log_lik", ...)
}

#' @rdname fp_occu-methods
#' @return `as_draws_df()` returns the posterior draws with the readable
#'   parameter names.
#' @export
as_draws_df.fp_occu <- function(x, ...) {
  x$draws
}

# helpers ---------------------------------------------------------------------

.fmt_formula <- function(f) paste(deparse(f, width.cutoff = 500L), collapse = "")

.print_diag <- function(d) {
  cat(sprintf(paste0("\n  Sampler: %s divergent transitions, max R-hat %.3f, ",
                     "min bulk/tail ESS %.0f/%.0f\n"),
              format(d$divergences), d$max_rhat, d$min_ess_bulk, d$min_ess_tail))
}
