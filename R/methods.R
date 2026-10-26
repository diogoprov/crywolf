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
      if (x$random_recorder) "+ (1 | recorder)" else "",
      if (isTRUE(x$random_det)) "+ (1 | site)" else "", "\n")
  cat("  false positive: ",
      if (is.null(f$fp)) "none (standard occupancy model)" else .fmt_formula(f$fp),
      if (isTRUE(x$random_fp)) "+ (1 | site)" else "", "\n")
  if (!is.null(f$ver) && length(all.vars(f$ver))) {
    cat("  verified:        ", .fmt_formula(f$ver), "\n")
  }
  if (!is.null(x$markov)) {
    cat("  runs:             previous detection (y_prev), occasions ordered by",
        x$markov, "\n")
  }
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
#'   missing, the sites used in the fit; for a false-positive rate with
#'   covariates measured on each occasion, the automated occasions, with the
#'   covariates of their sites.
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
  if (missing(newdata)) {
    d <- object$data
    newdata <- d$sites
    if (process == "fp" && isTRUE(object$stan_data$fp_by_row == 1)) {
      site_only <- setdiff(names(d$sites), names(d$auto))
      newdata <- cbind(d$auto, d$sites[match(d$auto[[d$site]], d$sites[[d$site]]),
                                       site_only, drop = FALSE])
    }
  }
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
#' @param reloo Logical. Refit the model without each site whose Pareto
#'   \eqn{k} exceeds `k_threshold`, and replace the approximate
#'   leave-one-out value of that site by the exact one. Each such site costs
#'   one fit, with the sampler settings of `x`. A site that is the only one
#'   with its recorder, or with a level of a factor covariate, cannot be
#'   predicted from a fit without it; it keeps its approximate value, with a
#'   warning.
#' @param k_threshold Pareto \eqn{k} above which a site is refitted.
#' @param integrate Logical. With a site random effect on the false-positive
#'   rate (`random_fp = TRUE`), on true detection (`random_det = TRUE`), or
#'   both, integrate each site's effects out of its log-likelihood (by
#'   Gauss-Hermite quadrature) before the leave-one-out computation. The
#'   default, `TRUE`, is the right choice for leaving a site out: a
#'   left-out site's effects are not known from its own data. With
#'   `FALSE`, the log-likelihood is conditional on the fitted site effects,
#'   which gives many sites with high Pareto \eqn{k}. Without a site
#'   random effect, `integrate` has no effect.
#' @param nodes Number of quadrature nodes for `integrate`, per site effect.
#'   The default, `NULL`, uses 100 with one site effect and 30 per effect
#'   with both (a grid of 900 nodes). With one effect and 100 nodes, the
#'   log-likelihood of a site was within 6e-4 of its value with 400 nodes in
#'   a test with 200 sites. The cost grows with the number of nodes in the
#'   grid, `nodes` with one effect and `nodes^2` with two.
#' @param seed Seed of the refits.
#' @return `loo()` returns a `psis_loo` object from 'loo'. The pointwise
#'   log-likelihood is computed per site, with the latent state marginalized,
#'   so the comparison is leave-one-site-out. Pareto \eqn{k} above 0.7 marks
#'   sites where the approximation is unreliable; `reloo = TRUE` computes
#'   their values exactly, and their \eqn{k} is then set to 0.
#' @export
loo.fp_occu <- function(x, reloo = FALSE, k_threshold = 0.7, seed = NULL,
                        integrate = TRUE, nodes = NULL, ...) {
  if ((isTRUE(x$random_fp) || isTRUE(x$random_det)) && integrate) {
    pars <- .ppc_pars(x$fit, x$stan_data, NULL)
    llm <- .integrate_site_effect(x$stan_data, pars, nodes)
    # iterations and chains as saved (thinning or a failed chain change them)
    dims <- dim(x$fit$draws("lp__"))[1:2]
    ll <- array(llm, dim = c(dims, ncol(llm)))
  } else {
    ll <- posterior::as_draws_array(x$fit$draws("log_lik"))
  }
  # each site scaled by its own maximum, so that no likelihood underflows
  site_max <- purrr::map_dbl(seq_len(dim(ll)[3]), function(j) max(ll[, , j]))
  r_eff <- loo::relative_eff(exp(sweep(ll, 3, site_max, "-")))
  l <- loo::loo(ll, r_eff = r_eff, ...)
  bad <- which(l$diagnostics$pareto_k > k_threshold)
  if (!reloo || !length(bad)) return(l)
  message("Refitting the model without each of ", length(bad), " site(s) ",
          "with Pareto k > ", k_threshold, ".")
  llm <- matrix(ll, ncol = dim(ll)[3])
  elpd <- purrr::map_dbl(bad, function(i) .elpd_without_site(x, i, seed, nodes))
  done <- bad[!is.na(elpd)]
  elpd <- elpd[!is.na(elpd)]
  if (!length(done)) return(l)
  lpd <- purrr::map_dbl(done, function(i) .log_mean_exp(llm[, i]))
  l$pointwise[done, "elpd_loo"] <- elpd
  l$pointwise[done, "p_loo"] <- lpd - elpd
  l$pointwise[done, "looic"] <- -2 * elpd
  # exact values: no Monte Carlo error from importance sampling
  l$pointwise[done, "mcse_elpd_loo"] <- 0
  l$pointwise[done, "influence_pareto_k"] <- 0
  l$diagnostics$pareto_k[done] <- 0
  n <- nrow(l$pointwise)
  cols <- c("elpd_loo", "p_loo", "looic")
  l$estimates[cols, "Estimate"] <- colSums(l$pointwise[, cols])
  l$estimates[cols, "SE"] <- sqrt(n * purrr::map_dbl(cols, function(cl) {
    stats::var(l$pointwise[, cl])
  }))
  l
}

#' Model weights by stacking
#'
#' Weights of several fitted models for combining their predictions, by
#' stacking of their leave-one-site-out predictive densities (Yao et al.
#' 2018), or by pseudo-BMA. A thin wrapper of [loo::loo_model_weights()]
#' that accepts fits or `loo()` results.
#'
#' @param ... Named `fp_occu` fits, or the results of [loo()] on them.
#' @param method `"stacking"` (default) or `"pseudobma"`.
#'
#' @return A vector of weights summing to one, named by the arguments in
#'   `...`.
#' @references Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using
#'   stacking to average Bayesian predictive distributions. *Bayesian
#'   Analysis*, 13, 917–1007. \doi{10.1214/17-BA1091}
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_fp(nsite = 200, nocc = 10, conf_effect = 1.5,
#'                    conf_cor = 0.7, seed = 1)
#' m0 <- fp_occu(sim$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE)
#' m1 <- fp_occu(sim$data, occ = ~ x, fp = ~ x, random_recorder = FALSE)
#' fp_loo_weights(constant = m0, fp_on_x = m1)
#' }
fp_loo_weights <- function(..., method = c("stacking", "pseudobma")) {
  method <- match.arg(method)
  fits <- list(...)
  if (is.null(names(fits)) || any(!nzchar(names(fits)))) {
    stop("Name every model, e.g. fp_loo_weights(constant = m0, fp_x = m1).",
         call. = FALSE)
  }
  loos <- purrr::map(fits, function(f) if (inherits(f, "fp_occu")) loo(f) else f)
  w <- loo::loo_model_weights(loos, method = method)
  stats::setNames(as.numeric(w), names(fits))
}

# exact leave-one-site-out log predictive density of site i
.elpd_without_site <- function(x, i, seed, nodes = NULL) {
  d <- x$data
  id <- d$sites[[d$site]][i]
  keep <- function(tb) if (is.null(tb)) NULL else tb[tb[[d$site]] != id, , drop = FALSE]
  d_i <- withCallingHandlers(
    fp_data(keep(d$auto), keep(d$sites), keep(d$verified), site = d$site,
            y = d$y, recorder = d$recorder, confirmed = d$confirmed),
    crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
  # fits saved before crywolf 1.5.0 lack the sampler settings: use defaults
  st <- purrr::compact(if (is.null(x$sampler)) list() else x$sampler)
  refit <- suppressWarnings(do.call(fp_occu, c(list(
    d_i, occ = x$formulas$occ, det = x$formulas$det, fp = x$formulas$fp,
    random_recorder = x$random_recorder, random_fp = isTRUE(x$random_fp),
    random_det = isTRUE(x$random_det),
    ver = if (is.null(x$formulas$ver)) ~ 1 else x$formulas$ver, markov = x$markov,
    priors = x$priors, seed = seed, refresh = 0, show_messages = FALSE), st)))
  # the left-out site must not carry a recorder or a factor level of its own,
  # or its prediction would need a parameter the refit lacks
  same <- intersect(c("recorders", "occ_names", "det_names", "fp_names", "ver_names"),
                    names(x$design))
  if (!identical(refit$design[same], x$design[same])) {
    warning("Site ", id, " is the only site with its recorder or with a level ",
            "of a factor covariate; its approximate value is kept.", call. = FALSE)
    return(NA_real_)
  }
  sd <- x$stan_data
  pars <- .ppc_pars(refit$fit, refit$stan_data, NULL)
  ll <- if (!is.null(pars$e_fp) || !is.null(pars$d_det)) {
    # the left-out site's effects are integrated over the fitted distribution
    .integrate_site_effect(sd, pars, nodes)[, i]
  } else {
    .ppc_core(sd, pars, "automated", "freeman-tukey", ll_only = TRUE)$ll[, i]
  }
  .log_mean_exp(ll)
}

# Log-likelihood of each site (S x nsite) with its random effects integrated
# out: log of the integral of p(y_i | e, d) over e ~ Normal(0, sigma_fp) (on
# the false-positive rate) and d ~ Normal(0, sigma_det) (on true detection),
# whichever the model has, for each draw, by Gauss-Hermite quadrature with a
# product grid of `nodes` points per effect. The terms are summed as they
# come (a running log-sum-exp), so memory does not grow with the grid.
.integrate_site_effect <- function(sd, pars, nodes = NULL) {
  effects <- c(e_fp = "sigma_fp", d_det = "sigma_det")
  effects <- effects[!purrr::map_lgl(pars[effects], is.null)]
  if (is.null(nodes)) nodes <- if (length(effects) > 1) 30 else 100
  gh <- .gauss_hermite(nodes)
  grid <- as.matrix(expand.grid(rep(list(seq_len(nodes)), length(effects))))
  S <- nrow(pars$beta)
  term <- function(g) {
    k <- grid[g, ]
    pars[names(effects)] <- purrr::map2(unname(effects), k, function(sigma, kj) {
      matrix(sqrt(2) * pars[[sigma]] * gh$x[kj], S, sd$nsite)
    })
    .ppc_core(sd, pars, "automated", "freeman-tukey", ll_only = TRUE)$ll +
      sum(log(gh$w[k])) - 0.5 * length(effects) * log(pi)
  }
  acc <- purrr::reduce(seq_len(nrow(grid))[-1], function(acc, g) {
    t <- term(g)
    m <- pmax(acc$m, t)
    list(m = m, s = acc$s * exp(acc$m - m) + exp(t - m))
  }, .init = list(m = term(1), s = 1))
  acc$m + log(acc$s)
}

# nodes and weights of Gauss-Hermite quadrature (weight function exp(-x^2)),
# from the eigen decomposition of the Jacobi matrix (Golub and Welsch)
.gauss_hermite <- function(n) {
  stopifnot(n >= 2)
  off <- sqrt(seq_len(n - 1) / 2)
  J <- matrix(0, n, n)
  J[cbind(seq_len(n - 1), 2:n)] <- off
  J[cbind(2:n, seq_len(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(x = e$values, w = sqrt(pi) * e$vectors[1, ]^2)
}

.log_mean_exp <- function(v) max(v) + log(mean(exp(v - max(v))))

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
