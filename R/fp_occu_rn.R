#' Fit a Royle-Nichols occupancy model with false positives
#'
#' @description
#' Fits the abundance-induced heterogeneity model of Royle and Nichols (2003)
#' with false positives in Stan. The number of individuals at site \eqn{i} is
#' \eqn{N_i \sim \mathrm{Poisson}(\lambda_i)}, each detected on an occasion
#' with probability \eqn{r}, so the automated method detects the species with
#' probability \eqn{1 - (1 - r)^{N_i}}, and the site is occupied when
#' \eqn{N_i > 0}. A site effect on \eqn{r} (`random_det`) then measures
#' variation in detection among sites beyond that due to abundance, which
#' [fp_occu()] cannot separate from it.
#'
#' @inheritParams fp_occu
#' @param occ One-sided formula for log abundance (\eqn{\log \lambda}); site
#'   covariates only. The occupancy probability is
#'   \eqn{\psi = 1 - e^{-\lambda}}.
#' @param det One-sided formula for the logit of the per-individual detection
#'   \eqn{r} of the automated method (site or occasion covariates).
#' @param fp One-sided formula for the logit of the false-positive rate, or
#'   `NULL` for the Royle-Nichols model without false positives. Site
#'   covariates only.
#' @param ver One-sided formula for the logit of the per-individual detection
#'   of the verified method.
#' @param random_fp Logical. Add a site random intercept to the logit of the
#'   false-positive rate.
#' @param random_det Logical. Add a site random intercept to the logit of the
#'   per-individual detection \eqn{r}: its standard deviation (`sd_det_site`)
#'   is the variation in detection among sites that abundance does not
#'   explain.
#' @param random_abund Logical. Add a site random intercept to log abundance
#'   (`sd_abund_site`), for variation in abundance beyond the covariates of
#'   `occ`.
#' @param n_max Largest abundance summed over. The fit warns when the
#'   Poisson mass above it is not negligible at some site.
#'
#' @details
#' **Model.** With \eqn{N_i \sim \mathrm{Poisson}(\lambda_i)},
#' \deqn{\log \lambda_i = X_i \beta \; [+ v_i],}
#' an automated occasion \eqn{k} of site \eqn{i} is a detection with
#' probability \deqn{1 - (1 - r_{ik})^{N_i} (1 - p_{10,i}),}
#' \deqn{\mathrm{logit}(r_{ik}) = Z_{ik} \alpha + u_{\mathrm{rec}[ik]} \; [+ d_i],}
#' \deqn{\mathrm{logit}(p_{10,i}) = \gamma_0 + X^{fp}_i \gamma \; [+ e_i],}
#' so false positives occur on every occasion, also at occupied sites, and an
#' unoccupied site (\eqn{N_i = 0}) is detected with probability
#' \eqn{p_{10,i}}. A verified occasion is a detection with probability
#' \eqn{1 - (1 - r^{v})^{N_i}}. The site effects \eqn{v_i}, \eqn{d_i} and
#' \eqn{e_i} are normal with standard deviations `sd_abund_site`,
#' `sd_det_site` and `sd_fp_site` and half-normal priors (`abund_site_sd`,
#' `det_site_sd` and `fp_site_sd` of [fp_priors()]). Unlike [fp_occu()], the
#' false-positive intercept \eqn{\gamma_0} is not tied to the detection
#' intercept by a gap: detection at an occupied site always exceeds
#' \eqn{p_{10}}, so the labels cannot switch, and \eqn{p_{10}} may exceed
#' the per-individual detection \eqn{r}. Its prior is the marginal prior the
#' gap implies in [fp_occu()]: Normal with mean `-gap_mean` and standard
#' deviation `sqrt(intercept_sd^2 + gap_sd^2)`. \eqn{N_i} is summed out
#' over 0 to `n_max`.
#'
#' **Observation-confirmation design.** Confirmed occasions (argument
#' `confirmed` of [fp_data()]) enter as in [fp_occu()], with probabilities of
#' at least one true observation \eqn{s_1 = 1 - (1 - r)^{N_i}} and at least
#' one false observation \eqn{s_0 = p_{10}}. They can be the same occasions
#' as the automated detections, so a fully annotated survey enters as one
#' row per occasion with both columns filled.
#'
#' **Not supported.** False-positive covariates measured on each occasion,
#' runs of detections (`markov` of [fp_occu()]), the outputs of [fp_occu()]
#' other than `print()`, `summary()` and `as_draws_df()` (for now).
#'
#' **Identifiability.** Simulations at 200 sites and at 21 sites with 17
#' occasions (in `dev/rn-fp/` of the package repository) recovered
#' `sd_det_site` near 1 and gave few false signals when it was 0; values
#' near 0.5 were not distinguishable from 0, and `sd_abund_site` was weakly
#' identified. With few sites the prior of `sd_det_site` carries much of the
#' posterior near 0: compare fits under several values of `det_site_sd`.
#'
#' @return An object of class `fp_occu_rn`, with the elements of an
#'   [fp_occu()] fit (`fit`, `draws`, `stan_data`, `formulas`, `design`,
#'   `data`, `priors`, `sampler`, `diagnostics`, `call`) and `n_max`. The
#'   draws hold the abundance coefficients (`abund_*`), the detection
#'   coefficients on the logit of \eqn{r} (`det_*`), the false-positive
#'   intercept and slopes, the verified detection, the derived
#'   quantities `lambda_ref` and `psi_ref` (at the reference of the
#'   covariates), `lambda_bar` and `psi_bar` (means over sites), `r_ref`,
#'   `p10_ref` and `r_verified`, and the standard deviations of the site
#'   effects.
#' @references
#' Royle, J. A. and Nichols, J. D. (2003). Estimating abundance from repeated
#' presence-absence data or point counts. *Ecology*, 84, 777–790.
#' \doi{10.1890/0012-9658(2003)084[0777:EAFRPA]2.0.CO;2}
#'
#' Chambert, T., Miller, D. A. W. and Nichols, J. D. (2015). Modeling false
#' positive detections in species occurrence data under different study
#' designs. *Ecology*, 96, 332–339. \doi{10.1890/14-1507.1}
#' @section Experimental:
#' This function is experimental (new after 2.3.0): its arguments and
#' results may change in a minor version. See [crywolf-api].
#' @seealso [simulate_fp_rn()], [fp_occu()]
#' @export
#' @examplesIf instantiate::stan_cmdstan_exists()
#' \donttest{
#' s <- simulate_fp_rn(nsite = 60, nocc = 10, sigma_r = 1, seed = 1)
#' m <- fp_occu_rn(s$data, random_det = TRUE, chains = 2,
#'                 iter_warmup = 300, iter_sampling = 300, seed = 1)
#' m
#' }
fp_occu_rn <- function(data, occ = ~ 1, det = ~ 1, fp = ~ 1,
                       random_recorder = !is.null(data$recorder),
                       ver = ~ 1, random_fp = FALSE, random_det = FALSE,
                       random_abund = FALSE, n_max = 20,
                       priors = fp_priors(), chains = 4,
                       parallel_chains = chains, threads_per_chain = "auto",
                       iter_warmup = 1000, iter_sampling = 1000,
                       adapt_delta = 0.95, seed = NULL, refresh = NULL,
                       prior_only = FALSE, ...) {
  .check_fp_data(data)
  .check_sampler(chains, iter_warmup, iter_sampling, adapt_delta, seed)
  .check_flag(random_abund, "random_abund")
  .check_flag(prior_only, "prior_only")
  if (!is.numeric(n_max) || length(n_max) != 1 || !is.finite(n_max) ||
      n_max < 1 || n_max != round(n_max)) {
    stop("`n_max` must be a single positive whole number.", call. = FALSE)
  }
  no_ver <- is.null(data$verified) || nrow(data$verified) == 0
  no_conf <- is.null(data$confirmed) || all(is.na(data$auto[[data$confirmed[1]]]))
  if (!is.null(fp) && no_ver && no_conf) {
    warning("No verified data: the false-positive rate is identified only by ",
            "the priors. Treat the estimates with caution.", call. = FALSE)
  }
  cores <- .sampler_cores(threads_per_chain, parallel_chains)
  sd_list <- .build_stan_data(data, occ, det, fp, random_recorder, priors,
                              random_fp = random_fp, ver = ver,
                              random_det = random_det)
  sd <- sd_list$stan_data
  if (sd$fp_by_row == 1) {
    stop("`fp` uses occasion-level covariates, which fp_occu_rn() does not ",
         "support; use site covariates only.", call. = FALSE)
  }
  sd <- .collapse_confirmed(sd)
  sd$Nmax <- as.integer(n_max)
  sd$estimate_abund_re <- as.integer(random_abund)
  sd$abund_site_sd <- if (is.null(priors$abund_site_sd)) 1 else priors$abund_site_sd
  sd$prior_only <- as.integer(prior_only)
  # fields of fp_occu.stan that this program does not declare
  sd[c("fp_by_row", "Xfp2", "Xfp3")] <- NULL
  sd_list$stan_data <- sd
  if (is.null(refresh)) {
    refresh <- if (interactive()) max(1, round((iter_warmup + iter_sampling) / 10)) else 0
  }
  mod <- .compile_model("fp_occu_rn")
  # as in fp_occu(): a fit as `init`, and output file names of our own
  # unless the user gives some in `...`
  run_sampler <- function(..., init = NULL, output_basename = .output_basename()) {
    mod$sample(..., init = .init_from_fit(init, chains),
               output_basename = output_basename)
  }
  fit <- run_sampler(data = sd, chains = chains,
                     parallel_chains = cores$parallel_chains,
                     threads_per_chain = cores$threads_per_chain,
                     iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                     adapt_delta = adapt_delta, seed = seed, refresh = refresh, ...)
  fit$draws()
  fit$sampler_diagnostics()
  draws <- .pretty_draws_rn(fit, sd_list)
  diag <- .diagnostics(fit, draws)
  p_above <- as.vector(fit$draws("p_above_nmax"))
  out <- structure(
    list(fit = fit, draws = draws, stan_data = sd,
         formulas = list(occ = occ, det = det, fp = fp, ver = ver),
         design = sd_list$design, data = data,
         random_recorder = isTRUE(sd$estimate_re == 1),
         random_fp = isTRUE(sd$estimate_fp_re == 1),
         random_det = isTRUE(sd$estimate_det_re == 1),
         random_abund = isTRUE(random_abund), n_max = as.integer(n_max),
         priors = priors,
         sampler = list(chains = chains, parallel_chains = cores$parallel_chains,
                        threads_per_chain = cores$threads_per_chain,
                        iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                        adapt_delta = adapt_delta, seed = fit$metadata()$seed),
         prior_only = prior_only, diagnostics = diag, call = match.call()),
    class = "fp_occu_rn"
  )
  .warn_diagnostics(diag)
  if (!prior_only && stats::median(p_above) > 0.01) {
    warning(sprintf(paste0("At some site the posterior median probability ",
                           "that abundance exceeds `n_max` (%d) is %.2f; ",
                           "raise `n_max`."), n_max, stats::median(p_above)),
            call. = FALSE)
  }
  out
}

#' Confirmed occasions of a site with the same recorder, detection
#' covariates and state collapse into one row of w3 occasions (exact, since
#' their terms are identical); they stay sorted by site
#' @noRd
.collapse_confirmed <- function(sd) {
  if (sd$N3 == 0) {
    sd$w3 <- integer(0)
    return(sd)
  }
  key <- do.call(paste, c(list(sd$site3, sd$rec3, sd$v3), as.data.frame(sd$Z3),
                          list(sep = "\r")))
  grp <- match(key, unique(key))
  first <- !duplicated(grp)
  sd$w3 <- tabulate(grp, nbins = sum(first))
  sd$v3 <- sd$v3[first]
  sd$site3 <- sd$site3[first]
  sd$rec3 <- sd$rec3[first]
  sd$Z3 <- sd$Z3[first, , drop = FALSE]
  sd$N3 <- sum(first)
  sd
}

#' Map from the Stan names of fp_occu_rn.stan to readable names
#' @noRd
.name_map_rn <- function(design, sd) {
  m <- c(
    stats::setNames(paste0("abund_", design$occ_names),
                    sprintf("beta[%d]", seq_along(design$occ_names))),
    stats::setNames(paste0("det_", design$det_names),
                    sprintf("alpha[%d]", seq_along(design$det_names))),
    if (sd$N1 > 0) c("delta0" = "verified_(Intercept)", "r_verified" = "r_verified"),
    if (sd$N1 > 0 && length(design$ver_names) > 1) {
      stats::setNames(paste0("verified_", design$ver_names[-1]),
                      sprintf("delta_s[%d]", seq_along(design$ver_names[-1])))
    },
    c("lambda_ref" = "lambda_ref", "psi_ref" = "psi_ref",
      "lambda_bar" = "lambda_bar", "psi_bar" = "psi_bar", "r_ref" = "r_ref")
  )
  if (sd$estimate_fp == 1) {
    slopes <- design$fp_names[-1]
    m <- c(m, c("gamma0" = "fp_(Intercept)"),
           if (length(slopes)) {
             stats::setNames(paste0("fp_", slopes),
                             sprintf("gamma_fp[%d]", seq_along(slopes)))
           },
           c("p10_ref" = "p10_ref"))
  }
  if (sd$estimate_re == 1) m <- c(m, c("sigma_u" = "sd_recorder"))
  if (sd$estimate_fp_re == 1) m <- c(m, c("sigma_fp[1]" = "sd_fp_site"))
  if (sd$estimate_det_re == 1) m <- c(m, c("sigma_det[1]" = "sd_det_site"))
  if (sd$estimate_abund_re == 1) m <- c(m, c("sigma_lam[1]" = "sd_abund_site"))
  m
}

#' Posterior draws of fp_occu_rn() with readable names
#' @noRd
.pretty_draws_rn <- function(fit, sd_list) {
  map <- .name_map_rn(sd_list$design, sd_list$stan_data)
  dr <- posterior::as_draws_df(fit$draws(variables = unique(sub("\\[.*", "", names(map)))))
  dr <- posterior::subset_draws(dr, variable = names(map))
  posterior::variables(dr) <- unname(map[posterior::variables(dr)])
  dr
}

#' Methods for Royle-Nichols fits
#'
#' @description
#' `print()`, `summary()` and `as_draws_df()` for fits of [fp_occu_rn()].
#'
#' @param x,object A fit of [fp_occu_rn()].
#' @param digits Digits to print.
#' @param prob Probability mass of the equal-tailed credible intervals.
#' @param ... Not used.
#' @return `print()` returns `x` invisibly; `summary()` an object of class
#'   `summary.fp_occu` (see [summary.fp_occu()]); `as_draws_df()` the
#'   posterior draws with readable names.
#' @name fp_occu_rn-methods
NULL

#' @rdname fp_occu_rn-methods
#' @export
print.fp_occu_rn <- function(x, digits = 3, ...) {
  sd <- x$stan_data
  cat("Bayesian Royle-Nichols occupancy model with ",
      if (isTRUE(sd$estimate_fp == 1)) "false positives" else "no false positives",
      if (isTRUE(sd$obs_confirm == 1)) " (observation-confirmation design)", "\n\n",
      sep = "")
  if (isTRUE(x$prior_only)) {
    cat("  Draws from the priors alone (prior_only = TRUE), not a fit to the data.\n\n")
  }
  f <- x$formulas
  cat("  log abundance:   ", .fmt_formula(f$occ),
      if (x$random_abund) "+ (1 | site)" else "", "\n")
  cat("  per-individual detection: ", .fmt_formula(f$det),
      if (x$random_recorder) "+ (1 | recorder)" else "",
      if (x$random_det) "+ (1 | site)" else "", "\n")
  cat("  false positive:  ", if (is.null(f$fp)) "none" else .fmt_formula(f$fp),
      if (x$random_fp) "+ (1 | site)" else "", "\n")
  if (!is.null(f$ver) && length(all.vars(f$ver))) {
    cat("  verified:        ", .fmt_formula(f$ver), "\n")
  }
  cat(sprintf("  data:             %d sites, %d automated (%d confirmed) and %d verified occasions; N summed over 0..%d\n\n",
              sd$nsite, sum(sd$n2) + sum(sd$w3), sum(sd$w3), sd$N1, x$n_max))
  print(as.data.frame(summary(x)$table), digits = digits, row.names = FALSE)
  .print_diag(x$diagnostics)
  invisible(x)
}

#' @rdname fp_occu_rn-methods
#' @export
summary.fp_occu_rn <- function(object, prob = 0.95, ...) {
  summary.fp_occu(object, prob = prob, ...)
}

#' @rdname fp_occu_rn-methods
#' @export
as_draws_df.fp_occu_rn <- function(x, ...) {
  x$draws
}
