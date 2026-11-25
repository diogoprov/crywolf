#' Simulate data for a Royle-Nichols model with false positives
#'
#' @description
#' Generates one dataset from the model of [fp_occu_rn()]: Poisson abundance
#' with a covariate and optional site effect, per-individual detection with an
#' optional site effect, false positives on every automated occasion, verified
#' occasions at a share of the sites, and, optionally, confirmed automated
#' occasions (observation-confirmation design).
#'
#' @param nsite Number of sites.
#' @param nocc Number of automated occasions per site: one number, or one per
#'   site (unbalanced designs).
#' @param psi Occupancy probability at the mean of the covariate `x`; the
#'   intercept of log abundance is \eqn{\log(-\log(1 - \psi))}.
#' @param b1 Slope of log abundance on the standardized covariate `x`.
#' @param r Per-individual detection probability of the automated method.
#' @param p10 False-positive probability of the automated method.
#' @param sigma_r Standard deviation of a site effect on the logit of `r`.
#' @param sigma_lam Standard deviation of a site effect on log abundance.
#' @param rv Per-individual detection probability of the verified method.
#' @param frac_verified Probability that a site has verified occasions.
#' @param nver Number of verified occasions per verified site.
#' @param frac_confirmed Probability that an automated occasion is confirmed;
#'   1 confirms every occasion, as in a fully annotated survey.
#' @param seed Optional integer seed.
#'
#' @details
#' On each automated occasion a true observation occurs with probability
#' \eqn{1 - (1 - r_i)^{N_i}} and a false one with probability `p10`,
#' independently; the occasion is a detection if either occurs. Confirmed
#' occasions get the columns `n_true` and `n_false` of `auto` (0 or 1, `NA`
#' elsewhere). A verified occasion is a detection with probability
#' \eqn{1 - (1 - r^v)^{N_i}}.
#'
#' @return A list of class `fp_sim` with elements `data` (an `fp_data`
#'   object; `sites` holds `x`) and `truth` (the generating values, the
#'   intercepts `b0` of log abundance and `a0` of the logit of `r`, the
#'   mean true-detection probability at occupied sites `p11`, the
#'   abundances `N`, the states `z` and the site effects).
#' @seealso [fp_occu_rn()], [simulate_fp()]
#' @export
#' @examples
#' s <- simulate_fp_rn(nsite = 50, nocc = 10, sigma_r = 1, seed = 1)
#' s
simulate_fp_rn <- function(nsite = 200, nocc = 10, psi = 0.6, b1 = 0.5,
                           r = 0.3, p10 = 0.05, sigma_r = 0, sigma_lam = 0,
                           rv = 0.5, frac_verified = 0.3, nver = 3,
                           frac_confirmed = 0, seed = NULL) {
  args <- list(nsite = nsite, psi = psi, b1 = b1, r = r, p10 = p10,
               sigma_r = sigma_r, sigma_lam = sigma_lam, rv = rv,
               frac_verified = frac_verified, nver = nver,
               frac_confirmed = frac_confirmed)
  bad <- names(args)[!vapply(args, function(v) is.numeric(v) && length(v) == 1 &&
                                       is.finite(v), logical(1))]
  if (length(bad)) {
    stop("Argument(s) ", paste0("`", bad, "`", collapse = ", "),
         " must be single finite numbers.", call. = FALSE)
  }
  if (!is.numeric(nocc) || !(length(nocc) %in% c(1, nsite)) || anyNA(nocc) ||
      any(nocc < 1) || any(nocc != round(nocc))) {
    stop("`nocc` must be one positive whole number or one per site.", call. = FALSE)
  }
  if (nver != round(nver)) stop("`nver` must be a whole number.", call. = FALSE)
  stopifnot(nsite >= 2, nver >= 1, psi > 0, psi < 1, r > 0, r < 1,
            p10 > 0, p10 < 1, rv > 0, rv <= 1, sigma_r >= 0, sigma_lam >= 0,
            frac_verified >= 0, frac_verified <= 1,
            frac_confirmed >= 0, frac_confirmed <= 1)
  if (!is.null(seed)) .set_seed(seed)
  nocc <- rep_len(as.integer(nocc), nsite)
  b0 <- log(-log(1 - psi))
  a0 <- stats::qlogis(r)

  x <- as.numeric(scale(stats::rnorm(nsite)))
  v <- stats::rnorm(nsite, 0, sigma_lam)
  lambda <- exp(b0 + b1 * x + v)
  N <- stats::rpois(nsite, lambda)
  d <- stats::rnorm(nsite, 0, sigma_r)
  r_i <- stats::plogis(a0 + d)

  # automated occasions, one row per site and occasion
  site <- rep(seq_len(nsite), nocc)
  occasion <- sequence(nocc)
  n_occ <- length(site)
  true_obs <- stats::rbinom(n_occ, 1, (1 - (1 - r_i)^N)[site])
  false_obs <- stats::rbinom(n_occ, 1, p10)
  auto <- data.frame(site = site, occasion = occasion,
                     y = pmax(true_obs, false_obs))
  if (frac_confirmed > 0) {
    checked <- stats::runif(n_occ) < frac_confirmed
    auto$n_true <- ifelse(checked, true_obs, NA_integer_)
    auto$n_false <- ifelse(checked, false_obs, NA_integer_)
  }

  # verified occasions
  is_ver <- stats::rbinom(nsite, 1, frac_verified) == 1
  vs <- rep(which(is_ver), each = nver)
  verified <- data.frame(site = vs, occasion = rep(seq_len(nver), sum(is_ver)),
                         y = stats::rbinom(length(vs), 1, (1 - (1 - rv)^N)[vs]))
  sites <- data.frame(site = seq_len(nsite), x = x)

  data <- fp_data(auto = auto, sites = sites, verified = verified,
                  confirmed = if (frac_confirmed > 0) c("n_true", "n_false"))
  truth <- list(regime = "Royle-Nichols", psi = psi, b0 = b0, b1 = b1, r = r,
                a0 = a0, p10 = p10, p10_marginal = p10, sigma_r = sigma_r,
                sigma_lam = sigma_lam, rv = rv,
                # mean automated detection (no false positives) at occupied sites
                p11 = if (any(N > 0)) mean((1 - (1 - r_i)^N)[N > 0]) else NA_real_,
                r11 = rv, frac_confirmed = frac_confirmed,
                n_verified_sites = sum(is_ver), N = N, z = as.integer(N > 0),
                lambda = lambda, site_effect_r = d, site_effect_lambda = v)
  structure(list(data = data, truth = truth), class = "fp_sim")
}
