#' Simulate data for a false-positive occupancy model
#'
#' Generates one dataset from known parameters under the multiple-methods
#' design: an automated method that returns true detections and false
#' positives at every site, and a verified method, free of false positives, at
#' a random subset of sites. The false-positive rate can be homogeneous, vary
#' at random among sites (unstructured), or track a latent gradient correlated
#' with the occupancy covariate (structured confusion, as when a soundalike
#' species becomes commoner along the same gradient as the target).
#'
#' @param nsite Number of sites.
#' @param nocc Number of automated occasions per site.
#' @param psi Occupancy probability at the mean of the covariate.
#' @param b1 Occupancy slope on the standardized covariate `x`.
#' @param p11 True-detection probability of the automated method at the mean
#'   of the detection covariates.
#' @param p10 False-positive probability of the automated method (its base
#'   value; the realized mean over sites is returned as `truth$p10_marginal`).
#' @param sigma_p10 Among-site standard deviation of the logit of the
#'   false-positive rate (unstructured, logit-normal heterogeneity). 0 for none.
#' @param p10_conc Optional concentration of a Beta distribution for the
#'   site-level false-positive rate, with mean `p10` (unstructured,
#'   aggregated heterogeneity; small values push rates towards 0 or 1).
#'   Overrides `sigma_p10`.
#' @param conf_effect Slope of the logit of the false-positive rate on the
#'   confusing gradient `w` (structured confusion). 0 for none.
#' @param conf_cor Correlation between `w` and the occupancy covariate `x`.
#' @param r11 Detection probability of the verified method.
#' @param frac_verified Probability that a site is verified.
#' @param nver Number of verified occasions per verified site.
#' @param a_bio,a_aci Effects of two standardized observation covariates
#'   (`bio`, `aci`) on true detection. 0 by default.
#' @param sd_recorder Standard deviation of a recorder random intercept on
#'   true detection. 0 by default.
#' @param nrec Number of recorders.
#' @param frac_confirmed Probability that an automated occasion is checked
#'   under the observation-confirmation design of Chambert et al. (2015). 0
#'   (default) for none. See Details.
#' @param seed Optional integer seed.
#'
#' @details
#' The occupancy covariate `x` is standardized, so `b1` is the true slope
#' exactly. The confusing gradient is \eqn{w = \rho x + \sqrt{1 - \rho^2}
#' \epsilon}, with \eqn{\epsilon} standard normal, so it has unit variance and
#' correlation \eqn{\rho} with `x`. In practice `w` is not observed; it is
#' returned in `sites` so that a model given the true structure can serve as a
#' positive control. With the same seed, the function returns the same
#' datasets as `simulate_fp()` of the fpcalib research compendium (v0.1.0,
#' \doi{10.5281/zenodo.23379100}), with `nocc`, `frac_verified`, `nver`,
#' `p10_conc`, `sd_recorder`, and `nrec` corresponding to its `K`,
#' `verified_coverage`, `n_verified`, `p10_beta_conc`, `sig_u`, and `nobser`
#' (checked by `dev/compare-simulate-fp.R` in the package repository).
#'
#' With `frac_confirmed > 0`, false observations occur at occupied sites too:
#' on each occasion, a true observation occurs with probability `p11` (and
#' its covariate effects) at occupied sites, and a false one with the site's
#' false-positive rate at any site, independently. The occasion is a
#' detection if either occurs. Checked occasions get the columns `n_true`
#' and `n_false` of `auto` (0 or 1; `NA` elsewhere). These draws come after
#' all others, so the default design reproduces the same random numbers.
#'
#' @return A list of class `fp_sim` with elements `data` (an `fp_data` object
#'   whose `sites` table holds `x` and `w`, and whose `auto` table holds `bio`,
#'   `aci`, and `recorder`) and `truth` (the generating values, the realized
#'   mean false-positive rate `p10_marginal`, the realized correlation
#'   between `w` and `x`, and the latent states `z`).
#' @seealso [fp_occu()], [fp_calibrate()]
#' @export
#' @examples
#' # homogeneous false positives
#' s0 <- simulate_fp(nsite = 100, nocc = 10, seed = 1)
#' s0
#'
#' # structured confusion
#' s1 <- simulate_fp(nsite = 100, nocc = 10, conf_effect = 1.5,
#'                   conf_cor = 0.7, seed = 1)
#' s1$truth$p10_marginal
simulate_fp <- function(nsite = 300, nocc = 20, psi = 0.5, b1 = -0.8,
                        p11 = 0.6, p10 = 0.05,
                        sigma_p10 = 0, p10_conc = NULL,
                        conf_effect = 0, conf_cor = 0,
                        r11 = 0.6, frac_verified = 0.3, nver = 3,
                        a_bio = 0, a_aci = 0, sd_recorder = 0, nrec = 5,
                        frac_confirmed = 0, seed = NULL) {
  stopifnot(nsite >= 2, nocc >= 1, nver >= 1,
            psi > 0, psi < 1, p11 > 0, p11 < 1, p10 > 0, p10 < 1,
            r11 > 0, r11 <= 1, frac_verified >= 0, frac_verified <= 1,
            abs(conf_cor) <= 1, sigma_p10 >= 0, sd_recorder >= 0,
            frac_confirmed >= 0, frac_confirmed <= 1)
  if (!is.null(seed)) .set_seed(seed)
  b0 <- stats::qlogis(psi)
  a0 <- stats::qlogis(p11)

  # occupancy
  x_raw <- stats::rnorm(nsite)
  x <- as.numeric(scale(x_raw))
  z <- stats::rbinom(nsite, 1, stats::plogis(b0 + b1 * x))

  # recorders
  rec <- sample(seq_len(nrec), nsite, replace = TRUE)
  u <- stats::rnorm(nrec, 0, sd_recorder)

  # automated occasions (site x occasion matrices)
  BIO <- matrix(stats::rnorm(nsite * nocc), nsite, nocc)
  ACI <- matrix(stats::rnorm(nsite * nocc), nsite, nocc)
  bio <- (BIO - mean(BIO)) / stats::sd(BIO)
  aci <- (ACI - mean(ACI)) / stats::sd(ACI)
  P11 <- stats::plogis(a0 + a_bio * bio + a_aci * aci + u[rec])

  # site-level false-positive rate
  clampl <- function(p) stats::qlogis(pmin(pmax(p, 1e-6), 1 - 1e-6))
  if (!is.null(p10_conc)) {
    p10_base <- stats::rbeta(nsite, p10 * p10_conc, (1 - p10) * p10_conc)
  } else if (sigma_p10 > 0) {
    p10_base <- stats::plogis(stats::qlogis(p10) + stats::rnorm(nsite, 0, sigma_p10))
  } else {
    p10_base <- rep(p10, nsite)
  }
  if (conf_effect != 0) {
    w <- conf_cor * x + sqrt(1 - conf_cor^2) * stats::rnorm(nsite)
    p10_i <- stats::plogis(clampl(p10_base) + conf_effect * w)
  } else {
    w <- rep(0, nsite)
    p10_i <- p10_base
  }

  Pdet <- matrix(z, nsite, nocc) * P11 + (1 - matrix(z, nsite, nocc)) *
    matrix(p10_i, nsite, nocc)
  Y2 <- matrix(stats::rbinom(nsite * nocc, 1, Pdet), nsite, nocc)

  # verified occasions
  is_ver <- stats::rbinom(nsite, 1, frac_verified) == 1
  M1 <- matrix(NA_integer_, nsite, nver)
  if (any(is_ver)) {
    nv <- sum(is_ver)
    pv <- ifelse(z[is_ver] == 1, r11, 0)
    M1[is_ver, ] <- matrix(stats::rbinom(nv * nver, 1, matrix(pv, nv, nver)),
                           nv, nver)
  }

  # observation confirmation (Chambert et al. 2015)
  if (frac_confirmed > 0) {
    Zm <- matrix(z, nsite, nocc)
    extra <- matrix(stats::rbinom(nsite * nocc, 1, matrix(p10_i, nsite, nocc)),
                    nsite, nocc)
    true_obs <- Zm * Y2
    false_obs <- (1L - Zm) * Y2 + Zm * extra
    Y2 <- pmax(true_obs, false_obs)
    checked <- stats::runif(nsite * nocc) < frac_confirmed
  }

  # long tables (column-major, as in the compendium)
  auto <- data.frame(
    site = rep(seq_len(nsite), times = nocc),
    occasion = rep(seq_len(nocc), each = nsite),
    y = as.vector(Y2),
    bio = as.vector(bio), aci = as.vector(aci),
    recorder = rep(rec, times = nocc)
  )
  if (frac_confirmed > 0) {
    auto$n_true <- ifelse(checked, as.vector(true_obs), NA_integer_)
    auto$n_false <- ifelse(checked, as.vector(false_obs), NA_integer_)
  }
  verified <- data.frame(
    site = rep(seq_len(nsite), times = nver),
    occasion = rep(seq_len(nver), each = nsite),
    y = as.vector(M1)
  )
  verified <- verified[!is.na(verified$y), , drop = FALSE]
  sites <- data.frame(site = seq_len(nsite), x = x, w = w)

  data <- fp_data(auto = auto, sites = sites, verified = verified,
                  recorder = "recorder",
                  confirmed = if (frac_confirmed > 0) c("n_true", "n_false"))

  regime <- if (conf_effect != 0) "structured"
            else if (!is.null(p10_conc) || sigma_p10 > 0) "unstructured"
            else "homogeneous"
  truth <- list(
    regime = regime, psi = psi, b0 = b0, b1 = b1, p11 = p11, p10 = p10,
    r11 = r11, p10_marginal = mean(p10_i), sigma_p10 = sigma_p10,
    p10_conc = if (is.null(p10_conc)) NA_real_ else p10_conc,
    conf_effect = conf_effect, conf_cor = conf_cor,
    conf_cor_realized = if (conf_effect != 0) stats::cor(w, x) else NA_real_,
    a_bio = a_bio, a_aci = a_aci, sd_recorder = sd_recorder,
    frac_confirmed = frac_confirmed,
    n_verified_sites = sum(is_ver), z = z, p10_site = p10_i
  )
  structure(list(data = data, truth = truth), class = "fp_sim")
}

#' @export
print.fp_sim <- function(x, ...) {
  t <- x$truth
  cat(sprintf("<fp_sim> %s false positives\n", t$regime))
  cat(sprintf("  truth: psi = %.2f, b1 = %.2f, p11 = %.2f, p10 = %.3f (realized mean %.3f), r11 = %.2f\n",
              t$psi, t$b1, t$p11, t$p10, t$p10_marginal, t$r11))
  if (t$regime == "structured") {
    cat(sprintf("  confusion: effect %.2f, correlation with x %.2f (realized %.2f)\n",
                t$conf_effect, t$conf_cor, t$conf_cor_realized))
  }
  cat(sprintf("  realized occupancy: %.2f\n", mean(t$z)))
  print(x$data)
  invisible(x)
}

# set.seed() with R's default generators whatever kind is active: the
# workers of 'furrr' use L'Ecuyer-CMRG, under which the same seed would give
# other draws than in a sequential run (and than the fpcalib compendium).
# The caller's generator kinds are restored when the caller returns.
.set_seed <- function(seed, envir = parent.frame()) {
  default <- c("Mersenne-Twister", "Inversion", "Rejection")
  old <- RNGkind()
  set.seed(seed, kind = default[1], normal.kind = default[2],
           sample.kind = default[3])
  if (!identical(old, default)) {
    restore <- call("RNGkind", old[1], old[2], old[3])
    do.call(base::on.exit, list(restore, add = TRUE), envir = envir)
  }
  invisible(NULL)
}
