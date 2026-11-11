#' Simulate data for a multispecies false-positive occupancy model
#'
#' Generates one dataset from known parameters for [fp_occu_ms()]: several
#' species on the same sites, occasions, recorders and verified occasions,
#' each with its own occupancy, true detection and false positives drawn
#' from distributions of the community, as in the model. **In
#' development**, as [fp_occu_ms()].
#'
#' The sites, the occupancy covariate `x`, the occasions with their
#' covariates `bio` and `aci`, the recorders, and which occasions are
#' verified (or checked, with `frac_confirmed > 0`) come from
#' [simulate_fp()] and are shared by the species. Each species \eqn{s} then
#' has occupancy intercept and slope on `x` \eqn{\beta_s = \mu_\beta +
#' \sigma_\beta z_s}, true-detection intercept \eqn{\alpha_s = \alpha +
#' \sigma_\alpha z_s}, and gap \eqn{g_s = \mu_g \exp(\sigma_g z_s)} between
#' its true-detection and false-positive intercepts, with \eqn{z_s} standard
#' normal, so that its false-positive rate is \eqn{\mathrm{logit}^{-1}
#' (\alpha_s - g_s)}. The slope of true detection on `bio`, the verified
#' detection probability and the recorder effects are shared. Given these
#' values, the detections of each species are drawn as in
#' [simulate.fp_occu()].
#'
#' @param nspecies Number of species.
#' @param nsite,nocc,frac_verified,nver,nrec,frac_confirmed As in
#'   [simulate_fp()].
#' @param psi,b1 Community mean occupancy probability at the mean of `x`, and
#'   community mean slope on `x`.
#' @param sd_psi,sd_b1 Standard deviations among species of the logit of
#'   occupancy at the mean of `x` and of the slope.
#' @param p11,p10 True-detection and false-positive probabilities of the
#'   community (at the mean detection intercept and the median gap);
#'   `p10` must be below `p11`.
#' @param sd_p11 Standard deviation among species of the true-detection
#'   intercept (logit scale).
#' @param sd_gap Standard deviation among species of the log of the gap.
#' @param a_bio Shared slope of true detection on `bio`.
#' @param r11 Detection probability of the verified method, below 1.
#' @param sd_recorder Standard deviation of the recorder effects on true
#'   detection, shared by the species. 0 by default.
#' @param seed Optional integer seed.
#'
#' @return A list of class `fp_sim_ms` with `data` (an `fp_data` object with
#'   the species column `species`; species are named `sp1`, `sp2`, ...) and
#'   `truth`: `pars`, the generating values named as the output of
#'   [fp_occu_ms()] with `occ = ~ x` and `det = ~ bio` (community values and
#'   the values of each species, `occ_x[sp1]`, ...); `z`, the latent states,
#'   a matrix of sites by species; and `p10_bar`, the false-positive rate of
#'   each species.
#' @seealso [fp_occu_ms()], [simulate_fp()]
#' @export
#' @examples
#' s <- simulate_fp_ms(nspecies = 4, nsite = 50, nocc = 8, seed = 1)
#' s
simulate_fp_ms <- function(nspecies = 5, nsite = 100, nocc = 10,
                           psi = 0.4, sd_psi = 1, b1 = -0.5, sd_b1 = 0.5,
                           p11 = 0.5, sd_p11 = 0.5, p10 = 0.05, sd_gap = 0.3,
                           a_bio = 0, r11 = 0.6, frac_verified = 0.3, nver = 3,
                           sd_recorder = 0, nrec = 5, frac_confirmed = 0, seed = NULL) {
  stopifnot(nspecies >= 1, psi > 0, psi < 1, p11 > 0, p11 < 1, p10 > 0, p10 < p11,
            r11 > 0, r11 < 1, sd_psi >= 0, sd_b1 >= 0, sd_p11 >= 0, sd_gap >= 0,
            sd_recorder >= 0)
  if (!is.null(seed)) .set_seed(seed)
  # the shared design; its detections are replaced below
  template <- simulate_fp(nsite = nsite, nocc = nocc, frac_verified = frac_verified,
                          nver = nver, nrec = nrec, frac_confirmed = frac_confirmed)$data
  use_re <- sd_recorder > 0
  sl <- .build_stan_data(template, ~ x, ~ bio, ~ 1, use_re, fp_priors())
  gen <- list(data = template, stan_data = sl$stan_data, design = sl$design,
              occasions = sl$occasions)

  mu_beta <- c(stats::qlogis(psi), b1)
  sigma_beta <- c(sd_psi, sd_b1)
  alpha0 <- stats::qlogis(p11)
  mu_gap <- alpha0 - stats::qlogis(p10)
  delta0 <- stats::qlogis(r11)
  u <- if (use_re) stats::rnorm(length(sl$design$recorders), 0, sd_recorder)
  sp <- sprintf("sp%d", seq_len(nspecies))
  beta <- matrix(mu_beta, 2, nspecies) + sigma_beta * matrix(stats::rnorm(2 * nspecies), 2)
  alpha <- alpha0 + sd_p11 * stats::rnorm(nspecies)
  gap <- mu_gap * exp(sd_gap * stats::rnorm(nspecies))

  res <- purrr::map(seq_len(nspecies), function(s) {
    p <- list(beta = beta[, s], alpha = c(alpha[s], a_bio), gamma0 = alpha[s] - gap[s],
              gamma_fp = numeric(0), delta0 = delta0, delta_s = numeric(0), u = u)
    r <- .simulate_from(gen, p)
    tag <- function(df) if (is.null(df)) NULL else transform(df, species = sp[s])
    r$data$auto <- tag(r$data$auto)
    r$data$verified <- tag(r$data$verified)
    r
  })
  names(res) <- sp
  sims <- .ms_combine_sims(res, utils::modifyList(template, list(species = "species")))

  occ <- c("occ_(Intercept)", "occ_x")
  pars <- c(stats::setNames(mu_beta, occ), stats::setNames(sigma_beta, paste0("sd_", occ)),
            `det_(Intercept)` = alpha0, det_bio = a_bio, `sd_det_(Intercept)` = sd_p11,
            `fp_(Intercept)` = alpha0 - mu_gap, fp_gap = mu_gap, sd_log_fp_gap = sd_gap,
            `verified_(Intercept)` = delta0,
            if (use_re) c(sd_recorder = sd_recorder),
            stats::setNames(as.vector(beta), sprintf("%s[%s]", occ, rep(sp, each = 2))),
            stats::setNames(alpha, sprintf("det_(Intercept)[%s]", sp)),
            stats::setNames(gap, sprintf("fp_gap[%s]", sp)),
            # as in fp_occu_ms(): under observation confirmation, p11 includes
            # the false observations at occupied sites
            stats::setNames(if (isTRUE(sl$stan_data$obs_confirm == 1))
                              1 - (1 - stats::plogis(alpha)) * (1 - stats::plogis(alpha - gap))
                            else stats::plogis(alpha), sprintf("p11_ref[%s]", sp)),
            stats::setNames(stats::plogis(alpha - gap), sprintf("p10_ref[%s]", sp)))
  structure(list(data = sims$data,
                 truth = list(pars = pars, z = sims$truth$z, p10_bar = sims$truth$p10_bar)),
            class = "fp_sim_ms")
}

#' @export
print.fp_sim_ms <- function(x, ...) {
  z <- x$truth$z
  cat(sprintf("<fp_sim_ms> %d species at %d sites\n", ncol(z), nrow(z)))
  print(data.frame(species = colnames(z), occupancy = colMeans(z),
                   p10 = unname(x$truth$p10_bar)),
        digits = 3, row.names = FALSE)
  print(x$data)
  invisible(x)
}
