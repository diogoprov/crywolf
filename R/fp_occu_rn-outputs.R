# Outputs of fp_occu_rn(): leave-one-site-out, site estimates, posterior
# predictive check and simulation. They rest on .rn_core(), which computes,
# for each posterior draw and site, the joint log probability of the site's
# data and each abundance 0..Nmax, as the Stan program does.

#' For each draw (rows) and site (columns): the marginal log-likelihood `ll`
#' (the log_lik of fp_occu_rn.stan), the probability that the site is
#' occupied given its data `pz`, the posterior mean abundance given its data
#' `n_mean`, and, with draw = TRUE, one abundance `N` drawn from its
#' conditional posterior (by the Gumbel-max trick, so that the values of N
#' are visited once and nothing of size Nmax is kept)
#' @noRd
.rn_core <- function(sd, pars, draw = FALSE) {
  S <- nrow(pars$beta)
  nsite <- sd$nsite
  rows <- function(v) rep(v, each = S)
  # log(1 - p) from the logit of p, and log(1 - exp(x)) for x <= 0
  lmiss <- function(eta) stats::plogis(-eta, log.p = TRUE)
  log1m_exp <- function(x) log(-expm1(x))

  llam <- pars$beta %*% t(sd$X)
  if (!is.null(pars$v_ab)) llam <- llam + pars$v_ab
  fp <- sd$estimate_fp == 1
  lp10 <- if (fp) .lp10_draws(sd, pars) else matrix(-Inf, S, nsite)
  lq0 <- stats::plogis(-lp10, log.p = TRUE)
  ls0 <- stats::plogis(lp10, log.p = TRUE)

  # log(1 - r) of each automated row, and log(1 - p10) of its site
  eta2 <- .rn_eta(pars, sd$Z, sd$rec2, sd$site2, sd$estimate_re == 1)
  lr2 <- lmiss(eta2)
  l10 <- if (fp) lq0[, sd$site2, drop = FALSE] else 0
  y2 <- sd$y2
  f2 <- sd$n2 - sd$y2
  det2 <- which(y2 > 0)
  if (sd$N1 > 0) {
    lr1 <- lmiss(.eta_verified(sd, pars))
    det1 <- which(sd$y1 == 1)
  }
  conf <- isTRUE(sd$obs_confirm == 1) && sd$N3 > 0
  if (conf) {
    lr3 <- lmiss(.rn_eta(pars, sd$Z3, sd$rec3, sd$site3, sd$estimate_re == 1))
    v3 <- sd$v3
    w3 <- rows(sd$w3)
    # log P(false observations | site), the same for every N
    odd <- v3 %% 2 == 1
    fpart3 <- lq0[, sd$site3, drop = FALSE]
    fpart3[, odd] <- ls0[, sd$site3[odd], drop = FALSE]
  }

  m <- matrix(-Inf, S, nsite)   # running maximum of the joint log probability
  s <- matrix(0, S, nsite)      # sum of exp(joint - m)
  t <- matrix(0, S, nsite)      # sum of n exp(joint - m)
  best <- matrix(-Inf, S, nsite)
  N <- if (draw) matrix(0L, S, nsite)
  joint0 <- NULL
  for (n in 0:sd$Nmax) {
    lm2 <- n * lr2 + l10
    t2 <- rows(f2) * lm2
    t2[, det2] <- t2[, det2] + rows(y2[det2]) * log1m_exp(lm2[, det2, drop = FALSE])
    em <- .site_sum(t2, sd$site2, nsite)
    if (sd$N1 > 0) {
      lm1 <- n * lr1
      t1 <- lm1
      t1[, det1] <- log1m_exp(lm1[, det1, drop = FALSE])
      em <- em + .site_sum(t1, sd$site1, nsite)
    }
    if (conf) {
      lq1 <- n * lr3
      tpart <- lq1
      tpart[, v3 >= 2] <- log1m_exp(lq1[, v3 >= 2, drop = FALSE])
      em <- em + .site_sum(w3 * (tpart + fpart3), sd$site3, nsite)
    }
    joint <- n * llam - exp(llam) - lgamma(n + 1) + em
    if (n == 0) {
      joint[, sd$known_present == 1] <- -Inf
      joint0 <- joint
    }
    m_new <- pmax(m, joint)
    scale <- ifelse(is.finite(m), exp(m - m_new), 0)
    add <- ifelse(is.finite(joint), exp(joint - m_new), 0)
    s <- s * scale + add
    t <- t * scale + n * add
    m <- m_new
    if (draw) {
      g <- joint - log(-log(stats::runif(S * nsite)))
      take <- g > best
      best[take] <- g[take]
      N[take] <- n
    }
  }
  ll <- m + log(s)
  list(ll = ll, pz = -expm1(joint0 - ll), n_mean = t / s, N = N)
}

#' logit of the per-individual detection of each row of a design matrix
#' (automated or confirmed rows), S draws x rows
#' @noRd
.rn_eta <- function(pars, Z, rec, site, has_re) {
  eta <- pars$alpha %*% t(Z)
  if (has_re) eta <- eta + pars$u[, rec, drop = FALSE]
  if (!is.null(pars$d_det)) eta <- eta + pars$d_det[, site, drop = FALSE]
  eta
}

#' @rdname fp_occu_rn-methods
#' @inheritParams loo.fp_occu
#' @param integrate Logical. With site random effects (`random_det`,
#'   `random_abund` or `random_fp`), integrate each site's effects out of its
#'   log-likelihood by Gauss-Hermite quadrature before the leave-one-out
#'   computation, as [loo.fp_occu()] does. The default number of nodes per
#'   effect is 100 with one effect, 30 with two and 12 with three; each node
#'   of the grid costs one pass over the abundances 0..`n_max`, so the time
#'   grows with `nodes`, `n_max` and the number of draws. In a test with 40
#'   sites, the log-likelihood of a site moved by at most 0.005 from 30 to
#'   60 nodes with two effects, and by at most 0.07 (median 1e-7) from 12
#'   to 20 nodes with three, where 20 nodes took about five times as long.
#'   Raise `nodes` when models are close; `integrate = FALSE` is quicker but
#'   gives many high Pareto \eqn{k}.
#' @export
loo.fp_occu_rn <- function(x, reloo = FALSE, k_threshold = 0.7, seed = NULL,
                           integrate = TRUE, nodes = NULL, ...) {
  loo.fp_occu(x, reloo = reloo, k_threshold = k_threshold, seed = seed,
              integrate = integrate, nodes = nodes, ...)
}

#' @rdname fp_occu_rn-methods
#' @inheritParams simulate.fp_occu
#' @export
simulate.fp_occu_rn <- function(object, nsim = 1, seed = NULL,
                                truth = c("mean", "posterior"), ...) {
  simulate.fp_occu(object, nsim = nsim, seed = seed, truth = truth, ...)
}

#' fp_sites() of a Royle-Nichols fit
#' @noRd
.fp_sites_rn <- function(object, prob, ndraws, seed) {
  sd <- object$stan_data
  if (!is.null(seed)) .set_seed(seed)
  pars <- .ppc_pars(object$fit, sd, ndraws)
  r <- .rn_core(sd, pars)
  a <- (1 - prob) / 2
  out <- .sites_table(object$data$sites, sd, pars, r$pz, prob,
                      psi = -expm1(-exp(pars$beta %*% t(sd$X))),
                      extra = c("abundance", "abundance_lower", "abundance_upper"))
  out$abundance <- colMeans(r$n_mean)
  out$abundance_lower <- .col_quantile(r$n_mean, a)
  out$abundance_upper <- .col_quantile(r$n_mean, 1 - a)
  structure(out, class = c("fp_sites", "data.frame"), prob = prob)
}

#' fp_ppc() of a Royle-Nichols fit: abundances drawn from their conditional
#' posterior, then the detections replicated given them
#' @noRd
.fp_ppc_rn <- function(object, method, fit_stat, ndraws, seed) {
  sd <- object$stan_data
  if (!is.null(seed)) .set_seed(seed)
  pars <- .ppc_pars(object$fit, sd, ndraws)
  N <- .rn_core(sd, pars, draw = TRUE)$N
  S <- nrow(N)
  nsite <- sd$nsite
  rows <- function(v) rep(v, each = S)
  # probability of at least one detection on an occasion, given N
  p_det <- function(eta, site, lq0) {
    -expm1(N[, site, drop = FALSE] * stats::plogis(-eta, log.p = TRUE) + lq0)
  }
  if (method == "automated") {
    lq0 <- if (sd$estimate_fp == 1)
      stats::plogis(-.lp10_draws(sd, pars), log.p = TRUE) else matrix(0, S, nsite)
    p2 <- p_det(.rn_eta(pars, sd$Z, sd$rec2, sd$site2, sd$estimate_re == 1),
                sd$site2, lq0[, sd$site2, drop = FALSE])
    n_row <- rows(sd$n2)
    y <- tabulate(rep(sd$site2, sd$y2), nbins = nsite)
    E <- .site_sum(n_row * p2, sd$site2, nsite)
    yrep <- .site_sum(matrix(stats::rbinom(length(p2), n_row, p2), S), sd$site2, nsite)
    keep <- sd$len2 > 0
    if (isTRUE(sd$obs_confirm == 1) && sd$N3 > 0) {
      # confirmed occasions enter as detections or not, as in fp_ppc() of
      # fp_occu(); each row holds w3 identical occasions
      p3 <- p_det(.rn_eta(pars, sd$Z3, sd$rec3, sd$site3, sd$estimate_re == 1),
                  sd$site3, lq0[, sd$site3, drop = FALSE])
      w <- rows(sd$w3)
      keep <- keep | tabulate(sd$site3, nbins = nsite) > 0
      y <- y + tabulate(rep(sd$site3, sd$w3 * (sd$v3 > 0)), nbins = nsite)
      E <- E + .site_sum(w * p3, sd$site3, nsite)
      yrep <- yrep + .site_sum(matrix(stats::rbinom(length(p3), w, p3), S),
                               sd$site3, nsite)
    }
  } else {
    p1 <- p_det(.eta_verified(sd, pars), sd$site1, 0)
    keep <- tabulate(sd$site1, nbins = nsite) > 0
    y <- tabulate(sd$site1[sd$y1 == 1], nbins = nsite)
    E <- .site_sum(p1, sd$site1, nsite)
    yrep <- .site_sum(matrix(stats::rbinom(length(p1), 1, p1), S), sd$site1, nsite)
  }
  y <- y[keep]
  E <- E[, keep, drop = FALSE]
  yrep <- yrep[, keep, drop = FALSE]
  Y <- matrix(y, S, length(y), byrow = TRUE)
  disc <- switch(fit_stat,
    "freeman-tukey" = function(obs) rowSums((sqrt(obs) - sqrt(E))^2),
    "chi-squared"   = function(obs) rowSums((obs - E)^2 / (E + 1e-6))
  )
  fit_obs <- disc(Y)
  fit_rep <- disc(yrep)
  structure(
    list(y = y, yrep = yrep, fit_obs = fit_obs, fit_rep = fit_rep,
         p_value = mean(fit_rep > fit_obs),
         site = object$data$sites[[object$data$site]][keep],
         method = method, fit_stat = fit_stat),
    class = "fp_ppc"
  )
}
