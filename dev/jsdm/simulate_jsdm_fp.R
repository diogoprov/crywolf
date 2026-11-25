# Simulator for the JSDM prototype (dev/jsdm/jsdm_fp.stan): occupancy of J
# species at S sites with K latent factors, automated detections with false
# positives, verified occasions at a share of the sites, and an optional site
# effect on detection shared by all species (correlated detection).
simulate_jsdm_fp <- function(S = 100, J = 10, nocc = 10, K = 2, load_sd = 1,
                             b0_mu = 0, b0_sd = 1, b1_mu = 0.5, b1_sd = 0.5,
                             p11 = c(0.2, 0.6), p10 = c(0.01, 0.05),
                             pv = 0.5, frac_verified = 0.3, nver = 3,
                             sigma_d = 0, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x <- as.numeric(scale(stats::rnorm(S)))
  b0 <- stats::rnorm(J, b0_mu, b0_sd)
  b1 <- stats::rnorm(J, b1_mu, b1_sd)
  Lambda <- matrix(0, J, K)
  if (K > 0) {
    Lambda[] <- stats::rnorm(J * K, 0, load_sd)
    Lambda[upper.tri(Lambda)] <- 0
    for (k in seq_len(K)) Lambda[k, k] <- abs(Lambda[k, k])
  }
  Z <- matrix(stats::rnorm(K * S), K, S)
  eta <- outer(x, b1) + matrix(b0, S, J, byrow = TRUE) +
    if (K > 0) t(Lambda %*% Z) else 0
  psi <- stats::plogis(eta)
  z <- matrix(stats::rbinom(S * J, 1, psi), S, J)
  a <- stats::qlogis(stats::runif(J, p11[1], p11[2]))
  q10 <- stats::runif(J, p10[1], p10[2])
  d <- stats::rnorm(S, 0, sigma_d)
  p_occ <- stats::plogis(outer(d, a, "+"))
  p <- ifelse(z == 1, p_occ, matrix(q10, S, J, byrow = TRUE))
  n <- matrix(nocc, S, J)
  y <- matrix(stats::rbinom(S * J, n, p), S, J)
  ver <- stats::rbinom(S, 1, frac_verified) == 1
  nv <- matrix(ifelse(ver, nver, 0L), S, J)
  v <- matrix(stats::rbinom(S * J, nv, pv * z), S, J)
  list(data = list(S = S, J = J, x = x, n = n, y = y, nv = nv, v = v),
       truth = list(b0 = b0, b1 = b1, a = a, p10 = q10, Lambda = Lambda,
                    LL = Lambda %*% t(Lambda), Z = Z, z = z, sigma_d = sigma_d,
                    occupancy = colMeans(z)))
}
