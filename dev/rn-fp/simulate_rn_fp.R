# Royle-Nichols occupancy with false positives: simulator for step 1
# (dev/roadmap-v2.md). N_i ~ Poisson(lambda_i), log lambda_i = b0 + b1 x_i;
# (plus a site effect N(0, sigma_lam) on log lambda); per-individual detection r_i, logit r_i = a0 + e_i, e_i ~ N(0, sigma_r);
# automated detection p11_i = 1 - (1 - r_i)^N_i at occupied sites (N_i > 0),
# false positives with probability p10 on every occasion; verified occasions
# (by ear, no false positives) at a share of the sites, with per-individual
# detection rv, so s_i = 1 - (1 - rv)^N_i.
simulate_rn_fp <- function(nsite = 200, nocc = 10, b0 = 0, b1 = 0.5,
                           a0 = -1, sigma_r = 0.5, sigma_lam = 0, p10 = 0.05, rv = 0.5,
                           frac_verified = 0.3, nver = 3, seed = 1) {
  set.seed(seed)
  x <- stats::rnorm(nsite)
  lambda <- exp(b0 + b1 * x + stats::rnorm(nsite, 0, sigma_lam))
  N <- stats::rpois(nsite, lambda)
  e <- stats::rnorm(nsite, 0, sigma_r)
  r <- stats::plogis(a0 + e)
  p11 <- 1 - (1 - r)^N
  p_auto <- 1 - (1 - p11) * (1 - p10)          # N = 0 gives p10
  k <- stats::rbinom(nsite, nocc, p_auto)       # automated detections
  is_ver <- stats::rbinom(nsite, 1, frac_verified) == 1
  nv <- ifelse(is_ver, nver, 0L)
  s <- 1 - (1 - rv)^N
  kv <- stats::rbinom(nsite, nv, s)
  list(data = list(nsite = nsite, nocc = nocc, x = x, k = k, nv = as.integer(nv),
                   kv = kv),
       truth = list(b0 = b0, b1 = b1, a0 = a0, sigma_r = sigma_r, sigma_lam = sigma_lam, p10 = p10,
                    rv = rv, N = N, z = as.integer(N > 0), r = r, p11 = p11,
                    psi = mean(1 - exp(-lambda))))
}
