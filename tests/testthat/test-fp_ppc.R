ppc_setup <- function(fp = ~ 1, S = 50, beta1 = 0, gamma0 = -3) {
  s <- simulate_fp(nsite = 40, nocc = 5, frac_verified = 0.5, seed = 3)
  sd <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, fp, FALSE,
                                   fp_priors())$stan_data
  pars <- list(beta = cbind(rep(beta1, S), rep(-0.8, S)),
               alpha = matrix(0.4, S, 1), delta0 = rep(0.5, S), u = NULL,
               gamma0 = if (is.null(fp)) NULL else rep(gamma0, S),
               gamma_fp = matrix(0, S, 0))
  list(sd = sd, pars = pars, S = S)
}

test_that(".ppc_core returns objects of the right shape", {
  set.seed(1)
  k <- ppc_setup()
  r <- crywolf:::.ppc_core(k$sd, k$pars, "automated", "freeman-tukey")
  expect_length(r$y, k$sd$nsite)
  expect_equal(dim(r$yrep), c(k$S, k$sd$nsite))
  expect_length(r$fit_obs, k$S)
  expect_length(r$fit_rep, k$S)
  expect_true(all(r$yrep <= 5))
  expect_equal(sum(r$y), sum(k$sd$y2))

  v <- crywolf:::.ppc_core(k$sd, k$pars, "verified", "chi-squared")
  expect_equal(ncol(v$yrep), sum(tabulate(k$sd$site1, k$sd$nsite) > 0))
  expect_equal(sum(v$y), sum(k$sd$y1))
})

test_that("z is fixed where the data determine it", {
  set.seed(2)
  # sites with a verified detection are occupied
  k <- ppc_setup()
  r <- crywolf:::.ppc_core(k$sd, k$pars, "automated", "freeman-tukey")
  expect_true(all(r$z[, k$sd$known_present == 1] == 1))

  # without false positives, any automated detection means occupied
  k0 <- ppc_setup(fp = NULL)
  r0 <- crywolf:::.ppc_core(k0$sd, k0$pars, "automated", "freeman-tukey")
  det <- tabulate(k0$sd$site2[k0$sd$y2 > 0], nbins = k0$sd$nsite) > 0
  expect_true(all(r0$z[, det] == 1))
})

test_that("the conditional probability of z matches a hand calculation", {
  # one site, one row of 4 occasions with 1 detection, no verified data
  sd <- list(nsite = 1, X = matrix(1, 1, 1), Z = matrix(1, 1, 1),
             estimate_re = 0, estimate_fp = 1, Kfp = 0,
             Xfp = matrix(0, 1, 0), site2 = 1L, rec2 = 1L, n2 = 4L, y2 = 1L,
             len2 = 1L, N1 = 0, y1 = integer(0), site1 = integer(0),
             known_present = 0L)
  S <- 20000
  pars <- list(beta = matrix(0, S, 1), alpha = matrix(qlogis(0.5), S, 1),
               delta0 = rep(0, S), u = NULL, gamma0 = rep(qlogis(0.1), S),
               gamma_fp = matrix(0, S, 0))
  set.seed(4)
  r <- crywolf:::.ppc_core(sd, pars, "automated", "freeman-tukey")
  f1 <- 0.5 * 0.5^4               # psi * p11 (1 - p11)^3
  f0 <- 0.5 * 0.1 * 0.9^3         # (1 - psi) * p10 (1 - p10)^3
  expect_equal(mean(r$z), f1 / (f1 + f0), tolerance = 0.02)
})

test_that("fp_ppc works on a fit", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 100, nocc = 6, seed = 5)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1,
                                random_recorder = FALSE, chains = 2,
                                iter_warmup = 300, iter_sampling = 300,
                                threads_per_chain = 1, seed = 1))
  p <- fp_ppc(m, ndraws = 200, seed = 1)
  expect_s3_class(p, "fp_ppc")
  expect_equal(dim(p$yrep), c(200, 100))
  expect_true(p$p_value >= 0 && p$p_value <= 1)
  expect_output(print(p), "Bayesian p-value")
  pv <- fp_ppc(m, method = "verified", fit_stat = "chi-squared", seed = 1)
  expect_equal(nrow(pv$yrep), 600)
  expect_length(pv$site, ncol(pv$yrep))
})

test_that("the R likelihood matches the Stan log_lik in every design", {
  skip_if_no_cmdstan()
  same_ll <- function(s, ...) {
    m <- suppressWarnings(fp_occu(s$data, occ = ~ x, ..., chains = 1,
                                  iter_warmup = 150, iter_sampling = 100,
                                  threads_per_chain = 1, seed = 3,
                                  show_messages = FALSE, sig_figs = 15))
    pars <- crywolf:::.ppc_pars(m$fit, m$stan_data, NULL)
    r <- crywolf:::.ppc_core(m$stan_data, pars, "automated", "freeman-tukey")
    ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
    expect_lt(max(abs(r$ll - ll) / pmax(1, abs(ll))), 1e-8)
    m
  }
  s1 <- simulate_fp(nsite = 80, nocc = 6, a_bio = 0.5, sd_recorder = 0.5,
                    frac_confirmed = 0.4, frac_verified = 0, seed = 6)
  m1 <- same_ll(s1, det = ~ bio, fp = ~ 1, random_recorder = TRUE)
  expect_output(suppressWarnings(print(m1)), "observation-confirmation design")
  expect_false("p_verified" %in% posterior::variables(m1$draws))
  expect_no_warning(fp_ppc(m1, ndraws = 50, seed = 1))
  s2 <- simulate_fp(nsite = 80, nocc = 6, conf_effect = 1, conf_cor = 0.5, seed = 6)
  same_ll(s2, det = ~ 1, fp = ~ w, random_recorder = FALSE)
  s3 <- simulate_fp(nsite = 80, nocc = 6, sigma_p10 = 1, seed = 7)
  m3 <- same_ll(s3, det = ~ 1, fp = ~ 1, random_recorder = FALSE,
                random_fp = TRUE)
  expect_true("sd_fp_site" %in% posterior::variables(m3$draws))
  expect_output(suppressWarnings(print(m3)), "\\(1 \\| site\\)")
})
