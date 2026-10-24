test_that("Gauss-Hermite quadrature integrates polynomials exactly", {
  gh <- crywolf:::.gauss_hermite(20)
  expect_equal(sum(gh$w), sqrt(pi))
  expect_equal(sum(gh$w * gh$x^2), sqrt(pi) / 2)
  expect_equal(sum(gh$w * gh$x^4), 3 * sqrt(pi) / 4)
  expect_equal(sum(gh$w * gh$x^3), 0)
})

test_that("loo() integrates the site effect on the false-positive rate", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 120, nocc = 8, sigma_p10 = 1, frac_verified = 0.3, seed = 1)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                                random_fp = TRUE, chains = 2, iter_warmup = 300,
                                iter_sampling = 300,
                                threads_per_chain = 1, seed = 1, refresh = 0,
                                show_messages = FALSE))
  set.seed(1)
  before <- .Random.seed
  l1 <- suppressWarnings(loo(m))
  expect_identical(.Random.seed, before)        # no random numbers drawn
  l0 <- suppressWarnings(loo(m, integrate = FALSE))
  expect_lt(sum(l1$diagnostics$pareto_k > 0.7), sum(l0$diagnostics$pareto_k > 0.7))
  expect_true(all(is.finite(l1$pointwise[, "elpd_loo"])))
  # the integral matches a Monte Carlo average over the site effect
  sd <- m$stan_data
  pars <- crywolf:::.ppc_pars(m$fit, sd, NULL)
  q <- crywolf:::.integrate_site_effect(sd, pars, 100)
  r <- 1
  E <- 20000
  p <- purrr::map(pars, function(v) {
    if (is.matrix(v)) v[rep(r, E), , drop = FALSE] else if (length(v)) rep(v[r], E) else v
  })
  set.seed(2)
  p$e_fp <- matrix(pars$sigma_fp[r] * stats::rnorm(E), E, sd$nsite)
  ll <- crywolf:::.ppc_core(sd, p, "automated", "freeman-tukey", ll_only = TRUE)$ll
  mc <- purrr::map_dbl(seq_len(sd$nsite), function(j) crywolf:::.log_mean_exp(ll[, j]))
  expect_lt(stats::median(abs(q[r, ] - mc)), 0.01)
})

test_that("loo() integrates site effects on detection and false positives together", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 80, nocc = 8, sigma_p10 = 1, frac_verified = 0.3, seed = 1)
  warns <- character(0)
  m <- withCallingHandlers(
    fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
            random_fp = TRUE, random_det = TRUE, chains = 2, iter_warmup = 200,
            iter_sampling = 100, threads_per_chain = 1, seed = 1,
            refresh = 0, show_messages = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  expect_match(warns, "mix slowly", all = FALSE)
  sd <- m$stan_data
  pars <- crywolf:::.ppc_pars(m$fit, sd, 5)
  q <- crywolf:::.integrate_site_effect(sd, pars, nodes = 20)
  # the two-dimensional rule against a Monte Carlo average, at one draw
  r <- 1
  E <- 20000
  p <- purrr::map(pars, function(v) {
    if (is.matrix(v)) v[rep(r, E), , drop = FALSE] else if (length(v)) rep(v[r], E) else v
  })
  set.seed(2)
  p$e_fp <- matrix(pars$sigma_fp[r] * stats::rnorm(E * sd$nsite), E)
  p$d_det <- matrix(pars$sigma_det[r] * stats::rnorm(E * sd$nsite), E)
  ll <- crywolf:::.ppc_core(sd, p, "automated", "freeman-tukey", ll_only = TRUE)$ll
  mc <- purrr::map_dbl(seq_len(sd$nsite), function(j) crywolf:::.log_mean_exp(ll[, j]))
  expect_lt(stats::median(abs(q[r, ] - mc)), 0.01)
  expect_true(all(is.finite(suppressWarnings(loo(m, nodes = 10))$pointwise[, "elpd_loo"])))
  expect_s3_class(simulate(m, seed = 1)[[1]]$data, "fp_data")
  expect_s3_class(suppressWarnings(fp_ppc(m, ndraws = 20, seed = 1)), "fp_ppc")
  expect_output(suppressWarnings(print(m)), "true detection: .*\\(1 \\| site\\)")
})

test_that("the integrated log-likelihood keeps its shape when draws are thinned", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, sigma_p10 = 1, seed = 2)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                                random_fp = TRUE, chains = 2, iter_warmup = 200,
                                iter_sampling = 200, thin = 2,
                                threads_per_chain = 1, seed = 1, refresh = 0,
                                show_messages = FALSE))
  l <- suppressWarnings(loo(m))
  # 2 chains x 200 iterations thinned by 2: 200 draws, not 400
  expect_equal(attr(l, "dims"), c(200L, 60L))
  expect_equal(nrow(l$pointwise), 60)
  expect_true(all(is.finite(l$pointwise[, "elpd_loo"])))
})
