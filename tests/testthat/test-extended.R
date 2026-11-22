# Long statistical tests, run only with CRYWOLF_EXTENDED_TESTS=true
# (see tests/README.md).

fit_ext <- function(data, occ = ~ x, ..., chains = 4, seed = 1) {
  suppressWarnings(fp_occu(data, occ = occ, fp = ~ 1, random_recorder = FALSE,
                           chains = chains, iter_warmup = 500, iter_sampling = 500,
                           threads_per_chain = 1, seed = seed, ...))
}
post <- function(m, par) posterior::extract_variable(m$draws, par)
mcse <- function(m, par) posterior::mcse_mean(post(m, par))

test_that("estimates agree with unmarked::occuFP() (extended)", {
  skip_if_not_extended()
  skip_if_not_installed("unmarked")
  s <- simulate_fp(nsite = 300, nocc = 10, frac_verified = 0.3, nver = 3, seed = 2026)
  umf <- fp_to_unmarked(s$data)
  m_um <- unmarked::occuFP(detformula = ~ METH, FPformula = ~ 1,
                           stateformula = ~ x, data = umf)
  cu <- unmarked::coef(m_um)
  se <- sqrt(diag(unmarked::vcov(m_um)))
  m_cw <- fit_ext(s$data)
  cw <- coef(m_cw)
  # posterior means within one maximum-likelihood standard error
  expect_lt(abs(cw[["occ_x"]] - cu[["psi(x)"]]), se[["psi(x)"]])
  expect_lt(abs(cw[["occ_(Intercept)"]] - cu[["psi(Int)"]]), se[["psi(Int)"]])
  expect_lt(abs(cw[["fp_(Intercept)"]] - cu[["fp(Int)"]]), se[["fp(Int)"]])
})

test_that("precision improves with the number of sites (extended)", {
  skip_if_not_extended()
  sd_x <- vapply(c(100, 400), function(n) {
    s <- simulate_fp(nsite = n, nocc = 10, seed = 41)
    stats::sd(post(fit_ext(s$data), "occ_x"))
  }, numeric(1))
  # four times the sites should about halve the posterior sd
  expect_lt(sd_x[2], 0.75 * sd_x[1])
})

test_that("run time grows about linearly with the number of sites (extended)", {
  skip_if_not_extended()
  fp_compile()
  secs <- vapply(c(150, 600), function(n) {
    s <- simulate_fp(nsite = n, nocc = 10, seed = 42)
    system.time(fit_ext(s$data, chains = 1, parallel_chains = 1))[["elapsed"]]
  }, numeric(1))
  # four times the sites: allow up to ten times the time (start-up included)
  expect_lt(secs[2] / secs[1], 10)
})

test_that("the scale of a covariate changes only the scale of its effect (extended)", {
  skip_if_not_extended()
  s <- simulate_fp(nsite = 300, nocc = 10, seed = 43)
  m1 <- fit_ext(s$data)
  d2 <- s$data
  d2$sites$x10 <- d2$sites$x * 10 + 5
  m2 <- fit_ext(d2, occ = ~ x10)
  tol <- 5 * sqrt(mcse(m1, "occ_x")^2 + (10 * mcse(m2, "occ_x10"))^2)
  expect_lt(abs(mean(post(m1, "occ_x")) - 10 * mean(post(m2, "occ_x10"))), tol)
  # the probabilities do not depend on the scale of the occupancy covariate
  expect_lt(abs(mean(post(m1, "p10_ref")) - mean(post(m2, "p10_ref"))),
            5 * sqrt(mcse(m1, "p10_ref")^2 + mcse(m2, "p10_ref")^2))
})

test_that("replicated detections match the observed scale (extended)", {
  skip_if_not_extended()
  s <- simulate_fp(nsite = 300, nocc = 10, seed = 44)
  pp <- fp_ppc(fit_ext(s$data), ndraws = 500, seed = 1)
  expect_equal(mean(pp$yrep), mean(pp$y), tolerance = 0.05)
})
