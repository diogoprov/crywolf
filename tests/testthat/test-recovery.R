fit_quick <- function(data, ..., seed = 1) {
  suppressWarnings(fp_occu(data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                           chains = 2, iter_warmup = 500, iter_sampling = 500,
                           threads_per_chain = 1, seed = seed, ...))
}
inside <- function(m, par, truth) {
  q <- stats::quantile(posterior::extract_variable(m$draws, par), c(0.025, 0.975))
  truth >= q[[1]] && truth <= q[[2]]
}

#' @srrstats {BS7.0, BS7.1} The prior is recovered by sampling without the likelihood (`prior_only = TRUE`).
test_that("sampling from the priors alone recovers the priors", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 4)
  pr <- fp_priors()
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                                chains = 4, iter_warmup = 500, iter_sampling = 1000,
                                threads_per_chain = 1, seed = 2, prior_only = TRUE))
  expect_true(m$prior_only)
  expect_output(print(m), "priors alone")
  expect_error(loo(m), "priors alone")
  dr <- posterior::as_draws_df(m$draws)
  expect_equal(sd(dr[["occ_(Intercept)"]]), pr$intercept_sd, tolerance = 0.1)
  expect_equal(sd(dr[["occ_x"]]), pr$slope_sd, tolerance = 0.1)
  expect_equal(mean(dr[["occ_x"]]), 0, tolerance = 0.25)
  # Normal(gap_mean, gap_sd) truncated at zero
  a <- -pr$gap_mean / pr$gap_sd
  gap_mean <- pr$gap_mean + pr$gap_sd * dnorm(a) / (1 - pnorm(a))
  expect_equal(mean(dr[["fp_gap"]]), gap_mean, tolerance = 0.1)
})

#' @srrstats {G5.6, G5.6a, BS7.2} Known parameters fall inside their 95% intervals.
test_that("fp_occu() recovers known parameters from simulated data", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 300, nocc = 20, seed = 21)
  m <- fit_quick(s$data)
  t <- s$truth
  expect_true(inside(m, "occ_x", t$b1))
  expect_true(inside(m, "p11_ref", t$p11))
  expect_true(inside(m, "p10_ref", t$p10))
})

#' @srrstats {BS5.0} The stored seed reproduces the fit.
test_that("the seed CmdStan used is stored and reproduces the fit", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 5)
  m1 <- suppressWarnings(fp_occu(s$data, occ = ~ x, random_recorder = FALSE, chains = 1,
                                 iter_warmup = 200, iter_sampling = 200,
                                 threads_per_chain = 1, seed = NULL))
  expect_true(is.numeric(m1$sampler$seed) && length(m1$sampler$seed) >= 1)
  m2 <- suppressWarnings(fp_occu(s$data, occ = ~ x, random_recorder = FALSE, chains = 1,
                                 iter_warmup = 200, iter_sampling = 200,
                                 threads_per_chain = 1, seed = m1$sampler$seed[1]))
  expect_identical(coef(m1), coef(m2))
})

#' @srrstats {G5.6b} Recovery over ten seeds (extended).
test_that("recovery holds over several datasets (extended)", {
  skip_if_not_extended()
  hits <- vapply(1:10, function(i) {
    s <- simulate_fp(nsite = 300, nocc = 20, seed = 100 + i)
    inside(fit_quick(s$data, seed = i), "occ_x", s$truth$b1)
  }, logical(1))
  # with 95% intervals, 10 datasets give at least 8 hits with probability 0.99
  expect_gte(sum(hits), 8)
})

#' @srrstats {G5.9, G5.9a, G5.9b} Noise of machine-epsilon scale and another seed change estimates by less than five MCSE (extended).
test_that("estimates are stable to noise in the data and to the seed (extended)", {
  skip_if_not_extended()
  s <- simulate_fp(nsite = 300, nocc = 20, seed = 31)
  m <- fit_quick(s$data, seed = 1)
  noisy <- s$data
  noisy$sites$x <- noisy$sites$x + stats::runif(nrow(noisy$sites), -1, 1) *
    sqrt(.Machine$double.eps)
  vars <- c("occ_(Intercept)", "occ_x", "p11_ref", "p10_ref")
  summ <- function(f) posterior::summarise_draws(posterior::subset_draws(f$draws, vars),
                                                 "mean", "mcse_mean")
  ref <- summ(m)
  for (other in list(fit_quick(noisy, seed = 1), fit_quick(s$data, seed = 2))) {
    o <- summ(other)
    # two independent runs differ by about sqrt(2) MCSE; allow five
    expect_true(all(abs(o$mean - ref$mean) < 5 * sqrt(o$mcse_mean^2 + ref$mcse_mean^2)))
  }
})
