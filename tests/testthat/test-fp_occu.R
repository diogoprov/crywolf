test_that("fp_occu fits and its methods work", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 150, nocc = 8, conf_effect = 1.5, conf_cor = 0.7,
                   seed = 11)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ w,
                                random_recorder = FALSE, chains = 2,
                                iter_warmup = 300, iter_sampling = 300,
                                threads_per_chain = 1, seed = 1))
  expect_s3_class(m, "fp_occu")
  cf <- coef(m)
  expect_true(all(c("occ_(Intercept)", "occ_x", "fp_w", "p10_bar") %in% names(cf)))
  expect_named(coef(m, "occ"), c("occ_(Intercept)", "occ_x"))
  expect_s3_class(summary(m)$table, "tbl_df")
  p <- predict(m, newdata = data.frame(x = c(-1, 0, 1)))
  expect_equal(nrow(p), 3)
  expect_true(all(p$estimate > 0 & p$estimate < 1))
  expect_s3_class(loo(m), "psis_loo")
  expect_output(print(m), "false positive")
})

test_that("fp_occu warns without verified data", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 50, nocc = 4, frac_verified = 0, seed = 2)
  w <- testthat::capture_warnings(
    fp_occu(s$data, occ = ~ x, chains = 1, iter_warmup = 100,
            iter_sampling = 100, threads_per_chain = 1, seed = 1, show_messages = FALSE)
  )
  expect_match(w, "No verified data", all = FALSE)
})

test_that("threads_per_chain is checked and does not change the fit", {
  s <- simulate_fp(nsite = 80, nocc = 6, seed = 5)
  expect_error(fp_occu(s$data, occ = ~ x, random_recorder = FALSE,
                       threads_per_chain = 0), "positive integer")
  skip_if_no_cmdstan()
  expect_true(crywolf:::.has_threads(fp_compile()))
  fit <- function(threads) {
    suppressWarnings(fp_occu(s$data, occ = ~ x, random_recorder = FALSE,
                             chains = 1, iter_warmup = 200,
                             iter_sampling = 200, seed = 3,
                             threads_per_chain = threads, refresh = 0))
  }
  m1 <- fit(1)
  m2 <- fit(2)
  # same seed and the same log density up to rounding: the first draws agree
  d1 <- posterior::as_draws_matrix(m1$draws)[1:20, ]
  d2 <- posterior::as_draws_matrix(m2$draws)[1:20, ]
  expect_equal(unclass(d1), unclass(d2), tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("threads default to \"auto\", but one per chain in repeated fits", {
  expect_identical(formals(fp_occu)$threads_per_chain, "auto")
  # the fits of a calibration or SBC run side by side: one thread each,
  # unless fit_args says otherwise
  got <- list()
  local_mocked_bindings(fp_occu = function(...) {
    got[[length(got) + 1]] <<- list(...)$threads_per_chain
    stop("not fitted")
  }, .package = "crywolf")
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  try(fp_sbc(s$data, reps = 1, ndraws = 9), silent = TRUE)
  try(fp_sbc(s$data, reps = 1, ndraws = 9,
             fit_args = list(chains = 4, iter_sampling = 100, threads_per_chain = 2)),
      silent = TRUE)
  try(fp_calibrate(data.frame(nsite = 20, nocc = 3), reps = 1), silent = TRUE)
  expect_equal(unlist(got), c(1, 2, 1))
})

test_that("R CMD check limits a fit to two cores", {
  old <- Sys.getenv("_R_CHECK_LIMIT_CORES_", NA)
  on.exit(if (is.na(old)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_") else
    Sys.setenv(`_R_CHECK_LIMIT_CORES_` = old))
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "FALSE")
  expect_false(crywolf:::.check_limits_cores())
  Sys.setenv(`_R_CHECK_LIMIT_CORES_` = "TRUE")
  expect_true(crywolf:::.check_limits_cores())
  expect_lte(crywolf:::.auto_threads(1), 2L)
  expect_equal(crywolf:::.auto_threads(4), 1L)
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  m <- suppressWarnings(fp_occu(s$data, random_recorder = FALSE, chains = 2,
                                threads_per_chain = 4, iter_warmup = 50,
                                iter_sampling = 50, seed = 1, refresh = 0,
                                show_messages = FALSE))
  expect_equal(m$sampler$threads_per_chain, 2L)
  expect_equal(m$sampler$parallel_chains, 1L)
})

test_that("threads_per_chain resolves \"auto\" and rejects other text", {
  expect_gte(crywolf:::.auto_threads(4), 1L)
  expect_equal(crywolf:::.auto_threads(1e6), 1L)
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  expect_error(suppressWarnings(fp_occu(s$data, threads_per_chain = "many")),
               "positive integer or \"auto\"")
  expect_error(suppressWarnings(fp_occu(s$data, threads_per_chain = 0)),
               "positive integer")
})

test_that("fits after the same seed get their own output files", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 30, nocc = 3, seed = 1)
  fit <- function() {
    set.seed(1)
    suppressWarnings(fp_occu(s$data, occ = ~ x, random_recorder = FALSE, chains = 1,
                             iter_warmup = 50, iter_sampling = 50,
                             threads_per_chain = 1, seed = 1,
                             refresh = 0, show_messages = FALSE))
  }
  f1 <- fit()
  f2 <- fit()
  expect_length(intersect(f1$fit$output_files(), f2$fit$output_files()), 0)
  b <- crywolf:::.output_basename()
  expect_false(identical(b, crywolf:::.output_basename()))
})

test_that("CmdStan versions below 2.40 are not usable", {
  expect_false(crywolf:::.cmdstan_recent("2.39.0"))
  expect_true(crywolf:::.cmdstan_recent("2.40.0"))
  expect_true(crywolf:::.cmdstan_recent("2.41.2"))
})
