test_that("fp_occu fits and its methods work", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 150, nocc = 8, conf_effect = 1.5, conf_cor = 0.7,
                   seed = 11)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ w,
                                random_recorder = FALSE, chains = 2,
                                iter_warmup = 300, iter_sampling = 300,
                                seed = 1))
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
            iter_sampling = 100, seed = 1, show_messages = FALSE)
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
