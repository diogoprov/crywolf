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
