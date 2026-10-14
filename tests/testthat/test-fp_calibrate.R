test_that("fp_calibrate rejects unknown design columns", {
  expect_error(fp_calibrate(data.frame(nope = 1), reps = 1), "Unknown")
})

test_that("fp_calibrate runs a tiny design", {
  skip_if_no_cmdstan()
  cal <- fp_calibrate(data.frame(nsite = 80, nocc = 5),
                      models = list(constant = ~ 1, fp_x = ~ x),
                      reps = 2,
                      fit_args = list(chains = 1, iter_warmup = 150,
                                      iter_sampling = 150))
  expect_s3_class(cal, "fp_calibration")
  s <- summary(cal)
  expect_true(all(c("coverage", "mc_se", "bias", "rmse") %in% names(s)))
  expect_setequal(unique(s$model), c("constant", "fp_x"))
})
