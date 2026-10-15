test_that("simulate_fp returns well-formed data", {
  s <- simulate_fp(nsite = 60, nocc = 4, nver = 2, seed = 3)
  expect_s3_class(s, "fp_sim")
  expect_s3_class(s$data, "fp_data")
  expect_equal(nrow(s$data$sites), 60)
  expect_equal(nrow(s$data$auto), 240)
  expect_equal(nrow(s$data$verified), 2 * s$truth$n_verified_sites)
  expect_equal(s$truth$regime, "homogeneous")
  expect_equal(s$truth$p10_marginal, 0.05)
})

test_that("simulate_fp is reproducible with a seed", {
  a <- simulate_fp(nsite = 30, nocc = 3, seed = 10)
  b <- simulate_fp(nsite = 30, nocc = 3, seed = 10)
  expect_identical(a$data$auto, b$data$auto)
})

test_that("verified detections never occur at unoccupied sites", {
  s <- simulate_fp(nsite = 200, nocc = 2, frac_verified = 1, seed = 4)
  v <- s$data$verified
  expect_true(all(v$y[s$truth$z[v$site] == 0] == 0))
})

test_that("structured confusion sets the gradient and raises the mean rate", {
  s <- simulate_fp(nsite = 2000, nocc = 1, conf_effect = 1.5, conf_cor = 0.7,
                   seed = 5)
  expect_equal(s$truth$regime, "structured")
  expect_equal(s$truth$conf_cor_realized, 0.7, tolerance = 0.05)
  expect_gt(s$truth$p10_marginal, 0.05)
})

test_that("unstructured heterogeneity is flagged", {
  expect_equal(simulate_fp(nsite = 20, nocc = 2, sigma_p10 = 1, seed = 1)$truth$regime,
               "unstructured")
  expect_equal(simulate_fp(nsite = 20, nocc = 2, p10_conc = 2, seed = 1)$truth$regime,
               "unstructured")
})

test_that("frac_confirmed adds confirmations without changing the default draws", {
  a <- simulate_fp(nsite = 60, nocc = 5, seed = 9)
  expect_identical(simulate_fp(nsite = 60, nocc = 5, frac_confirmed = 0, seed = 9), a)
  b <- simulate_fp(nsite = 60, nocc = 5, frac_confirmed = 0.4, seed = 9)
  unocc <- rep(a$truth$z, 5) == 0
  expect_identical(b$data$auto$y[unocc], a$data$auto$y[unocc])
  expect_identical(b$data$verified, a$data$verified)
  conf <- !is.na(b$data$auto$n_true)
  expect_true(any(conf))
  expect_true(all(b$data$auto$n_true[conf & unocc] == 0))
  expect_equal(b$data$confirmed, c("n_true", "n_false"))
})
