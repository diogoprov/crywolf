# Conditions of fp_data(), fp_data_wide() and the data checks they share
# that are not exercised elsewhere (test-fp_data.R, test-fp_occu_dyn.R,
# test-input-validation.R and test-fp_occu_ms.R cover the others).

#' @srrstats {G5.2b} Conditions raised by fp_data(), fp_data_wide() and the simulators.
test_that("fp_data() stops on non-tabular input and unusable detections", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 1, 2, 2), y = c(1, 0, 0, 1))
  expect_error(fp_data(as.matrix(auto), sites), "`auto` must be a data frame.",
               fixed = TRUE)
  expect_error(fp_data(auto, list(site = 1:2)), "`sites` must be a data frame.",
               fixed = TRUE)
  expect_error(fp_data(auto, sites, verified = "not a table"),
               "`verified` must be a data frame.", fixed = TRUE)
  # every detection missing: no occasion is left once they are dropped
  expect_error(fp_data(transform(auto, y = NA), sites),
               "`auto` has no automated occasion with a 0/1 detection in `y`.",
               fixed = TRUE)
  expect_error(fp_data(transform(auto, det = NA_real_), sites, y = "det"),
               "with a 0/1 detection in `det`.", fixed = TRUE)
})

test_that("fp_data() stops on missing recorders", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 1, 2, 2), y = c(1, 0, 0, 1),
                     rec = c("a", NA, "b", "b"))
  expect_error(fp_data(auto, sites, recorder = "rec"),
               "Column `rec` of `auto` has missing values.", fixed = TRUE)
  expect_s3_class(fp_data(transform(auto, rec = "a"), sites, recorder = "rec"),
                  "fp_data")
})

test_that("fp_data() stops on verified seasons outside the factor levels", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 1, 2, 2), y = c(1, 0, 0, 1),
                     t = factor(c("s1", "s2", "s1", "s2")))
  ver <- data.frame(site = 1:2, y = c(1, 0), t = c("s1", "s9"))
  expect_error(fp_data(auto, sites, verified = ver, season = "t"),
               "Seasons in `verified` are not levels of `auto$t` (e.g. s9).",
               fixed = TRUE)
  # missing seasons in `verified` stop too, not only those in `auto`
  expect_error(fp_data(auto, sites, verified = transform(ver, t = NA), season = "t"),
               "The season column `t` has missing values.", fixed = TRUE)
})

test_that("fp_data_wide() checks the sizes of its companions", {
  ya <- matrix(c(1, 0, 1, 0, 0, 1), 3, 2)
  expect_error(fp_data_wide(ya, site_covs = data.frame(x = 1:2)),
               "`site_covs` must have one row per row of `y_auto`.", fixed = TRUE)
  expect_error(fp_data_wide(ya, obs_covs = list(matrix(1, 3, 2))),
               "`obs_covs` must be a named list of matrices.", fixed = TRUE)
  expect_error(fp_data_wide(ya, obs_covs = list(bio = matrix(1, 3, 2),
                                                matrix(2, 3, 2))),
               "`obs_covs` must be a named list of matrices.", fixed = TRUE)
  expect_error(fp_data_wide(ya, recorder = c("a", "b")),
               "`recorder` must have one value per site.", fixed = TRUE)
  d <- fp_data_wide(ya, site_covs = data.frame(x = 1:3),
                    obs_covs = list(bio = matrix(1, 3, 2)), recorder = c("a", "b", "a"))
  expect_equal(d$auto$recorder, c("a", "b", "a", "a", "b", "a"))
})

test_that("the model functions take only fp_data objects", {
  sim <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  msg <- "`data` must be an `fp_data` object; see `fp_data()`."
  expect_error(fp_occu(sim$data$auto), msg, fixed = TRUE)
  expect_error(fp_prior_predictive(unclass(sim$data)), msg, fixed = TRUE)
  expect_error(fp_occu(sim), msg, fixed = TRUE)
})

test_that("fp_data_wide() rejects site covariates that are not a table", {
  ya <- matrix(c(1, 0, 0, 1, 1, 0), 3, 2)
  expect_error(fp_data_wide(ya, site_covs = c(1, 2, 3)),
               "`site_covs` must be a data frame or matrix with one row per row of `y_auto`.",
               fixed = TRUE)
  expect_s3_class(fp_data_wide(ya, site_covs = cbind(x = 1:3)), "fp_data")
})

test_that("simulate_fp_ms() takes single finite numbers and no verified sites", {
  expect_error(simulate_fp_ms(nsite = c(20, 30)),
               "Argument(s) `nsite` must be single finite numbers.", fixed = TRUE)
  expect_error(simulate_fp_ms(psi = "0.5", b1 = NA),
               "Argument(s) `psi`, `b1` must be single finite numbers.", fixed = TRUE)
  expect_error(simulate_fp_ms(ngroup = Inf),
               "Argument(s) `ngroup` must be single finite numbers.", fixed = TRUE)
  s <- simulate_fp_ms(nspecies = 2, nsite = 20, nocc = 4, frac_verified = 0, seed = 1)
  expect_true(is.null(s$data$verified) || nrow(s$data$verified) == 0)
  expect_setequal(unique(s$data$auto$species), c("sp1", "sp2"))
})
