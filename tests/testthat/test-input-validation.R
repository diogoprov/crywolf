sim_data <- function() {
  suppressWarnings(simulate_fp(nsite = 40, nocc = 5, seed = 3)$data)
}

#' @srrstats {G2.4e, G2.6} Factor detections are read by their labels.
test_that("a factor detection column is read by its labels", {
  d <- sim_data()
  a <- d$auto
  a$y <- factor(a$y, levels = c(1, 0))
  d2 <- suppressWarnings(fp_data(a, d$sites, d$verified))
  expect_identical(d2$auto$y, as.integer(d$auto$y))
  sl_f <- crywolf:::.build_stan_data(d2, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  sl_n <- crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  expect_false(is.null(sl_n$stan_data$y2))
  expect_identical(sl_f$stan_data$y2, sl_n$stan_data$y2)
})

#' @srrstats {G5.8b} Detections of unsupported types or values.
test_that("logical and character 0/1 detections are accepted; others are not", {
  d <- sim_data()
  a <- d$auto
  a$y <- as.logical(a$y)
  expect_identical(suppressWarnings(fp_data(a, d$sites))$auto$y, as.integer(d$auto$y))
  a$y <- as.character(d$auto$y)
  expect_identical(suppressWarnings(fp_data(a, d$sites))$auto$y, as.integer(d$auto$y))
  a$y <- ifelse(d$auto$y == 1, "yes", "no")
  expect_error(fp_data(a, d$sites), "coded 0/1")
  a$y <- d$auto$y + 1L
  expect_error(fp_data(a, d$sites), "coded 0/1")
})

#' @srrstats {G5.8, G5.8a, G5.8c} Zero-row and all-NA detections.
test_that("automated data without any detection record are rejected", {
  d <- sim_data()
  expect_error(fp_data(d$auto[0, ], d$sites), "no automated occasion")
  a <- d$auto
  a$y <- NA
  expect_error(fp_data(a, d$sites), "no automated occasion")
})

#' @srrstats {G2.16} Infinite covariates.
test_that("infinite covariates are rejected", {
  d <- sim_data()
  d$sites$x[1] <- Inf
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, fp_priors()),
               "Infinite values in the covariates of `occ`")
})

#' @srrstats {BS2.3, BS2.5} Lengths and values of prior parameters.
test_that("fp_priors() takes single positive numbers", {
  expect_error(fp_priors(intercept_sd = "a"), "`intercept_sd` must be a single positive")
  expect_error(fp_priors(slope_sd = c(1, 2)), "`slope_sd` must be a single positive")
  expect_error(fp_priors(gap_sd = Inf), "`gap_sd` must be a single positive")
  expect_error(fp_priors(recorder_sd = 0), "`recorder_sd` must be a single positive")
  expect_error(fp_priors(gap_mean = NA), "`gap_mean` must be a single finite")
  expect_s3_class(fp_priors(gap_mean = -1), "fp_priors")
})

#' @srrstats {BS2.2} Priors are validated before fitting.
test_that("priors and logical options are checked before fitting", {
  d <- sim_data()
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, list(slope_sd = 1)),
               "`priors` must come from fp_priors")
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, "yes", fp_priors()),
               "`random_recorder` must be TRUE or FALSE")
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, fp_priors(),
                                          random_fp = NA),
               "`random_fp` must be TRUE or FALSE")
  bad <- fp_priors()
  bad$gap_sd <- -1
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, bad),
               "`gap_sd` must be a single positive")
  old <- fp_priors()
  old$species_sd <- NULL
  expect_silent(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, old))
})

test_that("simulate_fp() takes single finite numbers", {
  expect_error(simulate_fp(nsite = c(20, 30)), "`nsite`")
  expect_error(simulate_fp(psi = "0.5"), "`psi`")
  expect_error(simulate_fp(b1 = NA), "`b1`")
})
