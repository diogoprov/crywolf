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

test_that("calibration fits a recorder effect only when one was simulated", {
  seen <- NULL
  local_mocked_bindings(fp_occu = function(...) {
    seen <<- list(...)$random_recorder
    stop("not fitted")
  })
  score <- function(sim, fit_args = list()) {
    crywolf:::.score_fit(sim, fp = ~ 1, occ = ~ x, det = ~ 1, model = "m",
                         prob = 0.95, seed = 1, fit_args = fit_args)
  }
  score(simulate_fp(nsite = 20, nocc = 2, seed = 1))
  expect_false(seen)
  score(simulate_fp(nsite = 20, nocc = 2, sd_recorder = 0.5, seed = 1))
  expect_true(seen)
  score(simulate_fp(nsite = 20, nocc = 2, seed = 1),
        fit_args = list(random_recorder = TRUE))
  expect_true(seen)
})

test_that("fp_calibrate checks dynamic designs against simulate_fp_dyn()", {
  expect_error(fp_calibrate(data.frame(nseason = 3), reps = 1),
               "Unknown simulate_fp\\(\\)")
  expect_error(fp_calibrate(data.frame(conf_effect = 1), dynamic = TRUE, reps = 1),
               "Unknown simulate_fp_dyn\\(\\)")
  expect_error(fp_calibrate(data.frame(nsite = 50), dynamic = NA, reps = 1),
               "TRUE or FALSE")
})

test_that("dynamic calibration scores the transition parameters", {
  seen <- NULL
  local_mocked_bindings(fp_occu_dyn = function(...) {
    seen <<- list(...)
    stop("not fitted")
  })
  score <- function(sim, fp = ~ 1, col = ~ 1) {
    crywolf:::.score_fit_dyn(sim, fp = fp, psi1 = ~ x, col = col, ext = ~ 1,
                             det = ~ 1, model = "m", prob = 0.95, seed = 1,
                             fit_args = list())
  }
  sim <- simulate_fp_dyn(nsite = 20, nseason = 3, nocc = 2, col = 0.3,
                         b_col = 1, seed = 1)
  r <- score(sim, col = ~ rain)
  expect_false(seen$random_recorder)
  expect_equal(seen$col, ~ rain)
  expect_setequal(r$parameter, c("psi1_(Intercept)", "psi1_x", "col_(Intercept)",
                                 "col_rain", "ext_(Intercept)", "p_verified",
                                 "p10_bar"))
  expect_equal(r$truth[r$parameter == "col_rain"], 1)
  expect_equal(r$truth[r$parameter == "col_(Intercept)"], stats::qlogis(0.3))
  expect_equal(r$truth[r$parameter == "p10_bar"], 0.05)
  expect_false(any(r$ok))
  # no false-positive model, no verified data
  sim0 <- suppressWarnings(simulate_fp_dyn(nsite = 20, nseason = 3, nocc = 2,
                                           frac_verified = 0, seed = 1))
  r0 <- suppressWarnings(score(sim0, fp = NULL))
  expect_false(any(c("p_verified", "p10_bar", "col_rain") %in% r0$parameter))
})

test_that("fp_calibrate runs a tiny dynamic design", {
  skip_if_no_cmdstan()
  cal <- fp_calibrate(data.frame(nsite = 60, nseason = 3), dynamic = TRUE,
                      reps = 1,
                      fit_args = list(chains = 1, iter_warmup = 150,
                                      iter_sampling = 150))
  expect_s3_class(cal, "fp_calibration")
  expect_true(attr(cal, "dynamic"))
  expect_true(all(cal$ok))
  expect_true(all(c("col_(Intercept)", "ext_(Intercept)") %in% cal$parameter))
})

fake_calibration <- function() {
  grid <- expand.grid(rep = 1:20, frac_verified = c(0.1, 0.5, 0.9),
                      model = c("constant", "fp_x"), stringsAsFactors = FALSE)
  grid$cell <- match(grid$frac_verified, c(0.1, 0.5, 0.9))
  set.seed(1)
  grid$parameter <- "occ_x"
  grid$truth <- -0.8
  grid$mean <- -0.8 + ifelse(grid$model == "constant", 0.4, 0) + rnorm(nrow(grid), 0, 0.1)
  grid$lower <- grid$mean - 0.2
  grid$upper <- grid$mean + 0.2
  grid$covered <- grid$truth >= grid$lower & grid$truth <= grid$upper
  grid$divergences <- 0L; grid$max_rhat <- 1; grid$ok <- TRUE
  grid$seconds <- 1; grid$seed <- seq_len(nrow(grid)); grid$conf_effect <- 1.5
  structure(tibble::as_tibble(grid),
            class = c("fp_calibration", "tbl_df", "tbl", "data.frame"),
            prob = 0.95)
}

test_that("plot.fp_calibration draws coverage and bias", {
  skip_if_not_installed("ggplot2")
  cal <- fake_calibration()
  p <- plot(cal)
  expect_s3_class(p, "ggplot")
  expect_equal(sort(unique(p$data$x)), c(0.1, 0.5, 0.9))
  expect_true(all(p$data$lo >= 0 & p$data$hi <= 1))
  pb <- plot(cal, type = "bias")
  expect_s3_class(pb, "ggplot")
  expect_true(all(pb$data$y[pb$data$model == "constant"] > 0.2))
  expect_error(plot(cal, parameter = "nope"), "must be one of")
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("plot.fp_calibration picks a dynamic parameter by default", {
  skip_if_not_installed("ggplot2")
  cal <- fake_calibration()
  cal$parameter <- "col_rain"
  attr(cal, "dynamic") <- TRUE
  p <- plot(cal)
  expect_match(p$labels$y, "col_rain")
  ps <- plot(cal, by = "frac_verified", verified = "sites")
  # defaults of simulate_fp_dyn(): 100 sites, 5 seasons
  expect_equal(sort(unique(ps$data$x)), c(0.1, 0.5, 0.9) * 100 * 5)
  expect_match(ps$labels$x, "site-seasons")
})

test_that("fp_verification_curve passes dynamic designs to fp_calibrate()", {
  seen <- NULL
  local_mocked_bindings(fp_calibrate = function(design, ...) {
    seen <<- c(list(design = design), list(...))
    cal <- fake_calibration()
    cal$parameter <- ifelse(cal$rep %% 2 == 0, "col_(Intercept)", "ext_x")
    cal
  })
  out <- capture.output(
    cal <- fp_verification_curve(nseason = 4, dynamic = TRUE, ext = ~ x,
                                 frac_verified = c(0.1, 0.5, 0.9), reps = 2)
  )
  expect_true(seen$dynamic)
  expect_equal(seen$ext, ~ x)
  expect_equal(seen$design$nseason, rep(4, 3))
  expect_equal(seen$design$frac_verified, c(0.1, 0.5, 0.9))
  expect_true(any(grepl("site-seasons", out)))
  expect_true(any(grepl("col_\\(Intercept\\)", out)))
  expect_true(any(grepl("ext_x", out)))
  expect_s3_class(cal, "fp_calibration")
})
