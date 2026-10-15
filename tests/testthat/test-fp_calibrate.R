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
