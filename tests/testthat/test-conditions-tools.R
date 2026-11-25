# Conditions of the calibration, threshold, import, export and plotting
# tools that no other test file triggers. None of these tests fits a model.

cond_calibration <- function() {
  grid <- expand.grid(rep = 1:5, frac_verified = c(0.1, 0.5),
                      model = "constant", stringsAsFactors = FALSE)
  grid$cell <- match(grid$frac_verified, c(0.1, 0.5))
  grid$parameter <- "occ_x"
  grid$truth <- -0.8
  grid$mean <- -0.8
  grid$lower <- -1
  grid$upper <- -0.6
  grid$covered <- TRUE
  grid$divergences <- 0L; grid$max_rhat <- 1; grid$ok <- TRUE
  grid$seconds <- 1; grid$seed <- seq_len(nrow(grid))
  structure(tibble::as_tibble(grid),
            class = c("fp_calibration", "tbl_df", "tbl", "data.frame"),
            prob = 0.95)
}

#' @srrstats {G5.2b} Conditions raised by the calibration, threshold, import, export and plotting tools.
test_that("plot.fp_calibration rejects an unknown axis and a misplaced `verified`", {
  skip_if_not_installed("ggplot2")
  cal <- cond_calibration()
  expect_error(plot(cal, by = "nope"), "`by` must be a column of the design.",
               fixed = TRUE)
  expect_error(plot(cal, by = "model", verified = "sites"),
               "`verified` applies only when the axis is `frac_verified`.",
               fixed = TRUE)
  # the same call on the right axis works
  expect_s3_class(plot(cal, by = "frac_verified", verified = "sites"), "ggplot")
})

test_that("missing optional packages stop the calibration tools", {
  local_mocked_bindings(requireNamespace = function(...) FALSE, .package = "base")
  expect_error(plot(cond_calibration()), "plot() needs the 'ggplot2' package.",
               fixed = TRUE)
  jobs <- data.frame(x = 1)
  expect_error(crywolf:::.run_jobs(jobs, function(x) x, parallel = TRUE),
               "`parallel = TRUE` needs the 'furrr' package.", fixed = TRUE)
})

test_that(".require_pkg() names the missing package", {
  expect_error(crywolf:::.require_pkg("crywolfNoSuchPackage"),
               "This function needs the 'crywolfNoSuchPackage' package.",
               fixed = TRUE)
  expect_true(crywolf:::.require_pkg("stats"))
})

test_that("fp_calibrate_fit() and fp_ppc() reject objects that are not fits", {
  expect_error(fp_calibrate_fit(list(), models = list(a = ~ 1)),
               "`fit` must be an `fp_occu` object.", fixed = TRUE)
  expect_error(fp_ppc(list()),
               "`object` must be an `fp_occu`, `fp_occu_dyn`, `fp_occu_ms` or `fp_occu_rn` fit.",
               fixed = TRUE)
})

test_that("fp_ppc() needs verified occasions for method = \"verified\"", {
  fake <- structure(list(stan_data = list(N1 = 0L)), class = "fp_occu")
  expect_error(fp_ppc(fake, method = "verified"),
               "The data have no verified occasions, so `method = \"verified\"`",
               fixed = TRUE)
})

test_that("fp_to_unmarked() checks the observation covariates and occasions", {
  skip_if_not_installed("unmarked")
  d <- simulate_fp(nsite = 20, nocc = 3, seed = 1)$data
  d1 <- d
  d1$auto$METH <- 1
  expect_error(fp_to_unmarked(d1, obs_covs = "METH"),
               "`METH` is the method covariate made by fp_to_unmarked()",
               fixed = TRUE)
  d2 <- d
  d2$auto$habitat <- "forest"
  d2$auto$noise <- 0
  expect_error(fp_to_unmarked(d2, obs_covs = c("noise", "habitat")),
               "Observation covariates must be numeric: habitat.", fixed = TRUE)
  d3 <- d
  d3$auto$visit <- 1L
  expect_error(fp_to_unmarked(d3, occasion = "visit"),
               "Two occasions of a site have the same `visit`.", fixed = TRUE)
})

test_that("fp_read_birdnet() stops without tables or without detections", {
  empty_dir <- tempfile()
  dir.create(empty_dir)
  expect_error(fp_read_birdnet(empty_dir), "No selection table found.", fixed = TRUE)
  tab <- file.path(empty_dir, "rec.txt")
  writeLines(paste("Selection", "Begin Time (s)", "End Time (s)", "Begin Path",
                   "Confidence", sep = "\t"), tab)
  expect_error(fp_read_birdnet(tab), "The selection tables hold no detection.",
               fixed = TRUE)
})

test_that("fp_aggregate() warns about reviewed clips outside the occasions", {
  occ <- data.frame(site = "A", file = "a1")
  det <- data.frame(file = "a1", species_code = "sp", confidence = 0.9)
  val <- data.frame(file = c("a1", "zz"), species_code = "sp",
                    confidence = c(0.9, 0.8), verdict = c(TRUE, TRUE))
  expect_warning(out <- fp_aggregate(occ, det, threshold = 0.5, validation = val),
                 "1 occasion(s) with reviewed clips are not in `occasions`",
                 fixed = TRUE)
  expect_equal(out$n_true, 1L)
  expect_equal(out$n_false, 0L)
})

test_that("plot.fp_occu() checks `pars` against the draws", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("bayesplot")
  fake <- structure(
    list(draws = posterior::as_draws_df(data.frame(`occ_(Intercept)` = 0,
                                                   check.names = FALSE))),
    class = "fp_occu")
  expect_error(plot(fake, type = "trace"),
               "The model has none of the default parameters; name them in `pars`.",
               fixed = TRUE)
  expect_error(plot(fake, type = "trace", pars = c("occ_(Intercept)", "nope")),
               "Unknown parameter(s): nope.", fixed = TRUE)
})

test_that("fp_effects() needs a false-positive formula for process = \"fp\"", {
  fake <- structure(list(formulas = list(occ = ~ x)), class = "fp_occu")
  expect_error(fp_effects(fake, "x", process = "fp"),
               "This model has no false-positive process to plot.", fixed = TRUE)
})

test_that("fp_sites() refuses to overwrite columns of the sites table", {
  sites <- data.frame(site = 1, psi = 0.3, p10 = 0.1)
  expect_error(
    crywolf:::.sites_table(sites, sd = list(X = matrix(1)),
                           pars = list(beta = matrix(0)), pz = matrix(0.5),
                           prob = 0.95),
    "The sites table already has column(s) psi, p10; rename them",
    fixed = TRUE)
})
