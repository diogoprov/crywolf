# Conditions that no other test raised, found by crossing the lines of each
# stop(), warning() and message() with line coverage. They use small fake
# objects or mocked internals, so that no model is fitted.

#' @srrstats {G5.2b} Every remaining error, warning and message is triggered and its text checked.
test_that("the calibration functions need named models", {
  msg <- paste("`models` must be a named list of formulas, one per model fitted by",
               "fp_calibrate().")
  expect_error(fp_calibrate(list(nsite = 20), models = list(~ 1)), msg, fixed = TRUE)
  fit <- structure(list(formulas = list(fp = ~ 1)), class = "fp_occu")
  msg <- paste("`models` must be a named list of formulas, one per model refitted by",
               "fp_calibrate_fit().")
  expect_error(fp_calibrate_fit(fit, models = list(~ 1)), msg, fixed = TRUE)
})

test_that("a recorder effect needs a recorder column", {
  d <- suppressWarnings(simulate_fp(nsite = 20, nocc = 3, seed = 1)$data)
  d$recorder <- NULL
  expect_error(crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, TRUE, fp_priors()),
               "`random_recorder = TRUE` but `data` has no recorder column.",
               fixed = TRUE)
})

test_that("fp_read_sheepdog() refuses other exports", {
  skip_if_not_installed("jsonlite")
  dir <- tempfile("sheepdog-other")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  writeLines('{"schema": "sheepdog-training"}', file.path(dir, "schema.json"))
  expect_error(fp_read_sheepdog(dir),
               "`path` holds a 'sheepdog-training' export; only 'sheepdog-occupancy' is supported.",
               fixed = TRUE)
})

test_that("fp_sbc() checks its ranks before fitting", {
  s <- suppressWarnings(simulate_fp(nsite = 20, nocc = 3, seed = 1))
  expect_error(fp_sbc(s$data, occ = ~ x, random_recorder = FALSE, reps = 1,
                      ndraws = 99, fit_args = list(chains = 1, iter_sampling = 50)),
               "`ndraws` exceeds the number of posterior draws of each fit (chains x iter_sampling).",
               fixed = TRUE)
  expect_error(crywolf:::.check_bins(99, 7),
               "`bins` must divide the number of possible ranks, ndraws + 1 = 100, so that every bin is equally likely.",
               fixed = TRUE)
})

test_that("loo() refuses draws from the priors alone", {
  x <- structure(list(prior_only = TRUE), class = "fp_occu")
  expect_error(loo::loo(x),
               "`x` holds draws from the priors alone (prior_only = TRUE); leave-one-out needs a fit to the data.",
               fixed = TRUE)
})

test_that("fp_effects() and plot.fp_sites() check their arguments", {
  m <- structure(list(formulas = list(occ = ~ x, fp = ~ 1)), class = "fp_occu")
  expect_error(fp_effects(m, "rain"), "`rain` is not in the occ formula.", fixed = TRUE)
  skip_if_not_installed("ggplot2")
  x <- structure(data.frame(site = 1:2, verified = c(TRUE, FALSE), detections = c(2, 0),
                            lon = 1:2, lat = 1:2),
                 class = c("fp_sites", "data.frame"))
  expect_error(plot(x, coords = "lon"), "`coords` must name two columns of `x`.",
               fixed = TRUE)
  expect_error(plot(x, coords = c("lon", "height")),
               "`coords` must name two columns of `x`.", fixed = TRUE)
})

# a fake fit with the fields fp_occu_mc() reads, and a given largest R-hat
fake_mc_fit <- function(rhat) {
  draws <- posterior::as_draws_array(array(stats::rnorm(10), c(10, 1, 1),
                                           dimnames = list(NULL, NULL, "a")))
  structure(list(draws = draws,
                 diagnostics = list(divergences = 0, max_rhat = rhat,
                                    min_ess_bulk = 100, min_ess_tail = 100)),
            class = "fp_occu")
}

test_that("fp_occu_mc() reports failed and unconverged fits", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  outcome <- list("ok", "error", "slow")
  fit_no <- 0L
  local_mocked_bindings(
    fp_compile = function(...) invisible(NULL),
    fp_occu = function(...) {
      fit_no <<- fit_no + 1L
      switch(outcome[[fit_no]],
             ok = fake_mc_fit(1),
             error = stop("the sampler broke"),
             slow = fake_mc_fit(1.5))
    },
    summary.fp_occu = function(object, ...) list(table = data.frame(variable = "a"))
  )
  w <- character()
  mc <- withCallingHandlers(
    fp_occu_mc(list(s$data, s$data, s$data)),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(startsWith(w, paste("1 fit(s) failed and were left out (version(s) 2);",
                                      "the first error: the sampler broke"))))
  expect_true(any(startsWith(w, paste("1 fit(s) did not converge (largest R-hat above 1.05)",
                                      "and were left out: version(s) 3."))))
  expect_equal(mc$diagnostics$kept, c(TRUE, FALSE))

  fit_no <- 0L
  outcome <- list("slow")
  expect_error(fp_occu_mc(list(s$data)),
               "No fit converged (largest R-hat above 1.05 in every fit); run longer chains.",
               fixed = TRUE)
})

test_that("fp_occu_mc() announces the parametrization of multispecies versions", {
  s <- simulate_fp_ms(nspecies = 2, nsite = 20, nocc = 3, seed = 1)
  local_mocked_bindings(
    .compile_model = function(...) invisible(NULL),
    fp_occu_ms = function(...) stop("not fitted")
  )
  expect_message(
    expect_error(fp_occu_mc(list(s$data)), "Every fit failed"),
    "Parametrization of every version \\(2 species at 20 sites\\): centered = .*see \\?fp_occu_ms"
  )
})
