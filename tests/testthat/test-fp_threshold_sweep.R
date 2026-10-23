test_that("threshold metrics match hand-counted confusion matrices including ties", {
  v <- data.frame(confidence = c(0.9, 0.8, 0.8, 0.4),
                  verdict = c(FALSE, TRUE, FALSE, TRUE))
  m <- crywolf:::.threshold_metrics(v)
  expect_equal(m$threshold, c(0.4, 0.8, 0.9))
  expect_equal(m$tp, c(2, 1, 0))
  expect_equal(m$fp, c(2, 2, 1))
  expect_equal(m$fn, c(0, 1, 2))
  expect_equal(m$precision, c(1/2, 1/3, 0))
  expect_equal(m$recall, c(1, 1/2, 0))
  expect_equal(m$f1, c(2/3, 2/5, 0))
})

test_that("validation filters species and handles incomplete or invalid reviews", {
  v <- data.frame(confidence = c(0.9, 0.5, NA, 0.2),
                  verdict = c(TRUE, FALSE, TRUE, NA), species_code = "a")
  expect_warning(m <- crywolf:::.threshold_metrics(v), "2 validation row")
  expect_equal(m$threshold, c(0.5, 0.9))
  w <- rbind(v[1:2, ], transform(v[1:2, ], species_code = "b"))
  expect_error(crywolf:::.threshold_metrics(w), "several species")
  expect_equal(crywolf:::.threshold_metrics(w, "b"),
               crywolf:::.threshold_metrics(v[1:2, ]))
  expect_error(crywolf:::.threshold_metrics(w, "missing"), "at least one")
  expect_error(crywolf:::.threshold_metrics(v[1:2, 1:2], "a"), "species_code")
  expect_error(crywolf:::.threshold_metrics(NULL), "validation")
  expect_error(crywolf:::.threshold_metrics(transform(v[1:2, ], verdict = FALSE)),
               "at least one")
  for (bad in list(Inf, -0.1, 1.1, "0.9")) {
    expect_error(crywolf:::.threshold_metrics(data.frame(confidence = bad, verdict = TRUE)),
                 "confidence")
  }
  expect_error(crywolf:::.threshold_metrics(data.frame(confidence = 0.5, verdict = "positive")),
               "verdict")
  expect_equal(crywolf:::.threshold_metrics(data.frame(confidence = c(0, 1), verdict = c(0, 1)))$recall,
               c(1, 1))
})

# Exercise selection, aggregation and output without running a Stan sampler.
test_that("named selection preserves metrics and fits each threshold once", {
  fitted <- list()
  local_mocked_bindings(fp_occu = function(data, ...) {
    fitted[[length(fitted) + 1L]] <<- data$auto$y
    structure(list(draws = posterior::as_draws_matrix(matrix(seq_len(40) / 40,
                    ncol = 1, dimnames = list(NULL, "occ_(Intercept)"))),
                   diagnostics = list(divergences = 0, max_rhat = 1)), class = "fp_occu")
  }, .package = "crywolf")
  occasions <- data.frame(file = letters[1:4], site = 1:4)
  detections <- data.frame(file = letters[1:4], confidence = c(0.9, 0.8, 0.7, 0.6))
  sites <- data.frame(site = 1:4)
  v <- data.frame(confidence = detections$confidence, verdict = c(TRUE, FALSE, TRUE, FALSE))
  out <- fp_threshold_sweep(occasions, detections, sites = sites, validation = v)
  selected <- attr(out, "threshold_selection")
  expect_equal(selected$method, c("MaxF1", "MaxPrecision"))
  expect_equal(selected$threshold, c(0.7, 0.9))
  expect_equal(out$threshold, c(0.7, 0.9))
  expect_equal(out$detection_rate, c(0.75, 0.25))
  expect_length(fitted, 2)
  expect_s3_class(out, "fp_threshold_sweep")
  expect_equal(nrow(attr(out, "validation_metrics")), 4L)

  # Both metrics tie at precision 1 across scores: retain all positives.
  fitted <- list()
  expect_warning(out <- fp_threshold_sweep(occasions, detections, sites = sites,
                            validation = transform(v, verdict = TRUE)), "Every automated")
  expect_equal(attr(out, "threshold_selection")$threshold, c(0.6, 0.6))
  expect_length(fitted, 1)
  fitted <- list()
  old <- fp_threshold_sweep(occasions, detections, c(0.9, 0.7, 0.9), sites)
  expect_equal(old$threshold, c(0.7, 0.9))
  expect_null(attr(old, "threshold_selection"))
  expect_length(fitted, 2)
  for (bad in list(numeric(), character(), NA_real_, Inf, -1, 2, list(0.5))) {
    expect_error(fp_threshold_sweep(occasions, detections, bad, sites), "thresholds")
  }
  expect_error(fp_threshold_sweep(occasions, detections, "typo", sites, validation = v),
               "arg")
})

test_that("named thresholds without validation give a clear error", {
  expect_error(fp_threshold_sweep(data.frame(file = "a", site = 1),
                                  data.frame(file = "a", confidence = 0.5),
                                  sites = data.frame(site = 1)),
               "derived from reviewed clips")
})

test_that("empty or odd thresholds stop, and unused validation warns", {
  occ <- data.frame(file = "a", site = 1)
  det <- data.frame(file = "a", confidence = 0.5)
  sites <- data.frame(site = 1)
  expect_error(fp_threshold_sweep(occ, det, thresholds = character(0), sites = sites),
               "nonempty")
  expect_error(fp_threshold_sweep(occ, det, thresholds = 2, sites = sites), "\\[0, 1\\]")
  val <- data.frame(confidence = 0.5, verdict = TRUE)
  # the sites table lacks `site`, so the call stops right after the warning
  expect_warning(try(fp_threshold_sweep(occ, det, thresholds = 0.5,
                                        sites = data.frame(x = 1), validation = val),
                     silent = TRUE),
                 "ignored with numeric")
})
