test_that("prior draws respect the priors and the gap", {
  s <- simulate_fp(nsite = 30, nocc = 4, seed = 1)
  sd <- crywolf:::.build_stan_data(s$data, ~ x, ~ bio, ~ x, TRUE, fp_priors(),
                                   random_fp = TRUE)$stan_data
  set.seed(1)
  d <- purrr::map(1:500, function(i) crywolf:::.prior_draw(sd, fp_priors()))
  expect_length(d[[1]]$beta, 2)
  expect_length(d[[1]]$u, sd$nrec)
  expect_true(all(purrr::map_lgl(d, function(p) p$gamma0 < p$alpha[1])))
  expect_true(all(purrr::map_dbl(d, "sigma_fp") > 0))
  expect_equal(sd(purrr::map_dbl(d, function(p) p$beta[2])), 2.5, tolerance = 0.1)
})

test_that("fp_prior_predictive returns summaries and the observed values", {
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 2)
  pp <- fp_prior_predictive(s$data, occ = ~ x, fp = ~ 1,
                            random_recorder = FALSE, ndraws = 20, seed = 1)
  expect_s3_class(pp, "fp_prior_predictive")
  expect_equal(nrow(pp), 20)
  expect_true(all(c("naive_occupancy", "automated_rate", "verified_rate",
                    "mean_psi", "p10_bar") %in% names(pp)))
  obs <- attr(pp, "observed")
  expect_equal(obs[["automated_rate"]], mean(s$data$auto$y))
  expect_true(all(pp$naive_occupancy >= 0 & pp$naive_occupancy <= 1))
  expect_output(print(pp), "share_below_observed")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(pp), "ggplot")
})

test_that("the calibration plot can show verification effort as a count", {
  skip_if_not_installed("ggplot2")
  cal <- tibble::as_tibble(data.frame(
    cell = rep(1:2, each = 4), model = "m", parameter = "occ_x", truth = -0.8,
    mean = -0.8, lower = -1, upper = -0.6, covered = TRUE, divergences = 0L,
    max_rhat = 1, ok = TRUE, seconds = 1, rep = rep(1:4, 2), seed = 1:8,
    nsite = 200, frac_verified = rep(c(0.1, 0.5), each = 4)))
  cal <- structure(cal, class = c("fp_calibration", class(cal)), prob = 0.95)
  p <- plot(cal, verified = "occasions", target = 0.9)
  expect_s3_class(p, "ggplot")
  expect_equal(sort(unique(p$data$x)), c(0.1, 0.5) * 200 * 3)
  expect_match(p$labels$x, "verified occasions")
})

test_that("plots and site estimates work on a fit, also after reloading it", {
  skip_if_no_cmdstan()
  skip_if_not_installed("ggplot2")
  s <- simulate_fp(nsite = 80, nocc = 6, conf_effect = 1, conf_cor = 0.5,
                   seed = 3)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ x,
                                random_recorder = FALSE, chains = 2,
                                iter_warmup = 200, iter_sampling = 200,
                                seed = 1, show_messages = FALSE))
  f <- tempfile(fileext = ".rds")
  saveRDS(m, f)
  unlink(m$fit$output_files())          # CmdStan's files are gone
  m <- readRDS(f)

  expect_s3_class(plot(m), "ggplot")
  skip_if_not_installed("bayesplot")
  expect_s3_class(plot(m, type = "trace"), "ggplot")

  e <- fp_effects(m, "x", transform = function(z) 2 * z + 10)
  expect_s3_class(e, "fp_effects")
  expect_equal(e$x, 2 * e$value + 10)
  expect_true(all(e$lower <= e$estimate & e$estimate <= e$upper))
  expect_error(fp_effects(m, "w"), "not in the occ formula")
  expect_s3_class(plot(e), "ggplot")

  st <- fp_sites(m, ndraws = 200, seed = 1)
  expect_equal(nrow(st), 80)
  expect_true(all(st$occupied[st$verified] == 1))
  expect_true(all(st$occupied >= 0 & st$occupied <= 1))
  expect_equal(sum(st$detections), sum(s$data$auto$y))
  expect_true(all(st$p10 > 0 & st$p10 < 1))
  st$lon <- runif(80)
  st$lat <- runif(80)
  expect_s3_class(plot(st), "ggplot")
  expect_s3_class(plot(st, coords = c("lon", "lat")), "ggplot")
  expect_error(plot(st, coords = "lon"), "two columns")

  p <- fp_ppc(m, ndraws = 100, seed = 1)
  expect_s3_class(plot(p), "ggplot")
  expect_s3_class(plot(p, type = "counts"), "ggplot")
  expect_s3_class(plot(p, type = "histogram"), "ggplot")
})

test_that("simulating from data does not repeat the design warnings", {
  s <- simulate_fp(nsite = 30, nocc = 4, seed = 4)
  sites <- rbind(s$data$sites, data.frame(site = 31, x = 0, w = 0))
  d <- suppressWarnings(fp_data(s$data$auto, sites, s$data$verified,
                                recorder = "recorder"))
  expect_no_warning(fp_prior_predictive(d, occ = ~ x, fp = ~ 1,
                                        random_recorder = FALSE, ndraws = 10,
                                        seed = 1))
})
