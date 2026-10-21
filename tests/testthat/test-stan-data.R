# These tests build the Stan data list without running Stan.

test_that("design matrices and indices are built correctly", {
  s <- simulate_fp(nsite = 40, nocc = 3, seed = 1)
  b <- crywolf:::.build_stan_data(s$data, occ = ~ x, det = ~ bio, fp = ~ x,
                                  random_recorder = TRUE, priors = fp_priors())
  sd <- b$stan_data
  expect_equal(sd$Kbeta, 2)
  expect_equal(sd$Kalpha, 2)
  expect_equal(sd$Kfp, 1)
  expect_equal(dim(sd$Xfp), c(40, 1))
  expect_equal(sd$N2, 120)          # bio varies by occasion: no collapsing
  expect_equal(sum(sd$n2), 120)
  expect_equal(sd$estimate_re, 1L)
  expect_equal(sd$nrec, length(unique(s$data$auto$recorder)))
  expect_true(all(sd$Z[, 1] == 1))
  expect_equal(b$design$fp_names, c("(Intercept)", "x"))
})

test_that("constant and no-false-positive models", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  b0 <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  expect_equal(b0$stan_data$Kfp, 0)
  expect_equal(b0$stan_data$estimate_re, 0L)
  bn <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, NULL, FALSE, fp_priors())
  expect_equal(bn$stan_data$estimate_fp, 0L)
})

test_that("detection formulas may use site covariates", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  b <- crywolf:::.build_stan_data(s$data, ~ 1, ~ x + bio, ~ 1, FALSE, fp_priors())
  expect_equal(b$stan_data$Kalpha, 3)
})

test_that("bad formulas fail early", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  expect_error(crywolf:::.build_stan_data(s$data, ~ nope, ~ 1, ~ 1, FALSE, fp_priors()),
               "not found")
  expect_error(crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ x - 1, FALSE, fp_priors()),
               "intercept")
  expect_error(crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, TRUE,
                                          fp_priors())$stan_data, NA)
})

test_that("parameter names map to the Stan output", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  b <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ w, TRUE, fp_priors())
  m <- crywolf:::.name_map(b$design, TRUE, TRUE)
  expect_true(all(c("occ_x", "fp_w", "fp_(Intercept)", "sd_recorder",
                    "p10_bar") %in% m))
  expect_equal(names(m)[m == "fp_w"], "gamma_fp[1]")
})

test_that("the constant-rate model has no slope names", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  b <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  m <- crywolf:::.name_map(b$design, TRUE, FALSE)
  expect_false(anyNA(names(m)))
  expect_false(any(grepl("gamma_fp", names(m))))
})

test_that("automated rows are sorted by site with correct segments", {
  sites <- data.frame(site = c("a", "b", "c"))
  auto <- data.frame(site = c("c", "a", "c", "a", "a"), y = c(1, 0, 0, 1, 0))
  # site b has no automated occasion on purpose (len2 = 0)
  d <- suppressWarnings(fp_data(auto, sites))
  b <- crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors())$stan_data
  # no detection covariates: one binomial row per site
  expect_equal(b$site2, c(1L, 3L))
  expect_equal(b$n2, c(3L, 2L))
  expect_equal(b$y2, c(1L, 1L))
  expect_equal(b$len2, c(1L, 0L, 1L))
  expect_equal(b$start2, c(1L, 2L, 2L))
})

test_that("collapsing keeps rows that differ in covariates or recorder", {
  sites <- data.frame(site = 1:2)
  auto <- data.frame(site = c(1, 1, 1, 2, 2), y = c(1, 0, 1, 0, 0),
                     v = c(0, 0, 1, 0, 0), rec = c("a", "a", "a", "a", "b"))
  d <- fp_data(auto, sites, recorder = "rec")
  b <- crywolf:::.build_stan_data(d, ~ 1, ~ v, ~ 1, TRUE, fp_priors())$stan_data
  expect_equal(b$N2, 4)
  expect_equal(b$n2, c(2L, 1L, 1L, 1L))
  expect_equal(b$y2, c(1L, 1L, 0L, 0L))
  expect_equal(b$Z[, 2], c(0, 1, 0, 0))
})

test_that("confirmed occasions leave the binomial rows", {
  s <- simulate_fp(nsite = 30, nocc = 4, frac_confirmed = 0.5, seed = 2)
  a <- s$data$auto
  b <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  sd <- b$stan_data
  conf <- !is.na(a$n_true)
  expect_equal(sd$obs_confirm, 1L)
  expect_equal(sd$N3, sum(conf))
  expect_equal(sum(sd$n2) + sd$N3, nrow(a))
  expect_equal(sum(sd$y2), sum(a$y[!conf]))
  expect_equal(sort(sd$v3), sort(2L * (a$n_true[conf] > 0) + (a$n_false[conf] > 0)))
  expect_equal(dim(sd$Z3), c(sum(conf), 1L))
  expect_error(crywolf:::.build_stan_data(s$data, ~ x, ~ 1, NULL, FALSE, fp_priors()),
               "needs a false-positive")
  m <- crywolf:::.name_map(b$design, TRUE, FALSE, has_verified = FALSE)
  expect_false("p_verified" %in% m)
})

test_that("the site random effect on false positives is switched on and named", {
  s <- simulate_fp(nsite = 30, nocc = 4, seed = 1)
  b0 <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors())
  expect_equal(b0$stan_data$estimate_fp_re, 0L)
  b1 <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE,
                                   fp_priors(fp_site_sd = 0.5), random_fp = TRUE)
  expect_equal(b1$stan_data$estimate_fp_re, 1L)
  expect_equal(b1$stan_data$fp_site_sd, 0.5)
  expect_error(crywolf:::.build_stan_data(s$data, ~ x, ~ 1, NULL, FALSE,
                                          fp_priors(), random_fp = TRUE),
               "needs a false-positive process")
  m <- crywolf:::.name_map(b1$design, TRUE, FALSE, use_fp_re = TRUE)
  expect_equal(m[["sigma_fp[1]"]], "sd_fp_site")
  expect_false("sigma_fp[1]" %in% names(crywolf:::.name_map(b1$design, TRUE, FALSE)))
  expect_output(print(fp_priors()), "sd_fp_site")
})
