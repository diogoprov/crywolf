test_that("the Markov lag follows the occasion order and resets at each site", {
  auto <- data.frame(site = c(1, 1, 1, 2, 2), t = c(3, 1, 2, 2, 1),
                     y = c(0, 1, 1, 1, 0))
  d <- fp_data(auto, data.frame(site = 1:2))
  b <- crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                  markov = "t")
  # site 1 in time order: y = 1, 1, 0 -> lag 0, 1, 1; site 2: 0, 1 -> 0, 0
  lag <- b$occasions$Z[, "y_prev"]
  expect_equal(unname(lag), c(0, 1, 1, 0, 0))
  expect_true("y_prev" %in% b$design$det_names)
  expect_error(crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                          markov = "nope"), "lacks column")
})

test_that("the Markov lag is not carried across gaps or ties", {
  # occasion 3 of site 1 is missing: occasion 4 starts a new run
  auto <- data.frame(site = c(1, 1, 1, 2, 2), t = c(1, 2, 4, 1, 2),
                     y = c(1, 1, 1, 1, 0))
  d <- fp_data(auto, data.frame(site = 1:2))
  b <- crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                  markov = "t")
  expect_equal(unname(b$occasions$Z[, "y_prev"]), c(0, 1, 0, 0, 1))
  expect_equal(b$occasions$has_prev, c(FALSE, TRUE, FALSE, FALSE, TRUE))
  # date-times: the step is the shortest interval (15 minutes)
  auto$t <- as.POSIXct("2024-01-01 00:00", tz = "UTC") + 900 * auto$t
  d <- fp_data(auto, data.frame(site = 1:2))
  b <- crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                  markov = "t")
  expect_equal(unname(b$occasions$Z[, "y_prev"]), c(0, 1, 0, 0, 1))
  auto$t <- c(1, 1, 2, 1, 2)
  d <- fp_data(auto, data.frame(site = 1:2))
  expect_error(crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                          markov = "t"), "same `markov` value")
})

test_that("simulated runs restart after a gap", {
  s <- simulate_fp(nsite = 200, nocc = 6, seed = 5)
  d <- s$data
  d$auto <- d$auto[d$auto$occasion != 3, ]
  sl <- crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, fp_priors(),
                                   markov = "occasion")
  gen <- list(data = d, formulas = list(occ = ~ x, det = ~ 1, fp = ~ 1, ver = ~ 1),
              stan_data = sl$stan_data, design = sl$design, occasions = sl$occasions)
  # an enormous run effect: after a detection the next adjacent occasion is
  # always a detection, but occasion 4 does not inherit occasion 2
  p <- list(beta = c(5, 0), alpha = c(-3, 50), delta0 = 0, delta_s = numeric(0),
            gamma0 = -20, gamma_fp = numeric(0))
  set.seed(1)
  a <- crywolf:::.simulate_from(gen, p)$data$auto
  a <- a[order(a$site, a$occasion), ]
  y <- matrix(a$y, ncol = 5, byrow = TRUE)          # occasions 1, 2, 4, 5, 6
  expect_true(all(y[y[, 1] == 1, 2] == 1))
  expect_true(any(y[, 2] == 1 & y[, 3] == 0))
  expect_true(all(y[y[, 3] == 1, 4] == 1))
})

test_that("occasion-level false-positive covariates switch to rows", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  s$data$auto$noise <- seq_len(nrow(s$data$auto))
  b0 <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ x, FALSE, fp_priors())
  expect_equal(b0$stan_data$fp_by_row, 0L)
  b1 <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ noise + x, FALSE, fp_priors())
  sd <- b1$stan_data
  expect_equal(sd$fp_by_row, 1L)
  expect_equal(dim(sd$Xfp2), c(sd$N2, 2L))
  expect_equal(ncol(sd$Xfp), 0L)
  # a distinct noise per occasion keeps every occasion as its own row
  expect_equal(sd$N2, nrow(s$data$auto))
})

test_that("verified covariates are sorted with their sites and named", {
  s <- simulate_fp(nsite = 30, nocc = 3, frac_verified = 0.5, seed = 2)
  s$data$verified$who <- rnorm(nrow(s$data$verified))
  b <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE, fp_priors(),
                                  ver = ~ who)
  sd <- b$stan_data
  expect_equal(sd$Kver_s, 1L)
  expect_false(is.unsorted(sd$site1))
  expect_equal(sort(sd$X1s[, 1]), sort(s$data$verified$who))
  m <- crywolf:::.name_map(b$design, TRUE, FALSE)
  expect_equal(m[["delta_s[1]"]], "verified_who")
  expect_error(crywolf:::.build_stan_data(s$data, ~ x, ~ 1, ~ 1, FALSE,
                                          fp_priors(), ver = ~ 0 + who),
               "keep the intercept")
})

test_that("fp_gap_prior recovers the gap from the precision", {
  set.seed(1)
  psi <- 0.4; p11 <- 0.3; p10 <- 0.03
  # share of true detections when false alarms also occur at occupied sites
  s <- (p11 - p10) / (1 - p10)
  prec <- psi * s / (psi * p11 + (1 - psi) * p10)
  n <- 1e5
  pr <- fp_gap_prior(round(prec * n), n - round(prec * n), p11 = p11, psi = psi)
  expect_s3_class(pr, "fp_priors")
  expect_equal(pr$gap_mean, qlogis(p11) - qlogis(p10), tolerance = 1e-3)
  expect_true(pr$gap_sd > 0)
  # a poor classifier implies a false-positive rate close to p11: a small gap
  low <- fp_gap_prior(1, 40, p11 = 0.1, psi = 0.1)
  expect_true(low$gap_mean > 0 && low$gap_mean < 1)
})

test_that("the new options fit, match the R likelihood, and refit exactly", {
  skip_if_no_cmdstan()
  set.seed(3)
  s <- simulate_fp(nsite = 80, nocc = 6, frac_verified = 0.4, seed = 3)
  d <- s$data
  d$auto$noise <- rnorm(nrow(d$auto))
  d$verified$who <- rnorm(nrow(d$verified))
  m <- suppressWarnings(fp_occu(d, occ = ~ x, fp = ~ noise, ver = ~ who,
                                markov = "occasion", random_recorder = FALSE,
                                chains = 1, iter_warmup = 200, iter_sampling = 150,
                                threads_per_chain = 1, seed = 1, show_messages = FALSE))
  expect_true(all(c("fp_noise", "verified_who", "det_y_prev") %in%
                    posterior::variables(m$draws)))
  pars <- crywolf:::.ppc_pars(m$fit, m$stan_data, NULL)
  r <- crywolf:::.ppc_core(m$stan_data, pars, "automated", "freeman-tukey")
  ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
  expect_lt(max(abs(r$ll - ll)), 1e-4)
  expect_output(print(m), "y_prev")
  expect_no_error(simulate(m, seed = 1))
  pf <- predict(m, process = "fp")
  expect_equal(nrow(pf), nrow(d$auto))
  expect_s3_class(fp_sites(m, ndraws = 50, seed = 1), "fp_sites")
  expect_s3_class(fp_ppc(m, method = "verified", ndraws = 50, seed = 1), "fp_ppc")

  # observation confirmation with an occasion covariate and runs
  s2 <- simulate_fp(nsite = 60, nocc = 5, frac_confirmed = 0.4, seed = 6)
  d2 <- s2$data
  d2$auto$noise <- rnorm(nrow(d2$auto))
  m2 <- suppressWarnings(fp_occu(d2, occ = ~ x, fp = ~ noise, markov = "occasion",
                                 random_recorder = FALSE, chains = 1,
                                 iter_warmup = 200, iter_sampling = 100,
                                 threads_per_chain = 1, seed = 1, show_messages = FALSE))
  pars2 <- crywolf:::.ppc_pars(m2$fit, m2$stan_data, NULL)
  r2 <- crywolf:::.ppc_core(m2$stan_data, pars2, "automated", "freeman-tukey")
  ll2 <- unclass(posterior::as_draws_matrix(m2$fit$draws("log_lik")))
  expect_lt(max(abs(r2$ll - ll2)), 1e-4)

  m0 <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1,
                                 random_recorder = FALSE, chains = 2,
                                 iter_warmup = 300, iter_sampling = 300,
                                 threads_per_chain = 1, seed = 1, show_messages = FALSE))
  l0 <- loo(m0)
  thr <- sort(l0$diagnostics$pareto_k, decreasing = TRUE)[2]
  l1 <- suppressMessages(loo(m0, reloo = TRUE, k_threshold = thr, seed = 2))
  i <- which.max(l0$diagnostics$pareto_k)
  expect_equal(l1$diagnostics$pareto_k[i], 0)
  expect_equal(unname(l1$pointwise[i, "mcse_elpd_loo"]), 0)
  # a fit saved before 1.5.0 has no sampler settings: defaults are used
  m_old <- m0
  m_old$sampler <- NULL
  expect_no_error(suppressMessages(loo(m_old, reloo = TRUE, k_threshold = thr, seed = 2)))
  expect_equal(l1$pointwise[i, "elpd_loo"], l0$pointwise[i, "elpd_loo"],
               tolerance = 0.05)
  w <- fp_loo_weights(a = m0, b = l0)
  expect_equal(sum(w), 1)
})

test_that("fp_sbc and fp_threshold_sweep run end to end", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 40, nocc = 4, seed = 4)
  fa <- list(chains = 1, iter_warmup = 150, iter_sampling = 150)
  sbc <- suppressWarnings(fp_sbc(s$data, occ = ~ x, fp = ~ 1,
                                 random_recorder = FALSE, reps = 3, ndraws = 19,
                                 fit_args = fa))
  expect_s3_class(sbc, "fp_sbc")
  expect_true(all(sbc$rank >= 0 & sbc$rank <= 19, na.rm = TRUE))
  expect_true("log_lik" %in% sbc$parameter)
  expect_s3_class(summary(sbc, bins = 4), "summary.fp_sbc")
  expect_error(summary(sbc, bins = 3), "must divide")
  expect_error(fp_sbc(s$data, fp = ~ 1, random_recorder = FALSE, ndraws = 500,
                      fit_args = fa), "exceeds")
  skip_if_not_installed("ggplot2")
  expect_s3_class(plot(sbc, bins = 4), "ggplot")

  a <- s$data$auto
  a$file <- paste(a$site, a$occasion)
  set.seed(1)
  det <- data.frame(file = a$file, confidence = ifelse(a$y == 1, runif(nrow(a), 0.4, 1),
                                                       runif(nrow(a), 0, 0.5)))
  sw <- suppressWarnings(fp_threshold_sweep(a[c("site", "file")], det, thresholds = c(0.3, 0.6),
                           sites = s$data$sites, verified = s$data$verified,
                           occ = ~ x, fp = ~ 1, chains = 1, iter_warmup = 150,
                           iter_sampling = 150, seed = 1, show_messages = FALSE))
  expect_s3_class(sw, "fp_threshold_sweep")
  expect_setequal(unique(sw$threshold), c(0.3, 0.6))
  expect_s3_class(plot(sw), "ggplot")
})
