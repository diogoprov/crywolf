known_generator <- function(s, fp = ~ 1) {
  sl <- crywolf:::.build_stan_data(s$data, ~ x, ~ 1, fp, FALSE, fp_priors())
  list(data = s$data, formulas = list(occ = ~ x, det = ~ 1, fp = fp),
       stan_data = sl$stan_data, design = sl$design, occasions = sl$occasions,
       draws = list(beta = rbind(c(0, -0.8), c(0.2, -0.6)),
                    alpha = matrix(qlogis(0.6), 2, 1),
                    delta0 = c(qlogis(0.6), 0), u = NULL,
                    gamma0 = c(qlogis(0.05), qlogis(0.1)),
                    gamma_fp = matrix(0, 2, 0)))
}

test_that("simulating from known values reproduces the detection rates", {
  s <- simulate_fp(nsite = 3000, nocc = 10, frac_verified = 0.3, seed = 1)
  gen <- known_generator(s)
  set.seed(2)
  p <- crywolf:::.pick_pars(gen$draws, "mean")
  expect_equal(p$gamma0, mean(gen$draws$gamma0))
  r <- crywolf:::.simulate_from(gen, list(beta = c(0, -0.8), alpha = qlogis(0.6),
                                          delta0 = qlogis(0.6), u = NULL,
                                          gamma0 = qlogis(0.05),
                                          gamma_fp = numeric(0)))
  a <- r$data$auto
  z <- r$truth$z
  expect_equal(mean(a$y[z[a$site] == 1]), 0.6, tolerance = 0.03)
  expect_equal(mean(a$y[z[a$site] == 0]), 0.05, tolerance = 0.1)
  v <- r$data$verified
  expect_equal(sum(v$y[z[v$site] == 0]), 0)
  expect_equal(dim(a), dim(s$data$auto))
  expect_equal(nrow(v), nrow(s$data$verified))
  expect_equal(r$truth$p10_bar, 0.05)
})

test_that("posterior draws keep clean names and targets fill absent slopes", {
  s <- simulate_fp(nsite = 50, nocc = 4, seed = 1)
  gen <- known_generator(s)
  set.seed(1)
  r <- crywolf:::.simulate_from(gen, crywolf:::.pick_pars(gen$draws, "posterior"))
  expect_setequal(names(r$truth$pars),
                  c("occ_(Intercept)", "occ_x", "det_(Intercept)",
                    "verified_(Intercept)", "fp_(Intercept)"))
  tg <- crywolf:::.fit_targets(r$truth, ~ x, gen)
  expect_equal(tg[["fp_x"]], 0)
  expect_true("p10_bar" %in% names(tg))
  expect_false("p10_bar" %in% names(crywolf:::.fit_targets(r$truth, NULL, gen)))
})

test_that("site effects on false positives are drawn anew from posterior means", {
  s <- simulate_fp(nsite = 4000, nocc = 4, seed = 1)
  gen <- known_generator(s)
  gen$draws$sigma_fp <- c(1, 1)
  gen$draws$e_fp <- matrix(0, 2, 4000)        # shrunk means would be zero
  set.seed(3)
  p <- crywolf:::.pick_pars(gen$draws, "mean")
  expect_null(p$e_fp)
  r <- crywolf:::.simulate_from(gen, p)
  # the new site effects have the posterior standard deviation (1)
  lp <- qlogis(r$truth$p10_site) - mean(gen$draws$gamma0)
  expect_equal(sd(lp), 1, tolerance = 0.05)
  set.seed(3)
  pp <- crywolf:::.pick_pars(gen$draws, "posterior")
  expect_length(pp$e_fp, 4000)
})

test_that("simulate() and fp_calibrate_fit() work on fits of both designs", {
  skip_if_no_cmdstan()
  fa <- list(chains = 1, iter_warmup = 150, iter_sampling = 150)
  s <- simulate_fp(nsite = 80, nocc = 6, seed = 5)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1,
                                random_recorder = FALSE, chains = 1,
                                iter_warmup = 200, iter_sampling = 200,
                                seed = 1, show_messages = FALSE))
  new <- simulate(m, nsim = 2, seed = 1)
  expect_length(new, 2)
  expect_s3_class(new[[1]]$data, "fp_data")
  cal <- suppressWarnings(fp_calibrate_fit(
    m, models = list(constant = ~ 1, fp_x = ~ x), reps = 2, fit_args = fa))
  expect_s3_class(cal, "fp_calibration")
  expect_setequal(unique(cal$model), c("constant", "fp_x"))
  expect_true("fp_x" %in% cal$parameter)
  expect_output(print(summary(cal)), "Occupancy slope")

  s2 <- simulate_fp(nsite = 80, nocc = 6, frac_confirmed = 0.4,
                    frac_verified = 0, seed = 6)
  m2 <- suppressWarnings(fp_occu(s2$data, occ = ~ x, fp = ~ 1,
                                 random_recorder = FALSE, chains = 1,
                                 iter_warmup = 200, iter_sampling = 200,
                                 seed = 1, show_messages = FALSE))
  d2 <- simulate(m2, seed = 2)[[1]]$data
  conf <- d2$auto[!is.na(d2$auto$n_true), ]
  expect_equal(nrow(conf), sum(!is.na(s2$data$auto$n_true)))
  expect_equal(conf$y, pmax(conf$n_true, conf$n_false))
})
