sim_data <- function(nsite = 40) {
  suppressWarnings(simulate_fp(nsite = nsite, nocc = 5, seed = 3)$data)
}
build <- function(d, occ = ~ x, fp = ~ 1) {
  crywolf:::.build_stan_data(d, occ, ~ 1, fp, FALSE, fp_priors())
}

test_that("list columns in a formula are refused with a clear message", {
  d <- sim_data()
  d$sites$l <- as.list(d$sites$x)
  expect_error(build(d, occ = ~ l), "`occ` are list columns: l")
})

test_that("collinear, constant and aliased covariates are flagged", {
  d <- sim_data()
  d$sites$x2 <- 2 * d$sites$x
  expect_warning(build(d, occ = ~ x + x2), "collinear: x2",
                 class = "crywolf_design_warning")
  d$sites$k <- 1
  expect_warning(build(d, fp = ~ k), "`fp` are collinear: k",
                 class = "crywolf_design_warning")
  d$sites$f <- factor(rep("a", nrow(d$sites)), levels = c("a", "b"))
  expect_warning(build(d, occ = ~ f), "droplevels",
                 class = "crywolf_design_warning")
  # more covariates than sites
  small <- sim_data(nsite = 4)
  for (k in 1:5) small$sites[[paste0("z", k)]] <- stats::rnorm(4)
  expect_warning(build(small, occ = ~ z1 + z2 + z3 + z4 + z5), "collinear",
                 class = "crywolf_design_warning")
  expect_silent(build(d, occ = ~ x))
})

test_that("ordered factors get polynomial contrasts, unordered ones treatment", {
  d <- sim_data()
  d$sites$h <- factor(rep(c("lo", "mid", "hi"), length.out = nrow(d$sites)),
                      levels = c("lo", "mid", "hi"))
  expect_true("hmid" %in% colnames(crywolf:::.design(~ h, d$sites, "occ")$X))
  d$sites$h <- factor(d$sites$h, ordered = TRUE)
  expect_true(all(c("h.L", "h.Q") %in% colnames(crywolf:::.design(~ h, d$sites, "occ")$X)))
})

test_that("sampler settings are checked before compiling", {
  d <- sim_data()
  run <- function(...) suppressWarnings(fp_occu(d, random_recorder = FALSE, ...))
  expect_error(run(chains = 0), "`chains` must be")
  expect_error(run(chains = 2.5), "`chains` must be")
  expect_error(run(iter_warmup = -1), "`iter_warmup` must be")
  expect_error(run(iter_sampling = 0), "`iter_sampling` must be")
  expect_error(run(adapt_delta = 1), "`adapt_delta` must be")
  expect_error(run(adapt_delta = c(0.8, 0.9)), "`adapt_delta` must be")
  expect_error(run(seed = "a"), "`seed` must be")
})

test_that("sampler warnings have their own class and cover tree depth and E-BFMI", {
  diag <- list(divergences = 0L, max_treedepth = 3L, low_ebfmi = 1L,
               max_rhat = 1, min_ess_bulk = 1000, min_ess_tail = 1000)
  w <- tryCatch(crywolf:::.warn_diagnostics(diag), warning = function(w) w)
  expect_s3_class(w, "crywolf_sampler_warning")
  expect_match(conditionMessage(w), "maximum tree depth")
  expect_match(conditionMessage(w), "E-BFMI")
  ok <- list(divergences = 0L, max_treedepth = 0L, low_ebfmi = 0L,
             max_rhat = 1, min_ess_bulk = 1000, min_ess_tail = 1000)
  expect_silent(crywolf:::.warn_diagnostics(ok))
  # muffled selectively
  expect_silent(withCallingHandlers(
    crywolf:::.warn_diagnostics(diag),
    crywolf_sampler_warning = function(w) invokeRestart("muffleWarning")
  ))
})
