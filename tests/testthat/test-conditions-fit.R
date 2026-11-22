# Conditions (errors and warnings) of the fitting functions and their data
# builders. Errors are raised before sampling; the warnings that fp_occu(),
# fp_occu_dyn() and fp_occu_ms() emit before compiling are reached by mocking
# the compiler so that the call stops right after them, with no CmdStan.

mock_stop <- function(...) stop("compiler reached", call. = FALSE)

#' @srrstats {G5.2b} Conditions raised by fp_occu(), fp_occu_dyn(), fp_occu_ms() and their data builders.
test_that("fp_occu() checks the sampler settings", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  expect_error(fp_occu(s$data, chains = 0),
               "`chains` must be a single whole number of at least 1.", fixed = TRUE)
  expect_error(fp_occu(s$data, chains = 1.5), "`chains` must be a single whole number",
               fixed = TRUE)
  expect_error(fp_occu(s$data, seed = -1),
               "`seed` must be NULL or a single non-negative whole number.", fixed = TRUE)
  expect_error(fp_occu(s$data, seed = c(1, 2)), "`seed` must be NULL", fixed = TRUE)
})

test_that("fp_occu() warns before compiling: no verified data, two site effects", {
  local_mocked_bindings(fp_compile = mock_stop)
  s0 <- simulate_fp(nsite = 20, nocc = 3, frac_verified = 0, seed = 2)
  w <- NULL
  expect_error(
    withCallingHandlers(fp_occu(s0$data, threads_per_chain = 1),
                        warning = function(x) {
                          w <<- c(w, conditionMessage(x))
                          invokeRestart("muffleWarning")
                        }),
    "compiler reached")
  expect_match(w, "No verified data: the false-positive rate is identified only by the priors",
               fixed = TRUE, all = FALSE)

  s <- simulate_fp(nsite = 20, nocc = 3, seed = 2)
  expect_error(
    expect_warning(fp_occu(s$data, random_fp = TRUE, random_det = TRUE,
                           threads_per_chain = 1),
                   "Site effects on both the false-positive rate and true detection mix slowly",
                   fixed = TRUE),
    "compiler reached")
  # neither warning with verified data and one site effect
  expect_error(expect_no_warning(fp_occu(s$data, random_fp = TRUE, threads_per_chain = 1)),
               "compiler reached")
})

test_that(".design() refuses malformed formulas", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  build <- function(...) crywolf:::.build_stan_data(s$data, ..., FALSE, fp_priors())
  expect_error(build("x", ~ 1, ~ 1), "`occ` must be a one-sided formula, e.g. ~ x.",
               fixed = TRUE)
  expect_error(build(y ~ x, ~ 1, ~ 1), "`occ` must be a one-sided formula", fixed = TRUE)
  expect_error(build(~ 1, ~ 1, "aci"), "`fp` must be a one-sided formula", fixed = TRUE)
  expect_error(build(~ x + nope + zip, ~ 1, ~ 1),
               "Variable(s) in `occ` not found: nope, zip.", fixed = TRUE)
  expect_error(build(~ 1, ~ nope, ~ 1), "Variable(s) in `det` not found: nope.",
               fixed = TRUE)
})

test_that(".build_stan_data() requires intercepts and a false-positive process", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  build <- function(d, occ = ~ 1, det = ~ 1, fp = ~ 1, ...) {
    crywolf:::.build_stan_data(d, occ, det, fp, FALSE, fp_priors(), ...)
  }
  expect_error(build(s$data, occ = ~ x - 1), "`occ` must keep the intercept.", fixed = TRUE)
  expect_error(build(s$data, det = ~ 0 + bio), "`det` must keep the intercept.", fixed = TRUE)
  expect_error(build(s$data, fp = ~ 0 + aci),
               "`fp` must keep the intercept; it is tied to the detection intercept",
               fixed = TRUE)
  expect_error(build(s$data, ver = ~ 0 + x), "`ver` must keep the intercept.", fixed = TRUE)

  s0 <- simulate_fp(nsite = 20, nocc = 3, frac_verified = 0, seed = 1)
  expect_error(build(s0$data, ver = ~ x),
               "`ver` has covariates but the data have no verified occasions.", fixed = TRUE)
  expect_equal(build(s0$data)$stan_data$Kfp, 0)

  sc <- simulate_fp(nsite = 20, nocc = 3, frac_confirmed = 0.5, seed = 1)
  expect_error(build(sc$data, fp = NULL),
               "The observation-confirmation design needs a false-positive process",
               fixed = TRUE)
  expect_equal(build(sc$data)$stan_data$estimate_fp, 1L)
})

test_that("the markov column is checked", {
  s <- simulate_fp(nsite = 10, nocc = 3, seed = 1)
  build <- function(d) {
    crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(), markov = "occasion")
  }
  d <- s$data
  d$auto$occasion <- as.character(d$auto$occasion)
  expect_error(build(d), "The `markov` column 'occasion' must be numeric, a date, or a date-time.",
               fixed = TRUE)
  d <- s$data
  d$auto$occasion[2] <- NA
  expect_error(build(d), "The `markov` column 'occasion' has missing values.", fixed = TRUE)
  d <- s$data
  d$auto$occasion[d$auto$site == 1] <- 1
  expect_error(build(d), "Two occasions of the same site share the same `markov` value",
               fixed = TRUE)
  # valid columns pass, dates included
  expect_equal(crywolf:::.markov_adjacent(c(1, 1, 1, 2), as.Date("2024-01-01") + c(0, 1, 3, 0),
                                          "t"),
               c(FALSE, TRUE, FALSE, FALSE))
})

test_that("fp_occu_dyn() checks its data and warns before compiling", {
  s <- simulate_fp(nsite = 20, nocc = 2, seed = 1)
  expect_error(fp_occu_dyn(s$data),
               "`data` must be an `fp_data` object with a season column; see `fp_data(season = )`.",
               fixed = TRUE)
  expect_error(fp_occu_dyn(list()), "with a season column", fixed = TRUE)

  local_mocked_bindings(.compile_model = mock_stop)
  dy <- simulate_fp_dyn(nsite = 20, nseason = 2, nocc = 2, seed = 1)
  d0 <- dy$data
  d0$verified <- NULL
  expect_error(
    expect_warning(fp_occu_dyn(d0, threads_per_chain = 1),
                   "No verified data: the false-positive rate of the dynamic model is identified only by the priors",
                   fixed = TRUE),
    "compiler reached")
  expect_error(
    expect_warning(fp_occu_dyn(dy$data, random_fp = TRUE, random_det = TRUE,
                               threads_per_chain = 1),
                   "of the dynamic model mix slowly", fixed = TRUE),
    "compiler reached")
  expect_error(expect_no_warning(fp_occu_dyn(dy$data, threads_per_chain = 1)),
               "compiler reached")
})

test_that(".build_dyn_stan_data() checks the transition formulas", {
  build <- function(d, psi1 = ~ 1, col = ~ 1, ext = ~ 1) {
    crywolf:::.build_dyn_stan_data(d, psi1, col, ext, ~ 1, ~ 1, FALSE, fp_priors())
  }
  d1 <- simulate_fp_dyn(nsite = 20, nseason = 1, nocc = 3, seed = 1)$data
  expect_error(build(d1, col = ~ x), "`col` has covariates but the data have one season.",
               fixed = TRUE)
  expect_error(build(d1, ext = ~ x), "`ext` has covariates but the data have one season.",
               fixed = TRUE)
  expect_equal(build(d1)$stan_data$nseason, 1)

  d2 <- simulate_fp_dyn(nsite = 20, nseason = 2, nocc = 2, seed = 1)$data
  msg <- "`psi1`, `col` and `ext` must keep the intercept."
  expect_error(build(d2, psi1 = ~ 0 + x), msg, fixed = TRUE)
  expect_error(build(d2, col = ~ x - 1), msg, fixed = TRUE)
  expect_error(build(d2, ext = ~ 0 + x), msg, fixed = TRUE)
})

test_that(".fp_sites_dyn() refuses site columns it would overwrite", {
  d <- simulate_fp_dyn(nsite = 10, nseason = 2, nocc = 2, seed = 1)$data
  d$sites$psi <- 1
  d$sites$p10 <- 2
  sd <- crywolf:::.build_dyn_stan_data(d, ~ 1, ~ 1, ~ 1, ~ 1, ~ 1, FALSE,
                                       fp_priors())$stan_data
  one <- matrix(0, 2, 1)
  local_mocked_bindings(.dyn_pars = function(...) list(beta = one, beta_col = one,
                                                       beta_ext = one),
                        .dyn_ffbs = function(...) NULL)
  obj <- list(stan_data = sd, fit = NULL, data = d)
  expect_error(crywolf:::.fp_sites_dyn(obj, 0.9, NULL, NULL),
               "The sites table already has column(s) psi, p10; rename them before calling fp_sites() on a dynamic fit.",
               fixed = TRUE)
})

test_that("fp_occu_ms() refuses unavailable options and checks fp_group", {
  ms <- simulate_fp_ms(nspecies = 2, nsite = 20, nocc = 3, seed = 1)$data
  quiet <- function(x) suppressMessages(x)
  expect_error(quiet(fp_occu_ms(ms, markov = "occasion")),
               "`markov` is not available in fp_occu_ms().", fixed = TRUE)
  expect_error(quiet(fp_occu_ms(ms, markov = "occasion", random_det = TRUE)),
               "`markov` and `random_det` are not available in fp_occu_ms().", fixed = TRUE)

  ms$sites$g <- rep(1:2, length.out = nrow(ms$sites))
  ms$sites$g[3] <- NA
  expect_error(quiet(fp_occu_ms(ms, fp_group = "g")),
               "The column `g` of the sites table has missing values.", fixed = TRUE)
})

test_that("fp_occu_ms() warns without verified data before compiling", {
  local_mocked_bindings(.compile_model = mock_stop)
  ms <- simulate_fp_ms(nspecies = 2, nsite = 20, nocc = 3, seed = 1)$data
  ms$verified <- NULL
  expect_error(
    expect_warning(suppressMessages(fp_occu_ms(ms, threads_per_chain = 1)),
                   "No verified data: the false-positive rates are identified only by the priors",
                   fixed = TRUE),
    "compiler reached")
})
