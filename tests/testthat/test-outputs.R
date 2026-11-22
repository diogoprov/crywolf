finite_num <- function(df) {
  num <- vapply(df, is.numeric, logical(1))
  all(is.finite(as.matrix(df[num])))
}

#' @srrstats {G5.3} No missing or infinite values in the results.
test_that("the results of a fit hold no missing or infinite values", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 120, nocc = 8, seed = 31)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                                chains = 2, iter_warmup = 300, iter_sampling = 300,
                                threads_per_chain = 1, seed = 3))
  expect_true(all(is.finite(coef(m))))
  expect_true(finite_num(summary(m)$table))
  expect_true(finite_num(predict(m, newdata = data.frame(x = seq(-2, 2, 0.5)))))
  expect_true(finite_num(fp_sites(m, ndraws = 200, seed = 1)))
  pp <- fp_ppc(m, ndraws = 200, seed = 1)
  expect_true(all(is.finite(pp$yrep)) && all(is.finite(pp$fit_obs)) &&
                all(is.finite(pp$fit_rep)))
  expect_false(anyNA(unlist(m$diagnostics)))
  st <- m$stan_data
  expect_false(anyNA(unlist(st[vapply(st, is.numeric, logical(1))])))
})

#' @srrstats {BS2.7, BS2.8} Starting values, and a run started from a previous fit.
test_that("a run can start from the end of a previous one", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 80, nocc = 6, seed = 32)
  fit <- function(...) suppressWarnings(
    fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE, chains = 2,
            iter_sampling = 200, threads_per_chain = 1, seed = 4, ...)
  )
  m1 <- fit(iter_warmup = 200)
  m2 <- fit(iter_warmup = 100, init = m1$fit)
  expect_s3_class(m2, "fp_occu")
  expect_equal(posterior::ndraws(m2$draws), 400)
  # starting values from a list per chain
  m3 <- fit(iter_warmup = 100, init = function() list(beta = c(0, 0)))
  expect_s3_class(m3, "fp_occu")
})

#' @srrstats {BS2.13} Progress silenced, warnings kept.
test_that("progress output can be silenced while sampler warnings remain", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 33)
  run <- function(...) fp_occu(s$data, occ = ~ x, random_recorder = FALSE,
                               chains = 1, threads_per_chain = 1, seed = 5,
                               refresh = 0, show_messages = FALSE,
                               show_exceptions = FALSE, ...)
  out <- utils::capture.output(m <- suppressWarnings(run(iter_warmup = 200,
                                                         iter_sampling = 200)))
  expect_length(out, 0)
  # a run far too short to converge still warns, with the sampler class
  w <- NULL
  utils::capture.output(invisible(withCallingHandlers(
    run(iter_warmup = 10, iter_sampling = 20),
    crywolf_sampler_warning = function(c) {
      w <<- c
      invokeRestart("muffleWarning")
    },
    warning = function(c) invokeRestart("muffleWarning")
  )))
  expect_s3_class(w, "crywolf_sampler_warning")
})

#' @srrstats {BS6.2} Warm-up draws in trace plots.
test_that("trace plots show the warm-up when the fit kept it", {
  skip_if_no_cmdstan()
  skip_if_not_installed("bayesplot")
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 34)
  fit <- function(...) suppressWarnings(
    fp_occu(s$data, occ = ~ x, random_recorder = FALSE, chains = 2,
            iter_warmup = 150, iter_sampling = 100, threads_per_chain = 1,
            seed = 6, refresh = 0, show_messages = FALSE, ...)
  )
  m <- fit(save_warmup = TRUE)
  w <- crywolf:::.warmup_draws(m)
  expect_equal(posterior::niterations(w$draws), 250)
  expect_equal(w$n_warmup, 150)
  expect_setequal(posterior::variables(w$draws), posterior::variables(m$draws))
  expect_equal(posterior::niterations(m$draws), 100)
  expect_s3_class(plot(m, type = "trace"), "ggplot")
  expect_null(crywolf:::.warmup_draws(fit()))
})
