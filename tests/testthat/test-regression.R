# Posterior regression test. References: dev/make-regression-refs.R.
# 1. The likelihood: log_lik recomputed from saved parameter draws must equal
#    the saved values. Deterministic, so any change to the likelihood shows.
# 2. The posterior: means from a fit with fixed seed and one thread must agree
#    with the saved ones within Monte Carlo error. Another platform gives
#    other draws, hence the tolerance; a change to priors, likelihood, or
#    sampling setup beyond about 8 Monte Carlo standard errors shows.
#    dev/check-regression-sensitivity.R shows which prior changes it flags.

refs_file <- test_path("fixtures", "regression-refs.rds")

test_that("the likelihood matches the reference values", {
  skip_if_no_cmdstan()
  refs <- readRDS(refs_file)
  designs <- .regression_designs()
  purrr::iwalk(refs, function(ref, name) {
    ll <- .regression_log_lik(.regression_stan_data(designs[[name]]), ref$params)
    expect_equal(ll, ref$log_lik, tolerance = 1e-8, ignore_attr = TRUE,
                 label = paste(name, "log_lik"))
  })
})

test_that("posterior means match the reference values", {
  skip_if_no_cmdstan()
  refs <- readRDS(refs_file)
  designs <- .regression_designs()
  purrr::iwalk(refs, function(ref, name) {
    fit <- .regression_fit(designs[[name]])
    # a run that did not converge is not evidence that the posterior moved
    expect_lt(fit$diagnostics$max_rhat, 1.05,
              label = paste(name, "largest R-hat (the fit did not converge)"))
    now <- as.data.frame(posterior::summarise_draws(fit$draws, "mean", "mcse_mean"))
    both <- merge(ref$summary, now, by = "variable", suffixes = c("_ref", "_now"))
    expect_setequal(both$variable, ref$summary$variable)
    # the difference of two independent runs has standard error
    # sqrt(mcse_ref^2 + mcse_now^2); allow about 8.5 of them
    tol <- 6 * sqrt(2) * pmax(both$mcse_mean_ref, both$mcse_mean_now) + 1e-8
    off <- abs(both$mean_now - both$mean_ref) > tol
    expect_false(any(off), label = paste0(
      name, ": ", paste(both$variable[off], collapse = ", "), " moved"))
  })
})
