test_that("fp_occu_mc() checks its input", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  expect_error(fp_occu_mc(s$data), "list of `fp_data`")
  expect_error(fp_occu_mc(list(s$data, 1)), "list of `fp_data`")
  other <- simulate_fp(nsite = 21, nocc = 3, seed = 2)
  expect_error(fp_occu_mc(list(s$data, other$data)), "same sites")
  moved <- s$data
  moved$sites$x <- moved$sites$x + 1
  expect_error(fp_occu_mc(list(s$data, moved)), "same sites table")
})

test_that("one version reproduces fp_occu(), and the variance splits", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 1)
  args <- list(occ = ~ x, fp = ~ 1, random_recorder = FALSE, chains = 2,
               parallel_chains = 1, iter_warmup = 150, iter_sampling = 150,
               threads_per_chain = 1, refresh = 0, show_messages = FALSE)
  m1 <- suppressWarnings(do.call(fp_occu, c(list(s$data, seed = 5), args)))
  mc1 <- suppressWarnings(do.call(fp_occu_mc, c(list(list(s$data), seed = 5), args)))
  expect_identical(unclass(posterior::as_draws_matrix(mc1$draws)),
                   unclass(posterior::as_draws_matrix(m1$draws)))
  expect_equal(summary(mc1)$table$share_between, rep(0, nrow(summary(mc1)$table)))

  v2 <- s$data
  set.seed(2)
  hit <- v2$auto$y == 1
  v2$auto$y[hit] <- stats::rbinom(sum(hit), 1, 0.7)
  mc2 <- suppressWarnings(do.call(fp_occu_mc, c(list(list(s$data, v2), seed = 5), args)))
  expect_equal(posterior::nchains(mc2$draws), 4)
  expect_equal(nrow(mc2$diagnostics), 2)
  t <- summary(mc2)$table
  # law of total variance with equal shares: pooled = within + between
  expect_equal(t$sd^2, t$sd_within^2 + t$sd_between^2, tolerance = 0.02)
  expect_true(all(t$share_between >= 0 & t$share_between <= 1))
  expect_output(print(mc2), "pooled over 2 fits")
  expect_identical(as_draws_df(mc2), mc2$draws)
})

test_that("identical versions put almost nothing between the fits", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 5, seed = 1)
  mc <- suppressWarnings(fp_occu_mc(list(s$data, s$data, s$data), occ = ~ x, fp = ~ 1,
                                    random_recorder = FALSE, chains = 2,
                                    parallel_chains = 1, iter_warmup = 150,
                                    iter_sampling = 150, threads_per_chain = 1,
                                    seed = 7, refresh = 0, show_messages = FALSE))
  t <- summary(mc)$table
  # only Monte Carlo error separates the fits of one dataset
  expect_lt(max(t$share_between), 0.1)
})

test_that("arguments set to NULL reach fp_occu(), also in parallel", {
  skip_if_no_cmdstan()
  skip_if_not_installed("furrr")
  s <- simulate_fp(nsite = 40, nocc = 4, seed = 1)
  old <- future::plan(future::sequential)
  on.exit(future::plan(old))
  mc <- suppressWarnings(fp_occu_mc(list(s$data), occ = ~ x, fp = NULL,
                                    random_recorder = FALSE, chains = 1,
                                    iter_warmup = 50, iter_sampling = 50, seed = 1,
                                    refresh = 0, show_messages = FALSE, parallel = TRUE,
                                    max_rhat = Inf))
  expect_false(any(grepl("^fp_", posterior::variables(mc$draws))))
})

test_that("fits that did not converge are left out, with a warning", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 40, nocc = 4, seed = 1)
  v2 <- s$data
  v2$auto$y[v2$auto$y == 1][1:3] <- 0L
  run <- function(max_rhat) {
    fp_occu_mc(list(s$data, v2), occ = ~ x, fp = ~ 1, random_recorder = FALSE,
               chains = 2, iter_warmup = 60, iter_sampling = 60,
               threads_per_chain = 1, seed = 3, refresh = 0,
               show_messages = FALSE, max_rhat = max_rhat)
  }
  all <- suppressWarnings(run(Inf))
  expect_true(all(all$diagnostics$kept))
  r <- all$diagnostics$max_rhat
  skip_if(r[1] == r[2], "both fits have the same R-hat")
  # identical fits again (same seeds, one thread): keep only the better one
  w <- character()
  one <- withCallingHandlers(run(mean(r)), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    invokeRestart("muffleWarning")
  })
  expect_match(w, "did not converge", all = FALSE)
  expect_equal(one$n_fits, 1)
  expect_equal(one$diagnostics$kept, r < mean(r))
  expect_output(print(one), "left out for lack of convergence")
  expect_equal(names(one$summaries), as.character(which(r < mean(r))))
  expect_error(suppressWarnings(run(min(r) - 1e-6)), "No fit converged")
  expect_error(run(NA), "single number")
})
