# fp_occu_rn() and simulate_fp_rn(): conditions without CmdStan (the
# compiler is mocked so that a call stops right after its checks), and, with
# CmdStan, the likelihood of the Stan program against one written in R
# occasion by occasion.

mock_stop_rn <- function(...) stop("compiler reached", call. = FALSE)

# log-likelihood of each site from the raw tables of an fp_data object, with
# N summed explicitly and every occasion a Bernoulli trial; r, p10 and rv are
# per site (vectors), lambda per site
rn_loglik_r <- function(data, lambda, r, p10, rv, n_max) {
  sites <- data$sites$site
  auto <- data$auto
  conf <- !is.na(auto$n_true)
  vapply(seq_along(sites), function(i) {
    a <- auto[auto$site == sites[i], , drop = FALSE]
    ca <- conf[auto$site == sites[i]]
    v <- data$verified[data$verified$site == sites[i], , drop = FALSE]
    ll_n <- vapply(0:n_max, function(n) {
      p_auto <- 1 - (1 - r[i])^n * (1 - p10[i])
      s1 <- 1 - (1 - r[i])^n
      ll <- sum(stats::dbinom(a$y[!ca], 1, p_auto, log = TRUE))
      ll <- ll + sum(log(ifelse(a$n_true[ca] > 0, s1, 1 - s1)) +
                       log(ifelse(a$n_false[ca] > 0, p10[i], 1 - p10[i])))
      ll + sum(stats::dbinom(v$y, 1, 1 - (1 - rv)^n, log = TRUE))
    }, numeric(1))
    lp <- stats::dpois(0:n_max, lambda[i], log = TRUE) + ll_n
    max(lp) + log(sum(exp(lp - max(lp))))
  }, numeric(1))
}

#' @srrstats {G5.2, G5.2a, G5.2b} Every condition of fp_occu_rn() and simulate_fp_rn() is triggered and its text checked.
test_that("fp_occu_rn() checks its arguments before compiling", {
  local_mocked_bindings(.compile_model = mock_stop_rn)
  s <- simulate_fp_rn(nsite = 20, nocc = 3, seed = 1)
  expect_error(fp_occu_rn(s$data, n_max = 0),
               "`n_max` must be a single positive whole number.", fixed = TRUE)
  expect_error(fp_occu_rn(s$data, n_max = 2.5),
               "`n_max` must be a single positive whole number.", fixed = TRUE)
  expect_error(fp_occu_rn(s$data, random_abund = NA),
               "`random_abund` must be TRUE or FALSE.", fixed = TRUE)
  s$data$auto$noise <- seq_len(nrow(s$data$auto))
  expect_error(fp_occu_rn(s$data, fp = ~ noise),
               "`fp` uses occasion-level covariates, which fp_occu_rn() does not support; use site covariates only.",
               fixed = TRUE)
  s0 <- simulate_fp_rn(nsite = 20, nocc = 3, frac_verified = 0, seed = 2)
  expect_warning(expect_error(fp_occu_rn(s0$data, threads_per_chain = 1), "compiler reached"),
                 "No verified data: the false-positive rate is identified only by the priors.",
                 fixed = TRUE)
  # without false positives no verified data are needed
  expect_error(expect_silent(fp_occu_rn(s0$data, fp = NULL, threads_per_chain = 1)),
               "compiler reached")
})

test_that("fp_occu_rn() passes the Royle-Nichols fields to Stan", {
  captured <- NULL
  local_mocked_bindings(.compile_model = function(...) {
    list(sample = function(data, ...) {
      captured <<- data
      stop("sampler reached", call. = FALSE)
    })
  })
  s <- simulate_fp_rn(nsite = 30, nocc = 4, frac_confirmed = 0.5, seed = 3)
  expect_error(fp_occu_rn(s$data, random_abund = TRUE, n_max = 15,
                          priors = fp_priors(abund_site_sd = 0.5),
                          threads_per_chain = 1), "sampler reached")
  expect_identical(captured$Nmax, 15L)
  expect_identical(captured$estimate_abund_re, 1L)
  expect_identical(captured$abund_site_sd, 0.5)
  expect_null(captured$Xfp2)
  # confirmed occasions collapse into rows whose weights add up to them all
  expect_identical(sum(captured$w3), sum(!is.na(s$data$auto$n_true)))
  expect_lt(captured$N3, sum(captured$w3))
  expect_false(is.unsorted(captured$site3))
})

test_that("simulate_fp_rn() checks its arguments and draws unbalanced designs", {
  expect_error(simulate_fp_rn(psi = c(0.2, 0.3)),
               "Argument(s) `psi` must be single finite numbers.", fixed = TRUE)
  expect_error(simulate_fp_rn(nsite = 5, nocc = c(2, 3)),
               "`nocc` must be one positive whole number or one per site.", fixed = TRUE)
  expect_error(simulate_fp_rn(nver = 2.5), "`nver` must be a whole number.", fixed = TRUE)
  expect_error(simulate_fp_rn(nsite = 5, nocc = 0),
               "`nocc` must be one positive whole number or one per site.", fixed = TRUE)
  nocc <- c(3, 5, 4, 6, 2)
  s <- simulate_fp_rn(nsite = 5, nocc = nocc, frac_confirmed = 1, seed = 4)
  expect_identical(as.integer(table(s$data$auto$site)), as.integer(nocc))
  expect_false(anyNA(s$data$auto$n_true))
  expect_identical(s$truth$z, as.integer(s$truth$N > 0))
  # a detection is a true or a false observation
  expect_identical(s$data$auto$y, pmax(s$data$auto$n_true, s$data$auto$n_false))
  expect_identical(simulate_fp_rn(nsite = 10, seed = 5)$data$auto,
                   simulate_fp_rn(nsite = 10, seed = 5)$data$auto)
})

test_that("fp_priors() takes abund_site_sd", {
  expect_identical(fp_priors(abund_site_sd = 2)$abund_site_sd, 2)
  expect_error(fp_priors(abund_site_sd = 0),
               "`abund_site_sd` must be a single positive number.", fixed = TRUE)
  expect_output(print(fp_priors()), "sd_abund_site   ~ half-Normal(0, 1) (fp_occu_rn())",
                fixed = TRUE)
})

#' @srrstats {G5.4, G5.4a} The marginal likelihood of fp_occu_rn.stan equals one written independently in R, occasion by occasion with the abundance summed explicitly.
test_that("the likelihood of fp_occu_rn() equals one written in R", {
  skip_if_no_cmdstan()
  s <- simulate_fp_rn(nsite = 25, nocc = c(rep(4, 12), rep(6, 13)), psi = 0.7,
                      sigma_r = 0.8, frac_verified = 0.5, frac_confirmed = 0.4,
                      seed = 6)
  m <- suppressWarnings(fp_occu_rn(s$data, occ = ~ x, random_det = TRUE,
                                   random_abund = TRUE, n_max = 25, chains = 1,
                                   iter_warmup = 150, iter_sampling = 20,
                                   threads_per_chain = 1, seed = 6, refresh = 0,
                                   show_messages = FALSE, sig_figs = 14))
  expect_s3_class(m, "fp_occu_rn")
  expect_output(suppressWarnings(print(m)),
                "Royle-Nichols occupancy model with false positives")
  expect_true(all(c("abund_x", "sd_det_site", "sd_abund_site", "psi_ref") %in%
                    posterior::variables(as_draws_df(m))))
  dm <- posterior::as_draws_matrix(m$fit$draws())
  d <- stats::setNames(as.numeric(dm[5, ]), colnames(dm))
  get <- function(v) unname(d[grepl(paste0("^", v, "\\["), names(d))])
  beta <- get("beta")
  lambda <- exp(beta[1] + beta[2] * s$data$sites$x + d[["sigma_lam[1]"]] * get("v_raw"))
  r <- stats::plogis(d[["alpha[1]"]] + d[["sigma_det[1]"]] * get("d_raw"))
  p10 <- rep(stats::plogis(d[["gamma0"]]), 25)
  rv <- stats::plogis(d[["delta0"]])
  ll_r <- rn_loglik_r(s$data, lambda, r, p10, rv, n_max = 25)
  # the collapsed rows are sums of Bernoulli terms, with no binomial constant
  expect_equal(get("log_lik"), ll_r, tolerance = 1e-8)
})

test_that("fp_occu_rn() warns when n_max cuts off the abundance", {
  skip_if_no_cmdstan()
  s <- simulate_fp_rn(nsite = 30, nocc = 5, psi = 0.99, seed = 7)
  expect_warning(
    withCallingHandlers(
      fp_occu_rn(s$data, n_max = 2, chains = 1, iter_warmup = 150, iter_sampling = 100,
                 threads_per_chain = 1, seed = 7, refresh = 0, show_messages = FALSE),
      crywolf_sampler_warning = function(w) invokeRestart("muffleWarning")),
    "probability that abundance exceeds `n_max` (2) is", fixed = TRUE)
})
