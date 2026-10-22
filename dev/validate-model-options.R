# Validation of the model options added in crywolf 1.5.0 (verified-method
# covariates, occasion-level false-positive covariates, Markov term).
# 1. With the options off, the draws are identical to those of the 1.4.0 Stan
#    program (same seed), in three designs.
# 2. With the options on, the log-likelihood computed in R (.ppc_core) equals
#    the Stan log_lik, in six combinations.
# 3. In simulations, the model recovers the effect of each option.
# Since crywolf 1.7.0, part 2 runs as tests/testthat/test-likelihood.R and
# the identity of draws across versions is guarded by
# tests/testthat/test-regression.R; this script keeps part 3 (recovery,
# too slow for the tests) and the record of the 1.5.0 validation.
# Run from the package root: Rscript dev/validate-model-options.R
# Output: dev/validate-model-options.txt
sink(file.path("dev", "validate-model-options.txt"), split = TRUE)
devtools::load_all(quiet = TRUE)
base_stan <- tempfile(fileext = ".stan")
writeLines(system2("git", c("show", "v1.4.0:inst/stan/fp_occu.stan"), stdout = TRUE), base_stan)
opts <- list(cpp_options = list(stan_threads = TRUE), quiet = TRUE)
old <- do.call(cmdstanr::cmdstan_model, c(list(base_stan), opts))
new <- fp_compile()
# 1. options off: draws identical to 1.4.0
strip <- function(sd) { sd[c("Kver_s", "X1s", "fp_by_row", "Xfp2", "Xfp3")] <- NULL; sd }
same <- function(s, ...) {
  m <- suppressWarnings(fp_occu(s$data, ..., seed = 7, chains = 2, iter_warmup = 300, iter_sampling = 300, refresh = 0, show_messages = FALSE))
  f0 <- old$sample(data = strip(m$stan_data), seed = 7, chains = 2, parallel_chains = 2, threads_per_chain = 1, iter_warmup = 300, iter_sampling = 300, adapt_delta = 0.95, refresh = 0, show_messages = FALSE)
  v <- c("beta", "alpha", "fp_gap", "delta0", "log_lik")
  max(abs(unclass(posterior::as_draws_matrix(m$fit$draws(v))) - unclass(posterior::as_draws_matrix(f0$draws(v)))))
}
s1 <- simulate_fp(nsite = 150, nocc = 8, a_bio = 0.5, sd_recorder = 0.5, conf_effect = 1, conf_cor = 0.5, seed = 3)
s2 <- simulate_fp(nsite = 120, nocc = 6, frac_confirmed = 0.4, frac_verified = 0.2, seed = 6)
s3 <- simulate_fp(nsite = 120, nocc = 6, seed = 8)
cat("OFF Miller RE fp~w:", same(s1, occ = ~ x, det = ~ bio, fp = ~ w, random_recorder = TRUE), "\n")
cat("OFF Chambert:", same(s2, occ = ~ x, fp = ~ 1, random_recorder = FALSE), "\n")
cat("OFF fp=NULL:", same(s3, occ = ~ x, fp = NULL, random_recorder = FALSE), "\n")
# 2. new options: R likelihood == Stan log_lik
llcheck <- function(label, d, ...) {
  m <- suppressWarnings(fp_occu(d, ..., seed = 3, chains = 1, iter_warmup = 200, iter_sampling = 100, refresh = 0, show_messages = FALSE))
  pars <- .ppc_pars(m$fit, m$stan_data, NULL)
  r <- .ppc_core(m$stan_data, pars, "automated", "freeman-tukey")
  ll <- unclass(posterior::as_draws_matrix(m$fit$draws("log_lik")))
  cat(sprintf("LL %-45s max|R - Stan| = %.2e  divergences %d\n", label, max(abs(r$ll - ll)), m$diagnostics$divergences))
  invisible(m)
}
set.seed(1)
d1 <- s1$data; d1$verified$listener <- rnorm(nrow(d1$verified))
llcheck("verified ~ listener", d1, occ = ~ x, det = ~ bio, fp = ~ 1, ver = ~ listener, random_recorder = TRUE)
d2 <- s1$data; d2$auto$noise <- rnorm(nrow(d2$auto))
llcheck("fp ~ noise (occasion) + x (site), Miller", d2, occ = ~ x, det = ~ bio, fp = ~ noise + x, random_recorder = TRUE)
d3 <- s2$data; d3$auto$noise <- rnorm(nrow(d3$auto))
llcheck("fp ~ noise, Chambert", d3, occ = ~ x, fp = ~ noise, random_recorder = FALSE)
llcheck("markov, det ~ bio (GLM)", s1$data, occ = ~ x, det = ~ bio, fp = ~ 1, markov = "occasion", random_recorder = TRUE)
llcheck("markov, det ~ 1, random_fp", s3$data, occ = ~ x, fp = ~ 1, markov = "occasion", random_recorder = FALSE, random_fp = TRUE)
llcheck("markov + fp ~ noise + ver, Chambert", { d3$verified$listener <- rnorm(nrow(d3$verified)); d3 }, occ = ~ x, fp = ~ noise, ver = ~ listener, markov = "occasion", random_recorder = FALSE)

# 3. recovery
set.seed(11)
s <- simulate_fp(nsite = 400, nocc = 20, frac_verified = 0.4, seed = 11)
d <- s$data
d$auto$noise <- rnorm(nrow(d$auto)); d$verified$listener <- rnorm(nrow(d$verified))
run <- function(label, truth, fp, ver, markov) {
  sl <- .build_stan_data(d, ~ x, ~ 1, fp, FALSE, fp_priors(), ver = ver, markov = markov)
  gen <- list(data = d, formulas = list(occ = ~ x, det = ~ 1, fp = fp, ver = ver),
              stan_data = sl$stan_data, design = sl$design, occasions = sl$occasions)
  p <- list(beta = c(0, -0.8), alpha = truth$alpha, delta0 = 0.4,
            delta_s = truth$delta_s, u = NULL, gamma0 = truth$alpha[1] - 2.5,
            gamma_fp = truth$gamma_fp)
  set.seed(5); sim <- .simulate_from(gen, p)
  m <- suppressWarnings(fp_occu(sim$data, occ = ~ x, fp = fp, ver = ver, markov = markov,
                                random_recorder = FALSE, seed = 1, refresh = 0, show_messages = FALSE))
  t <- m$draws; tab <- summary(m)$table
  for (v in names(truth$check)) {
    r <- tab[tab$variable == v, ]
    cat(sprintf("%-28s %-16s true %6.2f  est %6.2f [%6.2f, %6.2f]  rhat %.3f div %d\n",
                label, v, truth$check[[v]], r$mean, r$lower, r$upper, m$diagnostics$max_rhat, m$diagnostics$divergences))
  }
}
run("verified covariate", list(alpha = qlogis(0.4), delta_s = 1, gamma_fp = numeric(0),
    check = c(verified_listener = 1, occ_x = -0.8)), fp = ~ 1, ver = ~ listener, markov = NULL)
run("fp occasion covariate", list(alpha = qlogis(0.4), delta_s = numeric(0), gamma_fp = 1,
    check = c(fp_noise = 1, occ_x = -0.8)), fp = ~ noise, ver = ~ 1, markov = NULL)
run("markov runs", list(alpha = c(qlogis(0.2), 2), delta_s = numeric(0), gamma_fp = numeric(0),
    check = c(det_y_prev = 2, `det_(Intercept)` = qlogis(0.2), occ_x = -0.8)), fp = ~ 1, ver = ~ 1, markov = "occasion")
# the Markov term with missing occasions: runs are not carried across gaps
set.seed(12)
s_gap <- simulate_fp(nsite = 400, nocc = 20, frac_verified = 0.4, seed = 12)
d <- s_gap$data
d$auto$noise <- rnorm(nrow(d$auto)); d$verified$listener <- rnorm(nrow(d$verified))
d$auto <- d$auto[!(d$auto$occasion %in% c(5, 12)), ]
run("markov runs, gaps at 5 and 12", list(alpha = c(qlogis(0.2), 2), delta_s = numeric(0),
    gamma_fp = numeric(0), check = c(det_y_prev = 2, `det_(Intercept)` = qlogis(0.2),
    occ_x = -0.8)), fp = ~ 1, ver = ~ 1, markov = "occasion")
sink()
