# Times one gradient evaluation of the log density, at a fixed point, for the
# previous model and the current one (500 calls through the model methods of
# 'cmdstanr', one thread), on three datasets: two simulated ones, with rows
# collapsed (det ~ 1) and with occasion-level covariates, and the curiango
# data. Deterministic, unlike sampler timings, which vary with the load of
# the machine. Usage, from the package root:
#   Rscript dev/profile/benchmark-glm-gradient.R <previous fp_occu.stan>
# Result of 19/10/2026 in benchmark-glm-gradient.rds (gain of the current
# model): 1.1x with collapsed rows, 2.5x with 6,000 occasion-level rows, 3.3x
# on the curiango data (8.5 -> 2.6 ms; 4.0x in an earlier run of the same
# day, the spread between runs being about 15%). Two variants not adopted: keeping
# the site effects on the false-positive rate as parameters fixed at zero
# (2.3 ms) and a grainsize of nsite / 32 in reduce_sum() (no difference).
# Dropping the scalar intercept of sites with a constant recorder, for the
# per-row vector form, made the gradient 1.6 times slower.
Sys.setenv(STAN_NUM_THREADS = 1)
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
stan_data <- function(s, occ, det, fp, random_recorder) {
  crywolf:::.build_stan_data(s$data, occ, det, fp, random_recorder,
                             fp_priors())$stan_data
}
cases <- list(
  "simulated 300 x 20, det ~ 1 (collapsed rows)" =
    stan_data(simulate_fp(nsite = 300, nocc = 20, conf_effect = 1.5,
                          conf_cor = 0.7, seed = 1), ~ x, ~ 1, ~ x, FALSE),
  "simulated 300 x 20, det ~ bio + aci, recorder effect" =
    stan_data(simulate_fp(nsite = 300, nocc = 20, a_bio = 0.5, a_aci = -0.3,
                          sd_recorder = 0.5, seed = 2),
              ~ x, ~ bio + aci, ~ 1, TRUE),
  "curiango, det ~ bio + aci, recorder effect" = curiango_bench_data()
)
models <- c(previous = commandArgs(trailingOnly = TRUE)[1], current = fp_stan_file())
mods <- purrr::map(models, function(path) {
  cmdstanr::cmdstan_model(path, cpp_options = list(stan_threads = TRUE),
                          quiet = TRUE, compile_model_methods = TRUE,
                          force_recompile = TRUE)
})

time_gradient <- function(m, sd) {
  set.seed(1)
  init <- list(beta = rnorm(sd$Kbeta, 0, 0.3), alpha = rnorm(sd$Kalpha, 0, 0.3),
               u_raw = rnorm(sd$nrec), sigma_u = 0.5, fp_gap = 3,
               gamma_fp = rnorm(sd$Kfp, 0, 0.3), delta0 = 0.3)
  f <- m$sample(data = sd, chains = 1, iter_sampling = 1, iter_warmup = 0,
                fixed_param = TRUE, init = list(init), seed = 1, refresh = 0,
                show_messages = FALSE, threads_per_chain = 1)
  f$init_model_methods()
  up <- as.numeric(f$unconstrain_draws(format = "draws_matrix")[1, ])
  for (i in 1:50) f$grad_log_prob(up)   # warm up
  secs <- system.time(for (i in 1:500) f$grad_log_prob(up))[["elapsed"]]
  c(ms_per_grad = 1000 * secs / 500, lp = f$log_prob(up))
}
res <- purrr::imap(cases, function(sd, label) {
  purrr::imap(mods, function(m, nm) {
    data.frame(data = label, model = nm, rows = sd$N2, threads = 1L,
               t(time_gradient(m, sd)))
  }) |> purrr::list_rbind()
}) |> purrr::list_rbind()
res$speedup <- stats::ave(res$ms_per_grad, res$data, FUN = function(v) v[1] / v)
print(res, digits = 3)
saveRDS(res, "dev/profile/benchmark-glm-gradient.rds")
