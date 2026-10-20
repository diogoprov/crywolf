# Times one gradient evaluation of the log density on the curiango data, at a
# fixed point, for the previous model and the current one (200 calls through
# the model methods of 'cmdstanr', one thread). Deterministic, unlike sampler
# timings, which vary with the load of the machine. Usage, from the package
# root:
#   Rscript dev/profile/benchmark-glm-gradient.R <previous fp_occu.stan>
# Result of 19/10/2026 in benchmark-glm-gradient.rds: 8.8 ms -> 2.2 ms per
# gradient (4.0x). Two variants tried the same day and not adopted: keeping
# the site effects on the false-positive rate as parameters fixed at zero
# (2.3 ms) and a grainsize of nsite / 32 in reduce_sum() (no difference).
# Dropping the scalar intercept of sites with a constant recorder, for the
# per-row vector form, made the gradient 1.6 times slower.
Sys.setenv(STAN_NUM_THREADS = 1)
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
sd <- curiango_bench_data()
models <- c(previous = commandArgs(trailingOnly = TRUE)[1], current = fp_stan_file())

res <- purrr::imap(models, function(path, label) {
  m <- cmdstanr::cmdstan_model(path, cpp_options = list(stan_threads = TRUE),
                               quiet = TRUE, compile_model_methods = TRUE,
                               force_recompile = TRUE)
  set.seed(1)
  init <- list(beta = rnorm(sd$Kbeta, 0, 0.3), alpha = c(-2.5, 0.5, -1),
               u_raw = rnorm(sd$nrec), sigma_u = 0.7, fp_gap = 3,
               gamma_fp = numeric(0), delta0 = 0.3)
  f <- m$sample(data = sd, chains = 1, iter_sampling = 1, iter_warmup = 0,
                fixed_param = TRUE, init = list(init), seed = 1, refresh = 0,
                show_messages = FALSE, threads_per_chain = 1)
  f$init_model_methods()
  up <- as.numeric(f$unconstrain_draws(format = "draws_matrix")[1, ])
  g <- f$grad_log_prob(up)
  for (i in 1:20) f$grad_log_prob(up)   # warm up
  secs <- system.time(for (i in 1:200) f$grad_log_prob(up))[["elapsed"]]
  data.frame(model = label, threads = 1L, ms_per_grad = 1000 * secs / 200,
             lp = f$log_prob(up), gmax = max(abs(g)))
}) |> purrr::list_rbind()
res$speedup_vs_previous <- res$ms_per_grad[1] / res$ms_per_grad
print(res, digits = 4)
saveRDS(res, "dev/profile/benchmark-glm-gradient.rds")
