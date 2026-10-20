# Times short sampler runs on the curiango data (custom classifier, fp ~ 1,
# recorder effect) with the previous model and the current one, at one and
# four threads per chain: four chains, 150 warm-up and 100 sampling
# iterations, seed 1. The time per leapfrog step (one gradient evaluation)
# is comparable across runs even when the number of steps differs, but it
# varies with the load of the machine; benchmark-glm-gradient.R is the
# deterministic measure. Usage, from the package root:
#   Rscript dev/profile/benchmark-glm-sampler.R <previous fp_occu.stan>
# Result of 19/10/2026 in benchmark-glm-sampler.rds: 2.3x less time per
# leapfrog step with one thread, 2.6x with four.
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
sd <- curiango_bench_data()
previous <- cmdstanr::cmdstan_model(commandArgs(trailingOnly = TRUE)[1],
                                    cpp_options = list(stan_threads = TRUE),
                                    quiet = TRUE)
current <- fp_compile()

run <- function(label, mod, threads) {
  t0 <- Sys.time()
  fit <- mod$sample(data = sd, chains = 4, parallel_chains = 4,
                    threads_per_chain = threads, iter_warmup = 150,
                    iter_sampling = 100, adapt_delta = 0.95, seed = 1,
                    refresh = 0, show_messages = FALSE, save_warmup = TRUE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  leap <- sum(posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))$n_leapfrog__)
  chain_s <- mean(fit$time()$chains$total)
  out <- data.frame(model = label, threads = threads, wall_s = wall,
                    chain_s = chain_s, leapfrog = leap,
                    ms_per_leapfrog = 1000 * 4 * chain_s / leap,
                    beta2 = mean(fit$draws("beta[2]")),
                    alpha1 = mean(fit$draws("alpha[1]")))
  print(out, digits = 4)
  out
}
res <- purrr::list_rbind(list(
  run("previous", previous, 1), run("current", current, 1),
  run("previous", previous, 4), run("current", current, 4)
))
res$speedup_vs_previous_1thread <- res$ms_per_leapfrog[1] / res$ms_per_leapfrog
res$chains <- 4L
res$iter_warmup <- 150L
res$iter_sampling <- 100L
print(res, digits = 4)
saveRDS(res, "dev/profile/benchmark-glm-sampler.rds")
