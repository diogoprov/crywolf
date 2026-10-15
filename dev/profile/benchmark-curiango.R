# Times the curiango fit (custom classifier, fp ~ 1, recorder effect) with the
# previous model and with the reduce_sum() model at several threads per chain.
# Short runs; the time per leapfrog step (one gradient evaluation) is
# comparable across runs even when the number of steps differs.
# Usage, from the package root:
#   Rscript dev/profile/benchmark-curiango.R <old/fp_occu.stan> [warmup sampling]
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
args <- commandArgs(trailingOnly = TRUE)
iw <- if (length(args) > 1) as.integer(args[2]) else 150
is <- if (length(args) > 2) as.integer(args[3]) else 100
sd <- curiango_bench_data()
old <- cmdstanr::cmdstan_model(args[1], quiet = TRUE)
new <- fp_compile()

run <- function(label, mod, threads = NULL) {
  t0 <- Sys.time()
  fit <- mod$sample(data = sd, chains = 4, parallel_chains = 4,
                    threads_per_chain = threads, iter_warmup = iw,
                    iter_sampling = is, adapt_delta = 0.95, seed = 1,
                    refresh = 0, show_messages = FALSE, save_warmup = TRUE)
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  leap <- sum(posterior::as_draws_df(fit$sampler_diagnostics(inc_warmup = TRUE))$n_leapfrog__)
  out <- data.frame(model = label, threads_per_chain = if (is.null(threads)) 1 else threads,
                    wall_s = wall, chain_s = mean(fit$time()$chains$total),
                    leapfrog = leap, ms_per_leapfrog = 1000 * 4 * mean(fit$time()$chains$total) / leap,
                    beta2 = mean(fit$draws("beta[2]")), alpha1 = mean(fit$draws("alpha[1]")))
  print(out, digits = 4)
  out
}
res <- list(
  run("previous", old),
  run("reduce_sum", new, 1),
  run("reduce_sum", new, 2),
  run("reduce_sum", new, 4)
) |> purrr::list_rbind()
res$speedup <- res$ms_per_leapfrog[1] / res$ms_per_leapfrog
print(res, digits = 4)
saveRDS(res, "dev/profile/benchmark-curiango.rds")
