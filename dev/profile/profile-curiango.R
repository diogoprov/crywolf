# Profiles the likelihood of the curiango fit with Stan `profile` blocks
# (dev/profile/fp_occu_profile.stan, a copy of inst/stan/fp_occu.stan).
# Usage, from the package root: Rscript dev/profile/profile-curiango.R
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
sd <- curiango_bench_data()
cat("N2 =", sd$N2, " nsite =", sd$nsite, " N1 =", sd$N1, "\n")
mod <- cmdstanr::cmdstan_model("dev/profile/fp_occu_profile.stan", quiet = TRUE)
t0 <- Sys.time()
fit <- mod$sample(data = sd, chains = 4, parallel_chains = 4,
                  iter_warmup = 150, iter_sampling = 100, adapt_delta = 0.95,
                  seed = 1, refresh = 0)
cat("wall time:", format(Sys.time() - t0), "\n")
prof <- fit$profiles()
tab <- purrr::imap(prof, function(p, chain) cbind(chain = chain, p)) |>
  purrr::list_rbind()
agg <- stats::aggregate(cbind(total_time, forward_time, reverse_time) ~ name,
                        data = tab, FUN = sum)
agg$share <- round(agg$total_time / sum(agg$total_time), 3)
print(agg[order(-agg$total_time), ])
print(prof[[1]][, c("name", "total_time", "autodiff_calls", "no_autodiff_calls")])
saveRDS(list(profiles = prof, summary = agg, time = fit$time()),
        "dev/profile/profile-curiango.rds")
