# Pilot of classifier confusion (README.md). No latent factors and no shared
# detection effect in the truth; species 1, 3 and 5 receive automated
# detections from the presence of species 2, 4 and 6, with per-occasion
# probability `conf`. Verified occasions have no confusion. All fits use K = 1
# with the shared detection effect (the configuration the package would use):
#   D1  conf 0.15, verified data at 30% of the sites
#   D2  conf 0.30, verified data at 30% of the sites
#   D3  conf 0.30, verified data at 10% of the sites
# Spurious association: posterior mean of LL for the confused pairs (1-2, 3-4,
# 5-6) against the other pairs. Occupancy bias: posterior mean of the species
# intercept minus the truth, for confused and other species.
#   REPS=3 WORKERS=9 Rscript dev/jsdm/pilot_confusion.R
library(cmdstanr)
source("dev/jsdm/simulate_jsdm_fp.R")
reps <- as.integer(Sys.getenv("REPS", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "9"))
out_file <- Sys.getenv("OUT", "dev/jsdm/pilot_confusion.rds")
mod <- cmdstan_model("dev/jsdm/jsdm_fp.stan")

scen <- data.frame(id = c("D1", "D2", "D3"), conf = c(0.15, 0.3, 0.3),
                   frac_verified = c(0.3, 0.3, 0.1))
tasks <- expand.grid(rep = seq_len(reps), s = seq_len(nrow(scen)))
J <- 10
pairs <- cbind(c(1, 3, 5), c(2, 4, 6))

run <- function(rep, s) {
  sc <- scen[s, ]
  seed <- 4000 + 100 * s + rep
  cm <- matrix(0, J, J)
  cm[pairs] <- sc$conf
  sim <- simulate_jsdm_fp(J = J, K = 0, conf = cm,
                          frac_verified = sc$frac_verified, seed = seed)
  dat <- c(sim$data, list(K = 1, diag_pos = 0, det_re = 1))
  t0 <- Sys.time()
  fit <- mod$sample(data = dat, chains = 4, parallel_chains = 1,
                    iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
                    seed = seed, refresh = 0, show_messages = FALSE,
                    output_dir = tempdir())
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ll <- matrix(fit$summary("LL", "mean")$mean, J, J)
  conf_pair <- matrix(FALSE, J, J)
  conf_pair[pairs] <- TRUE
  conf_pair <- conf_pair | t(conf_pair)
  off <- lower.tri(ll)
  dr <- fit$draws(c("mu_b0", "sd_b0", "z_b0"), format = "draws_matrix")
  zb <- unclass(dr[, paste0("z_b0[", seq_len(J), "]")])
  b0 <- colMeans(as.numeric(dr[, "mu_b0"]) + as.numeric(dr[, "sd_b0"]) * zb)
  bias <- b0 - sim$truth$b0
  sm <- fit$summary()
  dg <- fit$diagnostic_summary(quiet = TRUE)
  data.frame(
    scenario = sc$id, rep = rep, seed = seed, minutes = secs / 60,
    divergent = sum(dg$num_divergent),
    rhat_LL = max(fit$summary("LL")$rhat, na.rm = TRUE),
    ess_LL = min(fit$summary("LL")$ess_bulk, na.rm = TRUE),
    LL_conf = mean(ll[off & conf_pair]),
    LL_other = mean(ll[off & !conf_pair]),
    absLL_other = mean(abs(ll[off & !conf_pair])),
    sigma_L = sm$mean[sm$variable == "sigma_L"],
    sigma_d = sm$mean[sm$variable == "sigma_d[1]"],
    b0_bias_conf = mean(bias[pairs[, 1]]),
    b0_bias_other = mean(bias[-pairs[, 1]]))
}

future::plan(future::multisession, workers = workers)
res <- do.call(rbind, furrr::future_pmap(tasks, run,
  .options = furrr::furrr_options(seed = NULL, scheduling = Inf)))
saveRDS(res, out_file)
print(res, digits = 3, row.names = FALSE)
