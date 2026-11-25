# Pilot of the JSDM prototype (README.md). Scenarios, each over REPS datasets
# of 100 sites x 10 species x 10 occasions with verified data (3 occasions)
# at 30% of the sites:
#   A1  two true factors, fit K = 2 with the diagonal free (sign fixed after
#       sampling, as jsdmstan)
#   A2  the same data, fit K = 2 with a positive diagonal
#   B1  no factors but detection correlated among species (shared site effect
#       on detection, sd 1), fit K = 1 without a detection effect
#   B2  the same data, fit K = 1 with the shared detection effect
#   C   no factors and no correlated detection, fit K = 1 (null)
# Recovery of Lambda Lambda' (off-diagonal): correlation and RMSE between the
# posterior mean and the truth; for the null scenarios, the size of the
# estimated off-diagonal entries. Convergence: largest R-hat, smallest bulk
# ESS, divergences.
#   REPS=3 WORKERS=10 Rscript dev/jsdm/pilot.R
library(cmdstanr)
source("dev/jsdm/simulate_jsdm_fp.R")
reps <- as.integer(Sys.getenv("REPS", "3"))
workers <- as.integer(Sys.getenv("WORKERS", "10"))
out_file <- Sys.getenv("OUT", "dev/jsdm/pilot.rds")
mod <- cmdstan_model("dev/jsdm/jsdm_fp.stan")

scen <- data.frame(
  id = c("A1", "A2", "B1", "B2", "C"),
  K_true = c(2, 2, 0, 0, 0), sigma_d = c(0, 0, 1, 1, 0),
  K_fit = c(2, 2, 1, 1, 1), diag_pos = c(0, 1, 0, 0, 0),
  det_re = c(0, 0, 0, 1, 0))
tasks <- expand.grid(rep = seq_len(reps), s = seq_len(nrow(scen)))

run <- function(rep, s) {
  sc <- scen[s, ]
  # scenarios sharing data (A1/A2, B1/B2) use the same seed
  seed <- 1000 * match(sc$K_true * 10 + sc$sigma_d, c(20, 1, 0)) + rep
  sim <- simulate_jsdm_fp(K = sc$K_true, sigma_d = sc$sigma_d, seed = seed)
  dat <- c(sim$data, list(K = sc$K_fit, diag_pos = sc$diag_pos, det_re = sc$det_re))
  t0 <- Sys.time()
  fit <- mod$sample(data = dat, chains = 4, parallel_chains = 1,
                    iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.95,
                    seed = seed, refresh = 0, show_messages = FALSE,
                    output_dir = tempdir())
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  J <- sim$data$J
  ll <- matrix(fit$summary("LL", "mean")$mean, J, J)
  off <- lower.tri(ll)
  key <- c("sigma_L", "mu_b0", "mu_a", "mu_lg", if (sc$det_re) "sigma_d[1]")
  sm <- fit$summary()
  dg <- fit$diagnostic_summary(quiet = TRUE)
  lam <- fit$summary("Lambda")
  data.frame(
    scenario = sc$id, rep = rep, seed = seed, minutes = secs / 60,
    divergent = sum(dg$num_divergent), treedepth = sum(dg$num_max_treedepth),
    rhat_max = max(sm$rhat, na.rm = TRUE),
    rhat_LL = max(fit$summary("LL")$rhat, na.rm = TRUE),
    rhat_Lambda = max(lam$rhat, na.rm = TRUE),
    ess_LL = min(fit$summary("LL")$ess_bulk, na.rm = TRUE),
    cor_LL = if (sc$K_true > 0) stats::cor(ll[off], sim$truth$LL[off]) else NA,
    rmse_LL = sqrt(mean((ll[off] - sim$truth$LL[off])^2)),
    mean_abs_LL = mean(abs(ll[off])),
    sigma_L = sm$mean[sm$variable == "sigma_L"],
    sigma_d = if (sc$det_re) sm$mean[sm$variable == "sigma_d[1]"] else NA,
    occ_mean = mean(sim$truth$occupancy))
}

future::plan(future::multisession, workers = workers)
res <- do.call(rbind, furrr::future_pmap(tasks, run,
  .options = furrr::furrr_options(seed = NULL, scheduling = Inf)))
saveRDS(res, out_file)
print(res, digits = 3, row.names = FALSE)
