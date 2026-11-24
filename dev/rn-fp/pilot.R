# Step 1 pilot: three data sets (site variation in r only, in lambda only,
# and in both), each fit with site effects on both; reports the posteriors of
# sigma_r and sigma_lam, diagnostics and run time. Run from the package root: Rscript dev/rn-fp/pilot.R
source("dev/rn-fp/simulate_rn_fp.R")
mod <- cmdstanr::cmdstan_model("dev/rn-fp/rn_fp.stan")
# Gauss-Hermite nodes and weights for a standard normal (Golub-Welsch)
gh <- function(nq) {
  J <- matrix(0, nq, nq)
  b <- sqrt(seq_len(nq - 1) / 2)
  J[cbind(1:(nq - 1), 2:nq)] <- b
  J[cbind(2:nq, 1:(nq - 1))] <- b
  ev <- eigen(J, symmetric = TRUE)
  list(nq = nq, gh_x = sqrt(2) * ev$values,
       gh_logw = log(ev$vectors[1, ]^2))
}
GH <- gh(9)
pars <- c("b0", "b1", "a0", "sigma_r_out", "sigma_lam_out", "p10", "rv")
fit_one <- function(sigma_r, sigma_lam, seed) {
  s <- simulate_rn_fp(nsite = 200, nocc = 10, sigma_r = sigma_r,
                      sigma_lam = sigma_lam, seed = seed)
  t0 <- Sys.time()
  f <- mod$sample(data = c(s$data, list(Nmax = 20L, site_effect = 1L,
                                        lam_effect = 1L), GH),
                  chains = 4, parallel_chains = 4, iter_warmup = 500,
                  iter_sampling = 500, seed = seed, refresh = 0,
                  show_messages = FALSE)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d <- f$diagnostic_summary(quiet = TRUE)
  sm <- f$summary(pars, "mean", ~quantile(.x, c(0.05, 0.95)), "rhat", "ess_bulk")
  sm$truth <- with(s$truth, c(b0, b1, a0, sigma_r, sigma_lam, p10, rv))
  list(summary = sm, secs = secs, divergent = sum(d$num_divergent),
       treedepth = sum(d$num_max_treedepth),
       sigma_r_draws = as.vector(f$draws("sigma_r_out")),
       sigma_lam_draws = as.vector(f$draws("sigma_lam_out")))
}
# A: only r varies; B: only lambda varies; C: both
res <- list(A = fit_one(1, 0, 21), B = fit_one(0, 0.7, 22), C = fit_one(1, 0.7, 23))
for (nm in names(res)) {
  cat("\n==", nm, "| tempo", round(res[[nm]]$secs), "s | divergências",
      res[[nm]]$divergent, "| treedepth", res[[nm]]$treedepth, "\n")
  print(as.data.frame(res[[nm]]$summary), digits = 3)
}
saveRDS(res, "dev/rn-fp/pilot.rds")
