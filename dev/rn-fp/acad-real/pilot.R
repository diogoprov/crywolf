# Pilot for the FrogNet request: realistic ACAD design (21 sites, 11-19
# occasions, every occasion confirmed, no verified data) fitted with
# fp_occu_rn(), under four priors on sigma_r: half-Normal(0, 0.5 / 1 / 2) and
# half-Student-t(3, 0, 1) (dev copy fp_occu_rn_t3.stan). Times one fit per
# prior on the slowest cell (psi 0.97, p10 0.15, sigma_r 1).
devtools::load_all(quiet = TRUE)
set.seed(500000)
nocc <- sample(11:19, 21, replace = TRUE)
t3 <- cmdstanr::cmdstan_model("dev/rn-fp/acad-real/fp_occu_rn_t3.stan",
                              cpp_options = list(stan_threads = TRUE))
fit_one <- function(s, prior) {
  pr <- fp_priors(det_site_sd = if (prior == "t3") 1 else as.numeric(prior))
  if (prior == "t3") {
    orig <- crywolf:::.compile_model
    assignInNamespace(".compile_model", function(...) t3, "crywolf")
    on.exit(assignInNamespace(".compile_model", orig, "crywolf"))
  }
  t0 <- Sys.time()
  m <- fp_occu_rn(s$data, occ = ~ x, random_det = TRUE, random_abund = TRUE,
                  n_max = 30, priors = pr, chains = 4, parallel_chains = 4,
                  threads_per_chain = 1, seed = 1, refresh = 0,
                  show_messages = FALSE)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  d <- m$fit$summary(c("sigma_det", "gamma0", "alpha", "sigma_lam"))
  list(prior = prior, minutes = el, summary = d,
       diag = m$fit$diagnostic_summary(quiet = TRUE))
}
s <- simulate_fp_rn(nsite = 21, nocc = nocc, psi = 0.97, r = 0.27, p10 = 0.15,
                    sigma_r = 1, sigma_lam = 0.7, frac_verified = 0,
                    frac_confirmed = 1, seed = 500001)
res <- lapply(c("0.5", "1", "2", "t3"), function(p) fit_one(s, p))
for (r in res) {
  cat("\nprior", r$prior, "minutes", round(r$minutes, 2),
      "divergences", sum(r$diag$num_divergent), "\n")
  print(r$summary[, c("variable", "median", "q5", "q95", "rhat", "ess_bulk")])
}
saveRDS(res, "dev/rn-fp/acad-real/pilot.rds")
