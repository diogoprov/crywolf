# Pre-computes the results of vignettes/markov.Rmd: the curiango data (custom
# classifier, threshold 0.40) fitted with and without the Markov term, their
# posterior predictive checks, LOO, and runs of detections in data simulated
# from each fit. About 30 minutes on 16 cores.
# Run from the package root after installing the package:
#   Rscript dev/precompute-markov.R
library(crywolf)
source(file.path("dev", "curiango-common.R"))
auto <- cur_occasions
auto$y <- as.integer(cur_detections$confidence >= 0.40)
d <- suppressWarnings(fp_data(auto, cur_sites, cur_verified, site = "site_night",
                              recorder = "recorder"))
fit <- function(markov) {
  t0 <- Sys.time()
  f <- fp_occu(d, occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1, markov = markov,
               random_recorder = TRUE, seed = 1, chains = 4, parallel_chains = 4,
               threads_per_chain = 4, refresh = 0, show_messages = FALSE)
  list(fit = f, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
r <- list(independent = fit(NULL), markov = fit("time"))

# runs of consecutive detections within a night (15-minute steps)
runs <- function(a) {
  a <- a[order(a$site_night, a$time), ]
  same <- c(FALSE, a$site_night[-1] == a$site_night[-nrow(a)])
  step <- c(NA, diff(as.numeric(a$time)))
  adjacent <- same & step == 900
  start <- a$y == 1 & !(adjacent & c(0, a$y[-nrow(a)]) == 1)
  id <- cumsum(start)
  len <- tabulate(id[a$y == 1])
  table(factor(pmin(len, 5), levels = 1:5, labels = c(1:4, "5+")))
}
obs_runs <- runs(d$auto)
sim_runs <- purrr::map(r, function(x) {
  sims <- simulate(x$fit, nsim = 50, seed = 1, truth = "posterior")
  purrr::map(sims, function(s) as.numeric(runs(s$data$auto))) |>
    do.call(what = rbind)
})
out <- list(
  summary = purrr::map(r, function(x) summary(x$fit)$table),
  diagnostics = purrr::map(r, function(x) x$fit$diagnostics),
  minutes = purrr::map_dbl(r, "minutes"),
  ppc = purrr::map(r, function(x) {
    p <- fp_ppc(x$fit, ndraws = 300, seed = 1)
    list(p_value = p$p_value, fit_obs = p$fit_obs, fit_rep = p$fit_rep)
  }),
  loo = as.data.frame(loo::loo_compare(purrr::map(r, function(x) loo(x$fit)))),
  runs_observed = obs_runs,
  runs_simulated = sim_runs,
  n_occasions = nrow(d$auto),
  n_gaps = sum(!crywolf:::.build_stan_data(d, ~ 1, ~ 1, ~ 1, FALSE, fp_priors(),
                                          markov = "time")$occasions$has_prev) -
    nrow(d$sites)
)
saveRDS(out, file.path("vignettes", "markov-results.rds"), compress = "xz")
