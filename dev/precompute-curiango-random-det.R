# The curiango data of vignettes/markov.Rmd (custom classifier, threshold
# 0.40) fitted with a site random effect on true detection (random_det), with
# and without the Markov term, next to the two fits of
# dev/precompute-markov.R refitted here so that one loo_compare() holds all
# four. Saves posterior summaries, diagnostics, posterior predictive checks,
# LOO (integrated over the site effects for random_det) and runs of
# detections in data simulated from each fit (runs() is in
# dev/curiango-common.R). Four fits of 4 chains x 4 threads, run one after
# another; the two without random_det took about 6 minutes each.
# Run from the package root after installing the package:
#   Rscript dev/precompute-curiango-random-det.R
library(crywolf)
source(file.path("dev", "curiango-common.R"))
auto <- cur_occasions
auto$y <- as.integer(cur_detections$confidence >= 0.40)
d <- suppressWarnings(fp_data(auto, cur_sites, cur_verified, site = "site_night",
                              recorder = "recorder"))
fit <- function(markov, random_det) {
  t0 <- Sys.time()
  f <- fp_occu(d, occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1, markov = markov,
               random_recorder = TRUE, random_det = random_det, seed = 1, chains = 4,
               parallel_chains = 4, threads_per_chain = 4, refresh = 0,
               show_messages = FALSE)
  list(fit = f, minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
r <- list(independent = fit(NULL, FALSE), markov = fit("time", FALSE),
          random_det = fit(NULL, TRUE), random_det_markov = fit("time", TRUE))

sim_runs <- purrr::map(r, function(x) {
  sims <- simulate(x$fit, nsim = 50, seed = 1, truth = "posterior")
  purrr::map(sims, function(s) as.numeric(runs(s$data$auto))) |>
    do.call(what = rbind)
})
loos <- purrr::map(r, function(x) loo(x$fit))
out <- list(
  summary = purrr::map(r, function(x) summary(x$fit)$table),
  diagnostics = purrr::map(r, function(x) x$fit$diagnostics),
  minutes = purrr::map_dbl(r, "minutes"),
  ppc = purrr::map(r, function(x) {
    p <- fp_ppc(x$fit, ndraws = 300, seed = 1)
    list(p_value = p$p_value, fit_obs = p$fit_obs, fit_rep = p$fit_rep)
  }),
  pareto_k = purrr::map(loos, function(l) l$diagnostics$pareto_k),
  loo = as.data.frame(loo::loo_compare(loos)),
  runs_observed = runs(d$auto),
  runs_simulated = sim_runs
)
saveRDS(crywolf:::.stamp(out), file.path("vignettes", "curiango-random-det-results.rds"),
        compress = "xz")
