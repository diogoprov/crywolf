# Prototype of a site random effect on the true-detection probability
# (random_det, crywolf 1.8.0), fitted with dev/random-det/fp_occu_random_det.stan.
#
# Data: 300 sites, 20 occasions, 30% of sites verified (3 occasions each),
# occupancy slope -0.8, p11 = 0.6 and p10 = 0.05 at the reference, with
# among-site variation on the logit scale of sd sigma_det (detection) and
# sigma_fp (false positives), each 0 or 1: four scenarios.
# Models: base, random_det, random_fp, and both. The results in
# vignettes/site-effects-scenarios.rds have 100 replicates each, 1600 fits:
# 20 in a first run (about 2 h 50 min on 10 workers) and 80 more with
# REP_FROM=21 REPS=100 (about 8 h on 14 workers).
#
# Questions: does random_det restore the occupancy slope when detection
# varies among sites? Does it cost anything when it does not? Can random_det
# and random_fp be fitted together, or do the two site effects compete for
# the same detections (mixing, R-hat, coverage)?
#
# Run from the package root, with the package installed:
#   OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 \
#     Rscript dev/random-det/simulate-scenarios.R
# Environment: WORKERS (default 10), REPS (default 20), REP_FROM (default 1:
# run replicates REP_FROM to REPS, e.g. REP_FROM=21 REPS=100 to add 80 to a
# first run of 20), SMOKE=1 for a quick check of the code (one tiny fit;
# SMOKE_MODEL picks the model).
# Each fit is saved in dev/random-det/jobs/ as soon as it ends, and fits
# already there are skipped, so an interrupted run resumes where it stopped.
# Output: dev/random-det/scenarios.rds, all the fits in dev/random-det/jobs/
# plus, for replicates below REP_FROM, those of the previous scenarios.rds.
library(crywolf)
smoke <- identical(Sys.getenv("SMOKE"), "1")
workers <- as.integer(Sys.getenv("WORKERS", "10"))
reps <- if (smoke) 1L else as.integer(Sys.getenv("REPS", "20"))
rep_from <- if (smoke) 1L else as.integer(Sys.getenv("REP_FROM", "1"))
job_dir <- file.path("dev", "random-det", "jobs")

stan_file <- file.path("dev", "random-det", "fp_occu_random_det.stan")
mod <- cmdstanr::cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE),
                               quiet = TRUE)   # compile once, before the workers

# data with among-site variation in detection and in false positives
sim_scenario <- function(seed, sigma_det, sigma_fp, nsite = 300, nocc = 20,
                         frac_verified = 0.3, nver = 3) {
  crywolf:::.set_seed(seed)
  x <- stats::rnorm(nsite)
  z <- stats::rbinom(nsite, 1, stats::plogis(-0.8 * x))
  p11 <- stats::plogis(stats::qlogis(0.6) + stats::rnorm(nsite, 0, sigma_det))
  p10 <- stats::plogis(stats::qlogis(0.05) + stats::rnorm(nsite, 0, sigma_fp))
  p <- ifelse(z == 1, p11, p10)
  auto <- data.frame(site = rep(seq_len(nsite), each = nocc),
                     occasion = rep(seq_len(nocc), nsite))
  auto$y <- stats::rbinom(nrow(auto), 1, p[auto$site])
  vsites <- sort(sample.int(nsite, round(frac_verified * nsite)))
  verified <- data.frame(site = rep(vsites, each = nver))
  verified$y <- stats::rbinom(nrow(verified), 1, 0.6 * z[verified$site])
  withCallingHandlers(
    fp_data(auto, data.frame(site = seq_len(nsite), x = x), verified),
    crywolf_design_warning = function(w) invokeRestart("muffleWarning"))
}

fit_proto <- function(d, random_det, random_fp, seed) {
  sd <- crywolf:::.build_stan_data(d, ~ x, ~ 1, ~ 1, FALSE, fp_priors(),
                                   random_fp = random_fp)$stan_data
  sd$estimate_det_re <- as.integer(random_det)
  sd$det_site_sd <- 1
  t0 <- Sys.time()
  fit <- mod$sample(data = sd, seed = seed, chains = if (smoke) 1 else 4,
                    parallel_chains = 1, threads_per_chain = 1,
                    iter_warmup = if (smoke) 100 else 1000,
                    iter_sampling = if (smoke) 50 else 1000, adapt_delta = 0.95,
                    refresh = 0, show_messages = FALSE,
                    output_basename = crywolf:::.output_basename())
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  vars <- c("beta", "alpha", "fp_gap", "sd_det",
            if (random_fp) "sigma_fp")
  s <- posterior::summarise_draws(fit$draws(vars), "mean",
                                  ~ stats::quantile(.x, c(0.025, 0.975)), "rhat")
  b <- s[s$variable == "beta[2]", ]
  get <- function(v, col = "mean") if (v %in% s$variable) s[[col]][s$variable == v] else NA_real_
  data.frame(slope = b$mean, lower = b$`2.5%`, upper = b$`97.5%`,
             covered = b$`2.5%` <= -0.8 & -0.8 <= b$`97.5%`,
             sd_det_est = get("sd_det"), sd_fp_est = get("sigma_fp[1]"),
             max_rhat = max(s$rhat, na.rm = TRUE),
             divergences = sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent),
             seconds = secs)
}

scenarios <- expand.grid(sigma_det = c(0, 1), sigma_fp = c(0, 1))
models <- data.frame(model = c("base", "random_det", "random_fp", "both"),
                     random_det = c(FALSE, TRUE, FALSE, TRUE),
                     random_fp = c(FALSE, FALSE, TRUE, TRUE))
jobs <- merge(merge(scenarios, models), data.frame(rep = seq_len(reps)))
jobs$seed <- seq_len(nrow(jobs))
jobs <- jobs[jobs$rep >= rep_from, ]
if (smoke) {
  smoke_model <- Sys.getenv("SMOKE_MODEL", "random_det")
  jobs <- jobs[jobs$sigma_det == 1 & jobs$sigma_fp == 1 & jobs$model == smoke_model, ][1, ]
}

run_one <- function(sigma_det, sigma_fp, model, random_det, random_fp, rep, seed) {
  f <- file.path(job_dir, sprintf("det%g-fp%g-%s-rep%03d.rds", sigma_det, sigma_fp,
                                  model, rep))
  if (!smoke && file.exists(f)) return(readRDS(f))
  # every model of a replicate sees the same data: the data seed is the
  # replicate's, the fit seed the job's
  d <- sim_scenario(1000 * rep + 10 * sigma_det + sigma_fp, sigma_det, sigma_fp,
                    nsite = if (smoke) 40 else 300, nocc = if (smoke) 4 else 20)
  out <- tryCatch(fit_proto(d, random_det, random_fp, seed),
                  error = function(e) data.frame(slope = NA, error = conditionMessage(e)))
  res <- data.frame(sigma_det, sigma_fp, model, rep, out)
  if (!smoke) saveRDS(res, f)
  res
}

if (!smoke) dir.create(job_dir, showWarnings = FALSE)
if (smoke || workers <= 1) {
  res <- purrr::pmap(jobs, run_one) |> purrr::list_rbind()
} else {
  future::plan(future::multisession, workers = workers)
  # one job at a time per worker, so that slow fits do not hold back a chunk
  res <- furrr::future_pmap(jobs, run_one,
                            .options = furrr::furrr_options(seed = TRUE, chunk_size = 1)) |>
    purrr::list_rbind()
}
if (smoke) {
  print(res)
} else {
  out_file <- file.path("dev", "random-det", "scenarios.rds")
  if (rep_from > 1) {
    old <- readRDS(out_file)
    # list_rbind() fills a column that one side lacks (error) with NA
    res <- purrr::list_rbind(list(old[old$rep < rep_from, ], res))
  }
  saveRDS(res, out_file)
  # summary per scenario and model; failed fits (NA) are left out
  agg <- purrr::map(split(res, res[c("sigma_det", "sigma_fp", "model")], drop = TRUE),
                    function(g) {
    m <- function(v) mean(g[[v]], na.rm = TRUE)
    data.frame(sigma_det = g$sigma_det[1], sigma_fp = g$sigma_fp[1], model = g$model[1],
               fits = sum(!is.na(g$slope)), coverage = m("covered"),
               bias = m("slope") + 0.8, sd_det_est = m("sd_det_est"),
               sd_fp_est = m("sd_fp_est"),
               rhat_over_1.01 = mean(g$max_rhat > 1.01, na.rm = TRUE),
               divergent = mean(g$divergences > 0, na.rm = TRUE),
               seconds = stats::median(g$seconds, na.rm = TRUE))
  }) |> purrr::list_rbind()
  print(agg, digits = 3)
}
