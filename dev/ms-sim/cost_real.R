# =============================================================================
# Cost of fp_occu_ms() on a dataset of real size: 20 species on the design of
# the curiango data (1753 site-nights x 56 one-minute recordings, about 1200
# verified recordings per species), simulated with simulate_fp_ms() (seed
# 7001, a_bio = 0.5, otherwise its defaults), fitted with centered = TRUE and
# centered_occ = TRUE.
#
# Three tasks, each saved in JOB_DIR as soon as it ends and skipped on a rerun:
#   threads   one chain of 20 + 20 iterations with THREADS threads per chain
#             (default 1, 2, 4, 8 and 16): seconds per gradient (leapfrog
#             step) in the sampling iterations, by threads.
#   full      4 chains of ITER + ITER iterations (default 1000),
#             adapt_delta = 0.95, WORKERS / 4 threads per chain: time,
#             leapfrog steps, divergences, treedepth, R-hat and ESS, and the
#             estimates of the parameters with a true value.
#   full-sdrec  the same fit on data simulated with recorder effects
#             (sd_recorder = 0.5, same seed). In "full" the true sd_recorder
#             is 0, at the boundary of the random effect the model
#             estimates; this task tells whether the divergences of "full"
#             come from that boundary.
#
# Run from the package root, with crywolf installed from this checkout:
#   WORKERS=32 Rscript dev/ms-sim/cost_real.R
# =============================================================================

library(crywolf)

workers <- as.integer(Sys.getenv("WORKERS", "32"))
iter <- as.integer(Sys.getenv("ITER", "1000"))
threads <- as.integer(strsplit(Sys.getenv("THREADS", "1,2,4,8,16"), ",")[[1]])
job_dir <- Sys.getenv("JOB_DIR", file.path("dev", "ms-sim", "cost-jobs"))
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)

sim <- function(sd_recorder) {
  simulate_fp_ms(nspecies = 20, nsite = 1753, nocc = 56, a_bio = 0.5,
                 frac_verified = 0.23, nver = 3, sd_recorder = sd_recorder,
                 seed = 7001)
}
s <- sim(0)
message("automated rows ", nrow(s$data$auto), ", verified rows ", nrow(s$data$verified))

fit <- function(chains, threads, warm, samp, s) {
  t0 <- Sys.time()
  m <- fp_occu_ms(s$data, occ = ~ x, det = ~ bio, fp = ~ 1, centered = TRUE,
                  centered_occ = TRUE, chains = chains, parallel_chains = chains,
                  threads_per_chain = threads, iter_warmup = warm,
                  iter_sampling = samp, adapt_delta = 0.95, seed = 1, refresh = 50)
  list(m = m, secs = as.numeric(difftime(Sys.time(), t0, units = "secs")))
}

f <- file.path(job_dir, "threads.rds")
if (!file.exists(f)) {
  res <- do.call(rbind, lapply(threads, function(th) {
    r <- fit(1, th, 20, 20, s)
    d <- r$m$fit$sampler_diagnostics(format = "df")
    data.frame(threads = th, seconds = r$secs,
               sampling_seconds = r$m$fit$time()$chains$sampling,
               leapfrog = sum(d$n_leapfrog__))
  }))
  res$sec_per_gradient <- res$sampling_seconds / res$leapfrog
  print(res)
  saveRDS(res, f)
}

full <- function(s) {
  r <- fit(4, max(1, workers %/% 4), iter, iter, s)
  d <- r$m$fit$sampler_diagnostics(format = "df")   # sampling iterations
  tab <- summary(r$m)$table
  out <- list(seconds = r$secs, threads_per_chain = max(1, workers %/% 4),
              iter = iter, chain_seconds = r$m$fit$time()$chains,
              sampler = as.data.frame(d), diagnostics = r$m$diagnostics,
              table = as.data.frame(tab), truth = s$truth$pars)
  print(out[c("seconds", "threads_per_chain", "chain_seconds", "diagnostics")])
  out
}

f <- file.path(job_dir, "full.rds")
if (!file.exists(f)) saveRDS(full(s), f)

f <- file.path(job_dir, "full-sdrec.rds")
if (!file.exists(f)) saveRDS(full(sim(0.5)), f)
