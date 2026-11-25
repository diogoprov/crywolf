# =============================================================================
# dev/rn-fp/sim/run_sbc.R
# Simulation-based calibration (Talts et al. 2018) of fp_occu_rn(): parameters
# drawn from the priors the fit uses, data from simulate_fp_rn(), and the rank
# of each true value among thinned posterior draws.
#
#   priors fp_priors(intercept_sd = 0.75, slope_sd = 0.3, gap_mean = 2,
#          gap_sd = 0.5, det_site_sd = 1, abund_site_sd = 0.3), narrower than
#          the package defaults so that abundance rarely reaches n_max = 50:
#          in a prior check of 4000 datasets, 0.1% had a site with
#          P(N > 50) > 0.001 (p_above_nmax and n_max_true are kept to check
#          this). The fit uses the same priors.
#   draws  beta[1] ~ N(0, 0.75) (log lambda at x = 0), beta[2] ~ N(0, 0.3),
#          alpha[1] ~ N(0, 0.75) (logit r), gamma0 ~ N(-2, sqrt(0.75^2 +
#          0.5^2)) (logit p10), delta0 ~ N(0, 0.75) (logit rv),
#          sigma_det ~ HN(0, 1), sigma_lam ~ HN(0, 0.3). A draw with psi
#          numerically 1 (beta[1] > 3.6) is drawn again.
#   designs A: 50 sites x 10 occasions, every occasion confirmed, no verified
#             data (the ACAD scheme; delta0 is not in the data and left out).
#           B: 50 sites x 10 occasions, 20% of occasions confirmed, verified
#             data (3 occasions) at 30% of sites.
#   fits   occ = ~ x, random_det, random_abund, n_max 50, 4 chains of
#          1000 + 1000, ADAPT_DELTA (default 0.99; the pilot
#          at 0.95 had divergences in most designs). Ranks over 100 draws thinned evenly
#          from the 4000 (0..100), with ESS kept to check thinning.
#
# REPS replicates per design (default 500); seed 700000 + 10000 x design + rep.
# Each replicate is saved in JOBS (default dev/rn-fp/sim/jobs-sbc) and skipped
# on a rerun. From the package root, with crywolf installed (LOAD_ALL = 1 for
# the checkout, local tests only):
#   WORKERS=64 Rscript dev/rn-fp/sim/run_sbc.R
# =============================================================================
if (Sys.getenv("LOAD_ALL") == "1") devtools::load_all(quiet = TRUE) else library(crywolf)
reps <- as.integer(Sys.getenv("REPS", "500"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.99"))
job_dir <- Sys.getenv("JOBS", "dev/rn-fp/sim/jobs-sbc")
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)

pr <- fp_priors(intercept_sd = 0.75, slope_sd = 0.3, gap_mean = 2, gap_sd = 0.5,
                det_site_sd = 1, abund_site_sd = 0.3)
designs <- data.frame(design = c("A", "B"), frac_confirmed = c(1, 0.2),
                      frac_verified = c(0, 0.3))
pars <- c("beta[1]", "beta[2]", "alpha[1]", "gamma0", "delta0",
          "sigma_det[1]", "sigma_lam[1]")

tasks <- expand.grid(rep = seq_len(reps), design = seq_len(nrow(designs)))
tasks$seed <- 700000 + 10000 * tasks$design + tasks$rep
tasks$file <- file.path(job_dir, sprintf("sbc%s-rep%04d.rds",
                                         designs$design[tasks$design], tasks$rep))
tasks <- tasks[!file.exists(tasks$file), ]
if (nzchar(Sys.getenv("MAX_TASKS"))) tasks <- head(tasks, as.integer(Sys.getenv("MAX_TASKS")))
message(nrow(tasks), " tasks to run on ", workers, " workers")
invisible(crywolf:::.compile_model("fp_occu_rn"))

draw_truth <- function() {
  repeat {
    th <- c(`beta[1]` = stats::rnorm(1, 0, 0.75), `beta[2]` = stats::rnorm(1, 0, 0.3),
            `alpha[1]` = stats::rnorm(1, 0, 0.75),
            gamma0 = stats::rnorm(1, -2, sqrt(0.75^2 + 0.5^2)),
            delta0 = stats::rnorm(1, 0, 0.75),
            `sigma_det[1]` = abs(stats::rnorm(1, 0, 1)),
            `sigma_lam[1]` = abs(stats::rnorm(1, 0, 0.3)))
    psi <- -expm1(-exp(th[["beta[1]"]]))
    if (psi < 1) return(th)
  }
}

run_task <- function(rep, design, seed, file) {
  dg <- designs[design, ]
  set.seed(seed)
  th <- draw_truth()
  s <- simulate_fp_rn(nsite = 50, nocc = 10, psi = -expm1(-exp(th[["beta[1]"]])),
                      b1 = th[["beta[2]"]], r = stats::plogis(th[["alpha[1]"]]),
                      p10 = stats::plogis(th[["gamma0"]]),
                      sigma_r = th[["sigma_det[1]"]], sigma_lam = th[["sigma_lam[1]"]],
                      rv = stats::plogis(th[["delta0"]]),
                      frac_verified = dg$frac_verified, nver = 3,
                      frac_confirmed = dg$frac_confirmed, seed = seed + 1)
  out_dir <- tempfile("sbc"); dir.create(out_dir)
  t0 <- Sys.time()
  m <- suppressWarnings(fp_occu_rn(
    s$data, occ = ~ x, random_det = TRUE, random_abund = TRUE, n_max = 50,
    priors = pr, chains = 4, parallel_chains = 1, threads_per_chain = 1,
    seed = seed, refresh = 0, show_messages = FALSE, output_dir = out_dir,
    adapt_delta = adapt_delta))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  use <- if (dg$frac_verified > 0) pars else setdiff(pars, "delta0")
  dr <- posterior::as_draws_matrix(m$fit$draws(use))
  thin <- dr[round(seq(40, nrow(dr), length.out = 100)), , drop = FALSE]
  ranks <- vapply(use, function(p) sum(thin[, p] < th[[p]]), numeric(1))
  sm <- as.data.frame(m$fit$summary(use, "median", ~quantile(.x, c(0.05, 0.95)),
                                    "rhat", "ess_bulk", "ess_tail"))
  d <- m$fit$diagnostic_summary(quiet = TRUE)
  res <- list(design = dg$design, adapt_delta = adapt_delta, rep = rep, seed = seed, truth = th, ranks = ranks,
              summary = sm, secs = secs,
              rhat_max = max(m$fit$summary()$rhat, na.rm = TRUE),
              divergent = sum(d$num_divergent), treedepth = sum(d$num_max_treedepth),
              p_above_nmax = stats::median(m$fit$draws("p_above_nmax")),
              n_max_true = max(s$truth$N), n_occupied = sum(s$truth$z))
  unlink(out_dir, recursive = TRUE)
  tmp <- paste0(file, ".tmp"); saveRDS(res, tmp); file.rename(tmp, file)
  invisible(NULL)
}

if (nrow(tasks) > 0) {
  if (workers > 1) future::plan(future::multisession, workers = workers) else future::plan(future::sequential)
  furrr::future_pwalk(tasks, run_task,
                      .options = furrr::furrr_options(seed = NULL, scheduling = Inf))
}
message("done: ", length(list.files(job_dir, "\\.rds$")), " job files")
