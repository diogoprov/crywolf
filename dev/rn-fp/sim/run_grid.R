# =============================================================================
# dev/rn-fp/sim/run_grid.R
# General simulation grid for fp_occu_rn(): bias and coverage of the
# parameters, and how well the site effect on detection (sigma_r) is told
# apart from the site effect on abundance (sigma_lam), across designs.
#
#   cells  nsite {25, 50, 100} x nocc {10, 20} x annotation {verified,
#          partial, full} x sigma_r {0, 1} x sigma_lam {0, 0.7}: 72 cells,
#          REPS replicates each (default 100). Seed 800000 + 1000 x cell + rep.
#          verified: no confirmed occasions, 3 verified occasions at 30% of
#                    sites (rv = 0.5).
#          partial:  20% of automated occasions confirmed, plus verified as
#                    above.
#          full:     every occasion confirmed, no verified data (the ACAD
#                    scheme).
#   fixed  psi 0.6, b1 0.5 on a standardized covariate x, r 0.3, p10 0.05.
#   fits   occ = ~ x, random_det, random_abund (both always estimated), the
#          package default priors, n_max 25, 4 chains of 1000 + 1000,
#          ADAPT_DELTA (default 0.99; the pilot
#          at 0.95 had divergences in most designs). Kept: summaries of the main quantities,
#          diagnostics, and 1000 thinned draws of sigma_r and sigma_lam.
#
# Each replicate is saved in JOBS (default dev/rn-fp/sim/jobs-grid) and
# skipped on a rerun. From the package root, with crywolf installed
# (LOAD_ALL = 1 for the checkout, local tests only):
#   WORKERS=64 Rscript dev/rn-fp/sim/run_grid.R
# GRID_CELLS (an R index such as "1,40,72" or "49:72", passed on by
# dev/dyn-sim/aws/launch_spot.sh) restricts the run to those cells.
# =============================================================================
if (Sys.getenv("LOAD_ALL") == "1") devtools::load_all(quiet = TRUE) else library(crywolf)
reps <- as.integer(Sys.getenv("REPS", "100"))
workers <- as.integer(Sys.getenv("WORKERS", "8"))
adapt_delta <- as.numeric(Sys.getenv("ADAPT_DELTA", "0.99"))
job_dir <- Sys.getenv("JOBS", "dev/rn-fp/sim/jobs-grid")
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)

cells <- expand.grid(sigma_lam = c(0, 0.7), sigma_r = c(0, 1),
                     annotation = c("verified", "partial", "full"),
                     nocc = c(10, 20), nsite = c(25, 50, 100),
                     stringsAsFactors = FALSE)
cells$frac_confirmed <- c(verified = 0, partial = 0.2, full = 1)[cells$annotation]
cells$frac_verified <- c(verified = 0.3, partial = 0.3, full = 0)[cells$annotation]
pars <- c("sigma_det[1]", "sigma_lam[1]", "alpha[1]", "beta[1]", "beta[2]",
          "gamma0", "delta0", "r_ref", "p10_ref", "psi_ref", "psi_bar", "p_above_nmax")

tasks <- expand.grid(rep = seq_len(reps), cell = seq_len(nrow(cells)))
if (nzchar(Sys.getenv("GRID_CELLS"))) {
  keep <- eval(str2lang(paste0("c(", Sys.getenv("GRID_CELLS"), ")")))
  if (!length(keep) || !all(keep %in% seq_len(nrow(cells)))) {
    stop("GRID_CELLS must select cells among 1..", nrow(cells), call. = FALSE)
  }
  tasks <- tasks[tasks$cell %in% keep, ]
}
tasks$seed <- 800000 + 1000 * tasks$cell + tasks$rep
tasks$file <- file.path(job_dir, sprintf("cell%02d-rep%04d.rds", tasks$cell, tasks$rep))
tasks <- tasks[!file.exists(tasks$file), ]
if (nzchar(Sys.getenv("MAX_TASKS"))) tasks <- head(tasks, as.integer(Sys.getenv("MAX_TASKS")))
message(nrow(tasks), " tasks to run on ", workers, " workers")
invisible(crywolf:::.compile_model("fp_occu_rn"))

run_task <- function(cell, rep, seed, file) {
  cl <- cells[cell, ]
  s <- simulate_fp_rn(nsite = cl$nsite, nocc = cl$nocc, psi = 0.6, b1 = 0.5,
                      r = 0.3, p10 = 0.05, sigma_r = cl$sigma_r,
                      sigma_lam = cl$sigma_lam, rv = 0.5,
                      frac_verified = cl$frac_verified, nver = 3,
                      frac_confirmed = cl$frac_confirmed, seed = seed)
  out_dir <- tempfile("rn"); dir.create(out_dir)
  t0 <- Sys.time()
  m <- suppressWarnings(fp_occu_rn(
    s$data, occ = ~ x, random_det = TRUE, random_abund = TRUE, n_max = 25,
    chains = 4, parallel_chains = 1, threads_per_chain = 1, seed = seed,
    refresh = 0, show_messages = FALSE, output_dir = out_dir,
    adapt_delta = adapt_delta))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  d <- m$fit$diagnostic_summary(quiet = TRUE)
  sm <- as.data.frame(m$fit$summary(pars, "mean", "median",
                                    ~quantile(.x, c(0.05, 0.95)), "rhat",
                                    "ess_bulk", "ess_tail"))
  dr <- posterior::as_draws_matrix(m$fit$draws(c("sigma_det[1]", "sigma_lam[1]")))
  dr <- dr[round(seq(4, nrow(dr), length.out = 1000)), , drop = FALSE]
  res <- list(cell = cell, adapt_delta = adapt_delta, rep = rep, seed = seed, design = cl,
              truth = s$truth[c("b0", "b1", "a0", "r", "p10", "sigma_r",
                                "sigma_lam", "psi", "rv")],
              n_occupied = sum(s$truth$z), n_verified_sites = s$truth$n_verified_sites,
              summary = sm, secs = secs,
              rhat_max = max(m$fit$summary()$rhat, na.rm = TRUE),
              divergent = sum(d$num_divergent), treedepth = sum(d$num_max_treedepth),
              ebfmi = d$ebfmi, sigma_draws = unclass(dr))
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
