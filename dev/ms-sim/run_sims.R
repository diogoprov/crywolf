# =============================================================================
# dev/ms-sim/run_sims.R
# PR 5 of the multispecies model (dev/proposal-mso.md, section 5):
# simulations for fp_occu_ms(), all from simulate_fp_ms() (100 sites except
# in the grid, 10 automated occasions per site, 30% of the sites verified with 3 occasions;
# occupancy ~ x and detection ~ bio by species, constant false-positive rate
# with a gap by species).
#
#   recovery  community of 5 or 20 species x 3 scenarios, 50 replicates
#             per cell: "base"; "rare", the first species with occupancy
#             0.05 at the mean of x; "unverified", the first species without
#             verified occasions. 50 species, "base" only, 20 replicates, to
#             show that the model scales. In the rare and unverified
#             scenarios the first species is also fitted alone with
#             fp_occu(), to measure how much the community lends it.
#             Pilot (cell 10): the datasets of the 50 species cell refitted
#             with centered = TRUE (model "ms_centered"), to see whether
#             the centred parametrization of the detection intercepts and
#             gaps mixes better than the default. Pilot (cell 11): the same
#             datasets with centered = TRUE and centered_occ = TRUE (model
#             "ms_centered_occ"), to see whether the centred occupancy
#             coefficients remove the failures left by cell 10.
#             Community: psi = 0.4 (sd 1 on the logit), b1 = -0.5 (sd 0.5),
#             p11 = 0.5 (sd 0.5), p10 = 0.05 (sd of log gap 0.3),
#             a_bio = 0.5, r11 = 0.6.
#   fs vs mixed  10 species, the 100 sites in 25 groups of 4 (as physical
#             sites with several months), 50 replicates per scenario:
#             "fs", random effects on logit(p10) by group x species
#             (sd_fp_group = 0.7); "mixed", these plus effects by site x
#             species (sd_fp_site = 0.5). Both models are fitted to both
#             datasets: "fs" (fp_group) and "mixed" (fp_group + random_fp).
#   grid      when do the centred parametrizations help? Base scenario
#             with 5 and 20 species at 30, 100 and 300 sites (10
#             occasions), 20 replicates, each dataset fitted four ways:
#             the default ("ms"), centered = TRUE ("ms_centered"),
#             centered_occ = TRUE ("ms_occ") and both ("ms_centered_occ").
#             At 100 sites the datasets are those of the recovery cells
#             (same seeds), whose default fits are not repeated; with 50
#             species at 100 sites only "ms_occ" is new (cells 7, 10, 11).
#   sbc       simulation-based calibration: the community parameters drawn
#             from fp_priors() (means from the priors of the intercepts,
#             slopes and gap; standard deviations among species from
#             half-Normal(0, species_sd)), 10 species, 80 sites, 8
#             occasions; the rank of each true value among 99 thinned
#             posterior draws. 200 replicates.
#
# Every fit: 4 chains of 1000 + 1000 iterations, adapt_delta = 0.95, one
# core (threads_per_chain = 1), random_recorder = FALSE. Seeds: 91001 on
# + 1000 x cell + replicate (recovery), 95001 on (sbc), 200001 on + 1000 x
# dataset (grid, at 30 and 300 sites), one per task.
#
# Each task is saved in dev/ms-sim/jobs/ as soon as it ends and skipped on
# a rerun. Run from the package root, with crywolf installed from this
# checkout:
#   WORKERS=16 Rscript dev/ms-sim/run_sims.R
# REPS sets the replicates of every cell (REPS = 1: one per cell, for a test); ITER
# sets the iterations (a pilot with fewer gives the time per iteration); PARTS
# (default "recovery,sbc,grid") selects the parts; CELLS (e.g. "2,5") keeps some
# recovery cells, for timing, and GRID_CELLS some cells of the grid. Do not edit
# this file while it runs.
# =============================================================================

library(crywolf)

reps <- as.integer(Sys.getenv("REPS", "0"))   # 0: the replicates above
workers <- as.integer(Sys.getenv("WORKERS", "16"))
parts <- strsplit(Sys.getenv("PARTS", "recovery,sbc,grid"), ",")[[1]]
cells_env <- Sys.getenv("CELLS", "")
grid_env <- Sys.getenv("GRID_CELLS", "")
job_dir <- Sys.getenv("JOB_DIR", file.path("dev", "ms-sim", "jobs"))   # pilots elsewhere
dir.create(job_dir, showWarnings = FALSE, recursive = TRUE)
n_rep <- function(full) if (reps > 0) reps else full

recovery_cells <- rbind(
  expand.grid(nspecies = c(5, 20), scenario = c("base", "rare", "unverified"),
              reps = 50, stringsAsFactors = FALSE),
  data.frame(nspecies = 50, scenario = "base", reps = 20),
  data.frame(nspecies = 10, scenario = c("fs", "mixed"), reps = 50),
  data.frame(nspecies = 50, scenario = "base", reps = 20),
  data.frame(nspecies = 50, scenario = "base", reps = 20)
)
recovery_cells$centered <- seq_len(nrow(recovery_cells)) %in% c(10, 11)
recovery_cells$centered_occ <- seq_len(nrow(recovery_cells)) == 11
# the pilots reuse the seeds, so the datasets, of cell 7
recovery_cells$seed_cell <- ifelse(recovery_cells$centered, 7, seq_len(nrow(recovery_cells)))
# grid: one fit per task; the datasets are shared by the four
# parametrizations through their seeds
grid_cells <- expand.grid(model = c("ms", "ms_centered", "ms_occ", "ms_centered_occ"),
                          nsite = c(30, 100, 300), nspecies = c(5, 20),
                          stringsAsFactors = FALSE)
grid_cells <- rbind(grid_cells[!(grid_cells$nsite == 100 & grid_cells$model == "ms"), ],
                    data.frame(model = "ms_occ", nsite = 100, nspecies = 50))
rownames(grid_cells) <- NULL
# seeds of the datasets: the recovery cells 1, 2 and 7 at 100 sites
grid_cells$seed_base <- ifelse(grid_cells$nsite == 100,
                               91000 + 1000 * c(`5` = 1, `20` = 2, `50` = 7)[
                                 as.character(grid_cells$nspecies)],
                               200000 + 1000 * as.integer(factor(paste(grid_cells$nspecies,
                                                                       grid_cells$nsite))))
keep_cells <- if (nzchar(cells_env)) as.integer(strsplit(cells_env, ",")[[1]]) else
  seq_len(nrow(recovery_cells))

tasks <- rbind(
  if ("recovery" %in% parts) {
    t <- do.call(rbind, lapply(seq_len(nrow(recovery_cells)), function(i)
      data.frame(rep = seq_len(n_rep(recovery_cells$reps[i])), cell = i)))
    t <- data.frame(part = "recovery", t,
                    seed = 91000 + 1000 * recovery_cells$seed_cell[t$cell] + t$rep)
    t[t$cell %in% keep_cells, ]
  },
  if ("grid" %in% parts) {
    t <- do.call(rbind, lapply(seq_len(nrow(grid_cells)), function(i)
      data.frame(rep = seq_len(n_rep(20)), cell = i)))
    t <- data.frame(part = "grid", t, seed = grid_cells$seed_base[t$cell] + t$rep)
    if (nzchar(grid_env)) t[t$cell %in% as.integer(strsplit(grid_env, ",")[[1]]), ] else t
  },
  if ("sbc" %in% parts) {
    t <- data.frame(rep = seq_len(n_rep(200)), cell = 1L)
    data.frame(part = "sbc", t, seed = 95000 + seq_len(nrow(t)))
  }
)
# longest tasks first (most species), so that the last ones to end are short
size <- ifelse(tasks$part == "recovery",
               100 * recovery_cells$nspecies[pmin(tasks$cell, nrow(recovery_cells))],
               ifelse(tasks$part == "grid",
                      grid_cells$nspecies[pmin(tasks$cell, nrow(grid_cells))] *
                        grid_cells$nsite[pmin(tasks$cell, nrow(grid_cells))],
                      800))
tasks <- tasks[order(-size), ]

# ITER (default 1000) sets warmup and sampling iterations, for a timing pilot
iter <- as.integer(Sys.getenv("ITER", "1000"))
fit_args <- list(random_recorder = FALSE, chains = 4, parallel_chains = 1,
                 threads_per_chain = 1, iter_warmup = iter, iter_sampling = iter,
                 adapt_delta = 0.95, refresh = 0, show_messages = FALSE)

# summary rows of the parameters with a true value
summary_rows <- function(m, truth) {
  s <- summary(m)$table
  s <- s[s$variable %in% names(truth), c("variable", "mean", "lower", "upper",
                                         "rhat", "ess_bulk")]
  s$truth <- unname(truth[s$variable])
  as.data.frame(s)
}

# one model fitted to one dataset -> rows, or one row with the error
fit_rows <- function(model, expr) {
  t0 <- Sys.time()
  res <- tryCatch(expr, error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (inherits(res, "error")) {
    return(data.frame(model = model, ok = FALSE, msg = conditionMessage(res),
                      variable = NA_character_, mean = NA_real_, lower = NA_real_,
                      upper = NA_real_, rhat = NA_real_, ess_bulk = NA_real_,
                      truth = NA_real_, divergences = NA_integer_,
                      max_treedepth = NA_integer_, seconds = secs))
  }
  data.frame(model = model, ok = TRUE, msg = NA_character_, res$rows,
             divergences = res$m$diagnostics$divergences,
             max_treedepth = res$m$diagnostics$max_treedepth, seconds = secs)
}

run_recovery <- function(cell, seed) {
  # the job rows keep the columns of the first run, so that they bind
  centered <- recovery_cells$centered[cell]
  centered_occ <- recovery_cells$centered_occ[cell]
  cfg <- recovery_cells[cell, c("nspecies", "scenario", "reps")]
  re <- cfg$scenario %in% c("fs", "mixed")
  s <- simulate_fp_ms(nspecies = cfg$nspecies, nsite = 100, nocc = 10, a_bio = 0.5,
                      psi_rare = if (cfg$scenario == "rare") 0.05,
                      ngroup = if (re) 25, sd_fp_group = if (re) 0.7 else 0,
                      sd_fp_site = if (cfg$scenario == "mixed") 0.5 else 0, seed = seed)
  d <- s$data
  if (cfg$scenario == "unverified") {
    d <- suppressWarnings(fp_data(d$auto, d$sites, d$verified[d$verified$species != "sp1", ],
                                  species = "species"))
  }
  truth <- s$truth$pars
  fit_ms <- function(model, ...) fit_rows(model, {
    m <- suppressWarnings(do.call(fp_occu_ms, c(list(d, occ = ~ x, det = ~ bio, fp = ~ 1,
                                                     seed = seed, ...), fit_args)))
    list(m = m, rows = summary_rows(m, truth))
  })
  if (re) {
    rows <- rbind(fit_ms("fs", fp_group = "group", centered = FALSE, centered_occ = FALSE),
                  fit_ms("mixed", fp_group = "group", random_fp = TRUE, centered = FALSE,
                         centered_occ = FALSE))
    return(cbind(cfg, n_verified_sites = length(unique(d$verified$site)), rows))
  }
  if (centered) {
    rows <- if (centered_occ) fit_ms("ms_centered_occ", centered = TRUE, centered_occ = TRUE) else
      fit_ms("ms_centered", centered = TRUE, centered_occ = FALSE)
    return(cbind(cfg, n_verified_sites = length(unique(d$verified$site)), rows))
  }
  rows <- fit_ms("ms", centered = FALSE, centered_occ = FALSE)
  if (cfg$scenario != "base") {
    # the first species alone, its values named as in the multispecies output
    one <- suppressWarnings(fp_data(d$auto[d$auto$species == "sp1", ], d$sites,
                                    d$verified[d$verified$species == "sp1", ]))
    t1 <- truth[grepl("\\[sp1\\]$", names(truth))]
    names(t1) <- sub("\\[sp1\\]$", "", names(t1))
    single <- fit_rows("single", {
      m <- suppressWarnings(do.call(fp_occu, c(list(one, occ = ~ x, det = ~ bio, fp = ~ 1,
                                                    seed = seed), fit_args)))
      r <- summary_rows(m, t1)
      r$variable <- paste0(r$variable, "[sp1]")
      list(m = m, rows = r)
    })
    rows <- rbind(rows, single)
  }
  cbind(cfg, n_verified_sites = length(unique(d$verified$site)), rows)
}

# one parametrization of the grid fitted to one base dataset
run_grid <- function(cell, seed) {
  g <- grid_cells[cell, ]
  s <- simulate_fp_ms(nspecies = g$nspecies, nsite = g$nsite, nocc = 10, a_bio = 0.5,
                      seed = seed)
  args <- switch(g$model, ms = list(centered = FALSE, centered_occ = FALSE),
                 ms_centered = list(centered = TRUE, centered_occ = FALSE),
                 ms_occ = list(centered = FALSE, centered_occ = TRUE),
                 ms_centered_occ = list(centered = TRUE, centered_occ = TRUE))
  rows <- fit_rows(g$model, {
    m <- suppressWarnings(do.call(fp_occu_ms, c(list(s$data, occ = ~ x, det = ~ bio,
                                                     fp = ~ 1, seed = seed),
                                                args, fit_args)))
    list(m = m, rows = summary_rows(m, s$truth$pars))
  })
  cbind(nspecies = g$nspecies, nsite = g$nsite, rows)
}

# community parameters drawn from fp_priors(), as in the model block of
# fp_occu_ms.stan: means Normal(0, intercept_sd) and Normal(0, slope_sd),
# median gap Normal(gap_mean, gap_sd) truncated at 0, standard deviations
# half-Normal(0, species_sd), verified intercept Normal(0, intercept_sd)
run_sbc <- function(seed) {
  pr <- fp_priors()
  set.seed(seed)
  half <- function() abs(stats::rnorm(1, 0, pr$species_sd))
  repeat {
    gap <- stats::rnorm(1, pr$gap_mean, pr$gap_sd)
    if (gap > 0) break
  }
  b0 <- stats::rnorm(1, 0, pr$intercept_sd)
  alpha <- stats::rnorm(1, 0, pr$intercept_sd)
  s <- simulate_fp_ms(nspecies = 10, nsite = 80, nocc = 8,
                      psi = stats::plogis(b0), sd_psi = half(),
                      b1 = stats::rnorm(1, 0, pr$slope_sd), sd_b1 = half(),
                      p11 = stats::plogis(alpha), sd_p11 = half(),
                      p10 = stats::plogis(alpha - gap), sd_gap = half(),
                      a_bio = stats::rnorm(1, 0, pr$slope_sd),
                      r11 = stats::plogis(stats::rnorm(1, 0, pr$intercept_sd)),
                      seed = seed)
  keep <- c("occ_(Intercept)", "occ_x", "sd_occ_(Intercept)", "sd_occ_x",
            "det_(Intercept)", "det_bio", "sd_det_(Intercept)", "fp_gap",
            "sd_log_fp_gap", "verified_(Intercept)", "occ_(Intercept)[sp1]",
            "det_(Intercept)[sp1]", "p10_ref[sp1]")
  truth <- s$truth$pars[keep]
  t0 <- Sys.time()
  res <- tryCatch({
    m <- suppressWarnings(do.call(fp_occu_ms, c(list(s$data, occ = ~ x, det = ~ bio,
                                                     fp = ~ 1, seed = seed, centered = FALSE,
                                                     centered_occ = FALSE), fit_args)))
    dr <- posterior::as_draws_matrix(as_draws_df(m))
    thin <- dr[round(seq(1, posterior::ndraws(dr), length.out = 99)), names(truth)]
    list(ranks = colSums(sweep(unclass(thin), 2, truth, "<")), diag = m$diagnostics)
  }, error = function(e) e)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ok <- !inherits(res, "error")
  data.frame(variable = names(truth), truth = unname(truth),
             rank = if (ok) unname(res$ranks[names(truth)]) else NA_real_,
             ok = ok, msg = if (ok) NA_character_ else conditionMessage(res),
             max_rhat = if (ok) res$diag$max_rhat else NA_real_,
             divergences = if (ok) res$diag$divergences else NA_integer_,
             seconds = secs)
}

run_one <- function(part, rep, cell, seed) {
  f <- file.path(job_dir, sprintf("%s-cell%02d-rep%04d.rds", part, cell, rep))
  if (file.exists(f)) return(invisible(NULL))
  out <- switch(part, recovery = run_recovery(cell, seed), grid = run_grid(cell, seed),
                sbc = run_sbc(seed))
  saveRDS(cbind(part = part, rep = rep, cell = cell, seed = seed, out), f)
  invisible(NULL)
}

crywolf::fp_compile()                          # once, before the workers
invisible(crywolf:::.compile_model("fp_occu_ms"))
future::plan(future::multisession, workers = workers)
t0 <- Sys.time()
furrr::future_pwalk(tasks, run_one,
                    .options = furrr::furrr_options(seed = TRUE, chunk_size = 1))
message("done: ", nrow(tasks), " tasks in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "hours")), 2), " h")
