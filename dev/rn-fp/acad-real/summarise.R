# Summaries of the realistic ACAD grid (run_grid.R): per cell and prior on
# sigma_r, the share of fits that converged (largest R-hat < 1.01, no
# divergences), the median bias and the coverage of the 90% interval of
# sigma_r and p10, how often the lower 5% bound of sigma_r clears 0.1 (a false
# signal where sigma_r = 0, power where it is 1), and a Savage-Dickey Bayes
# factor BF01 for sigma_r = 0 (posterior density at 0 from a kernel estimate
# reflected at the boundary, over the prior density at 0; approximate).
# Fits refitted at adapt_delta 0.99 (run_rerun.R, folder RERUN_JOBS) replace
# the original ones; the column adapt_delta says which run each fit comes from.
#   JOBS=<folder> Rscript dev/rn-fp/acad-real/summarise.R
jobs <- Sys.getenv("JOBS", "dev/rn-fp/acad-real/jobs")
rerun_jobs <- Sys.getenv("RERUN_JOBS", "dev/rn-fp/acad-real/jobs-ad99")
files <- list.files(jobs, "\\.rds$", full.names = TRUE)
prior0 <- c(hn0.5 = 2 * dnorm(0, 0, 0.5), hn1 = 2 * dnorm(0, 0, 1),
            hn2 = 2 * dnorm(0, 0, 2), t3 = 2 * dt(0, 3))
dens0 <- function(x) 2 * stats::density(c(x, -x), from = 0, to = 0, n = 1)$y
rows <- lapply(files, function(f) {
  r <- readRDS(f)
  ad <- setNames(rep(0.95, length(r$fits)), names(r$fits))
  for (p in names(r$fits)) {
    g <- file.path(rerun_jobs, sprintf("cell%02d-rep%04d-%s.rds", r$cell, r$rep, p))
    if (file.exists(g)) { r$fits[[p]] <- readRDS(g)$fits[[p]]; ad[[p]] <- 0.99 }
  }
  do.call(rbind, lapply(r$fits, function(ft) {
    s <- ft$summary
    g <- function(v, col) s[s$variable == v, col]
    data.frame(cell = r$cell, rep = r$rep, psi = r$design$psi, p10 = r$design$p10,
               sigma_r = r$design$sigma_r, prior = ft$prior, adapt_delta = ad[[ft$prior]],
               secs = ft$secs, divergent = ft$divergent, rhat_max = ft$rhat_max,
               ok = ft$rhat_max < 1.01 && ft$divergent == 0,
               sr_med = g("sigma_det[1]", "median"), sr_lo = g("sigma_det[1]", "5%"),
               sr_hi = g("sigma_det[1]", "95%"),
               p10_med = g("p10_ref", "median"), p10_lo = g("p10_ref", "5%"),
               p10_hi = g("p10_ref", "95%"),
               bf01 = dens0(ft$sigma_r_draws) / prior0[[ft$prior]])
  }))
})
d <- do.call(rbind, rows)
cover <- function(lo, hi, t) mean(lo <= t & t <= hi)
out <- do.call(rbind, lapply(split(d, list(d$cell, d$prior), drop = TRUE), function(x) data.frame(
  psi = x$psi[1], p10 = x$p10[1], sigma_r = x$sigma_r[1], prior = x$prior[1],
  n = nrow(x), converged = mean(x$ok), div_any = mean(x$divergent > 0), minutes = round(median(x$secs) / 60, 2),
  sr_bias = median(x$sr_med - x$sigma_r), sr_cover = cover(x$sr_lo, x$sr_hi, x$sigma_r),
  sr_lo_gt_0.1 = mean(x$sr_lo > 0.1),
  bf01_med = median(x$bf01), bf01_gt3 = mean(x$bf01 > 3), bf01_lt1_3 = mean(x$bf01 < 1/3),
  p10_bias = median(x$p10_med - x$p10), p10_cover = cover(x$p10_lo, x$p10_hi, x$p10))))
out <- out[order(-out$psi, out$p10, out$sigma_r, out$prior), ]
print(out, digits = 2, row.names = FALSE)
saveRDS(list(fits = d, cells = out), "dev/rn-fp/acad-real/grid.rds")
