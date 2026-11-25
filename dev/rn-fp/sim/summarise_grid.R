# Summaries of the general grid (run_grid.R), per cell: the share of fits that
# converged (largest R-hat < 1.01, no divergences), median bias and coverage of
# the 90% interval of sigma_r, sigma_lam, p10, r and psi at x = 0, how often
# the lower 5% bound of sigma_r clears 0.1 (a false signal where sigma_r = 0,
# power where it is 1), and a Savage-Dickey Bayes factor BF01 for sigma_r = 0
# (posterior density at 0 from a kernel estimate reflected at the boundary,
# over the half-Normal(0, 1) prior density at 0; approximate).
#   JOBS=<folder> Rscript dev/rn-fp/sim/summarise_grid.R
jobs <- Sys.getenv("JOBS", "dev/rn-fp/sim/jobs-grid")
files <- list.files(jobs, "\\.rds$", full.names = TRUE)
dens0 <- function(x) 2 * stats::density(c(x, -x), from = 0, to = 0, n = 1)$y
prior0 <- 2 * stats::dnorm(0)
d <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  s <- r$summary
  g <- function(v, col) s[s$variable == v, col]
  q <- function(v, nm, truth) {
    out <- data.frame(g(v, "median") - truth, g(v, "5%") <= truth & truth <= g(v, "95%"))
    names(out) <- paste0(nm, c("_err", "_cov")); out
  }
  cbind(data.frame(cell = r$cell, rep = r$rep, nsite = r$design$nsite, nocc = r$design$nocc,
                   annotation = r$design$annotation, sigma_r = r$design$sigma_r,
                   sigma_lam = r$design$sigma_lam, minutes = r$secs / 60,
                   divergent = r$divergent, rhat_max = r$rhat_max,
                   ok = r$rhat_max < 1.01 && r$divergent == 0,
                   sr_lo = g("sigma_det[1]", "5%"),
                   bf01 = dens0(r$sigma_draws[, "sigma_det[1]"]) / prior0),
        q("sigma_det[1]", "sr", r$truth$sigma_r), q("sigma_lam[1]", "sl", r$truth$sigma_lam),
        q("p10_ref", "p10", r$truth$p10), q("r_ref", "r", r$truth$r),
        q("psi_ref", "psi", r$truth$psi))
}))
cells <- do.call(rbind, lapply(split(d, d$cell), function(x) data.frame(
  cell = x$cell[1], nsite = x$nsite[1], nocc = x$nocc[1], annotation = x$annotation[1],
  sigma_r = x$sigma_r[1], sigma_lam = x$sigma_lam[1], n = nrow(x),
  converged = mean(x$ok), div_any = mean(x$divergent > 0),
  minutes = stats::median(x$minutes),
  sr_bias = stats::median(x$sr_err), sr_cover = mean(x$sr_cov),
  sr_lo_gt_0.1 = mean(x$sr_lo > 0.1), bf01_gt3 = mean(x$bf01 > 3),
  bf01_lt1_3 = mean(x$bf01 < 1/3),
  sl_bias = stats::median(x$sl_err), sl_cover = mean(x$sl_cov),
  p10_bias = stats::median(x$p10_err), p10_cover = mean(x$p10_cov),
  r_bias = stats::median(x$r_err), r_cover = mean(x$r_cov),
  psi_bias = stats::median(x$psi_err), psi_cover = mean(x$psi_cov))))
print(cells, digits = 2, row.names = FALSE)
saveRDS(list(fits = d, cells = cells), "dev/rn-fp/sim/grid.rds")
