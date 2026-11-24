# Summaries of the step 1 grid (run_grid.R): per cell, the share of fits that
# converged (largest R-hat < 1.01, no divergences), the bias and the coverage
# of the 90% intervals of sigma_r and sigma_lam, and, where the true value is
# 0, how often the lower 5% bound of sigma_r clears 0.1 (a false signal);
# where it is positive, how often it does (power). Run from the package root:
#   Rscript dev/rn-fp/summarise_grid.R
files <- list.files(file.path("dev", "rn-fp", "jobs"), "\\.rds$", full.names = TRUE)
res <- lapply(files, readRDS)
row <- function(r) {
  s <- r$summary
  get <- function(v, col) s[s$variable == v, col]
  data.frame(cell = r$cell, rep = r$rep, r$design, secs = r$secs,
             divergent = r$divergent, rhat_max = max(s$rhat),
             sr_mean = get("sigma_r_out", "mean"), sr_lo = get("sigma_r_out", "5%"),
             sr_hi = get("sigma_r_out", "95%"),
             sl_mean = get("sigma_lam_out", "mean"), sl_lo = get("sigma_lam_out", "5%"),
             sl_hi = get("sigma_lam_out", "95%"))
}
d <- do.call(rbind, lapply(res, row))
d$ok <- d$rhat_max < 1.01 & d$divergent == 0
cover <- function(lo, hi, truth) mean(lo <= truth & truth <= hi)
out <- do.call(rbind, lapply(split(d, d$cell), function(x) data.frame(
  cell = x$cell[1], sigma_r = x$sigma_r[1], sigma_lam = x$sigma_lam[1],
  frac_verified = x$frac_verified[1], n = nrow(x), converged = mean(x$ok),
  div_any = mean(x$divergent > 0), minutes = round(median(x$secs) / 60, 1),
  sr_bias = mean(x$sr_mean - x$sigma_r), sr_cover = cover(x$sr_lo, x$sr_hi, x$sigma_r),
  sr_lo_gt_0.1 = mean(x$sr_lo > 0.1),
  sl_bias = mean(x$sl_mean - x$sigma_lam), sl_cover = cover(x$sl_lo, x$sl_hi, x$sigma_lam))))
print(out, digits = 2, row.names = FALSE)
saveRDS(list(fits = d, cells = out), file.path("dev", "rn-fp", "grid.rds"))
