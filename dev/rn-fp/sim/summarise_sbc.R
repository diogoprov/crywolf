# Summaries of the SBC (run_sbc.R): per design and parameter, a chi-squared
# test of uniformity of the ranks (0..100) in 20 bins, the share of fits with
# divergences or R-hat >= 1.01, and rank histograms (sbc-ranks.png). The test
# is repeated without the fits that had divergences, to see whether they drive
# any departure.
#   JOBS=<folder> Rscript dev/rn-fp/sim/summarise_sbc.R
jobs <- Sys.getenv("JOBS", "dev/rn-fp/sim/jobs-sbc")
res <- lapply(list.files(jobs, "\\.rds$", full.names = TRUE), readRDS)
d <- do.call(rbind, lapply(res, function(r) data.frame(
  design = r$design, rep = r$rep, parameter = names(r$ranks), rank = r$ranks,
  divergent = r$divergent, rhat_max = r$rhat_max,
  ess_min = min(r$summary$ess_bulk), p_above_nmax = r$p_above_nmax,
  n_max_true = r$n_max_true, minutes = r$secs / 60)))
unif <- function(rank) {
  bins <- tabulate(pmin(rank %/% 5, 19) + 1, 20)  # 0..100 -> 20 bins; 100 joins the last
  # expected counts: 5 rank values per bin, 6 in the last
  p <- c(rep(5, 19), 6) / 101
  suppressWarnings(stats::chisq.test(bins, p = p)$p.value)
}
out <- do.call(rbind, lapply(split(d, list(d$design, d$parameter), drop = TRUE), function(x) {
  ok <- x$divergent == 0
  data.frame(design = x$design[1], parameter = x$parameter[1], n = nrow(x),
             p_uniform = unif(x$rank), n_nodiv = sum(ok), p_uniform_nodiv = unif(x$rank[ok]),
             mean_rank = mean(x$rank))
}))
fits <- unique(d[, c("design", "rep", "divergent", "rhat_max", "ess_min",
                     "p_above_nmax", "n_max_true", "minutes")])
diag <- do.call(rbind, lapply(split(fits, fits$design), function(x) data.frame(
  design = x$design[1], n = nrow(x), div_any = mean(x$divergent > 0),
  rhat_ge_1.01 = mean(x$rhat_max >= 1.01), ess_min_med = stats::median(x$ess_min),
  p_above_nmax_gt_0.001 = mean(x$p_above_nmax > 0.001),
  minutes_med = stats::median(x$minutes))))
print(diag, digits = 3, row.names = FALSE)
print(out[order(out$design, out$parameter), ], digits = 3, row.names = FALSE)
saveRDS(list(ranks = d, uniformity = out, diagnostics = diag), "dev/rn-fp/sim/sbc.rds")

png("dev/rn-fp/sim/sbc-ranks.png", width = 1600, height = 700, res = 120)
pars <- unique(d$parameter)
op <- par(mfrow = c(2, length(pars)), mar = c(3, 3, 2, 0.5), mgp = c(1.8, 0.5, 0))
for (dg in c("A", "B")) for (p in pars) {
  r <- d$rank[d$design == dg & d$parameter == p]
  if (!length(r)) { plot.new(); next }
  h <- graphics::hist(r, breaks = seq(-0.5, 100.5, length.out = 21), plot = FALSE)
  n <- length(r)
  plot(h, main = paste0(dg, ": ", p), xlab = "rank", ylab = "count", col = "grey80",
       border = "white", ylim = c(0, max(h$counts, stats::qbinom(0.995, n, 1/20))))
  abline(h = stats::qbinom(c(0.005, 0.5, 0.995), n, 1/20), lty = c(2, 1, 2), col = "grey30")
}
par(op); invisible(dev.off())
