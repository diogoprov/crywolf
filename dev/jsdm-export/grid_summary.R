# Summary of dev/jsdm-export/grid/ (grid.R): convergence of the fp_occu_ms()
# fits, then, per cell, the mean over replicates of the recovery of LL' by
# gllvm. Writes grid_summary.rds.
#   Rscript dev/jsdm-export/grid_summary.R
files <- list.files("dev/jsdm-export/grid", full.names = TRUE)
res <- lapply(files, readRDS)
diag <- do.call(rbind, lapply(res, function(r) data.frame(
  S = r$S, load_sd = r$load_sd, seed = r$seed, minutes = r$mins_fp_occu_ms,
  draws_ok = r$draws_ok, true_mean_abs = r$true_mean_abs,
  divergences = sum(r$diagnostics$divergences), max_rhat = r$diagnostics$max_rhat,
  min_ess_bulk = r$diagnostics$min_ess_bulk)))
long <- do.call(rbind, lapply(res, function(r) data.frame(
  S = r$S, load_sd = r$load_sd, seed = r$seed, input = rownames(r$scores),
  cor = r$scores[, "cor"], mean_abs = r$scores[, "mean_abs"], row.names = NULL)))
agg <- aggregate(cbind(cor, mean_abs) ~ input + load_sd + S, long, mean)
agg_sd <- aggregate(cor ~ input + load_sd + S, long, stats::sd)
agg$cor_sd <- agg_sd$cor
# how often the draws beat the naive states, per dataset
w <- reshape(long[c("S", "load_sd", "seed", "input", "cor")], idvar = c("S", "load_sd", "seed"),
             timevar = "input", direction = "wide")
w$draws_beat_naive <- w$cor.draws_mean > w$cor.naive
beat <- aggregate(draws_beat_naive ~ load_sd + S, w, sum)
true_abs <- aggregate(true_mean_abs ~ load_sd + S, diag, mean)
saveRDS(list(diagnostics = diag, by_dataset = long, by_cell = agg, beat = beat,
             true_mean_abs = true_abs), "dev/jsdm-export/grid_summary.rds")
cat("Convergence of fp_occu_ms():\n")
print(summary(diag[c("minutes", "divergences", "max_rhat", "min_ess_bulk", "draws_ok")]))
cat("\nRecovery of LL' (mean over 5 replicates):\n")
agg$input <- factor(agg$input, c("oracle", "naive", "draws_mean", "draw_single_median"))
print(agg[order(agg$S, agg$load_sd, agg$input), ], digits = 2, row.names = FALSE)
cat("\nDatasets where the mean over draws beat the naive states (of 5):\n")
print(beat, row.names = FALSE)
print(true_abs, digits = 2, row.names = FALSE)
