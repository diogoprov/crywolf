# =============================================================================
# dev/ms-sim/summarise.R
# Reads the job files of run_sims.R (dev/ms-sim/jobs/) and writes
# vignettes/multispecies-results.rds (read by vignettes/multispecies.Rmd) and
# the tables of dev/ms-sim/results.md. Species-level values ("[sp3]") are
# pooled over species as "[sp]"; in the rare and unverified scenarios the
# first species is also kept apart ("[sp1]"), to compare the multispecies fit
# with fp_occu() alone. Fits count as converged with R-hat < 1.01 for every
# parameter with a true value and no divergent transitions. Cells that have
# not ended yet are left out. Run from the package root:
#   Rscript dev/ms-sim/summarise.R
# =============================================================================

job_dir <- file.path("dev", "ms-sim", "jobs")
read_part <- function(part) {
  f <- list.files(job_dir, pattern = paste0("^", part, "-"), full.names = TRUE)
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, readRDS))
}
recovery <- read_part("recovery")
sbc <- read_part("sbc")

# one row per parameter (and grouping), with bias, RMSE, 95% coverage and width
metrics <- function(d, by) {
  key <- interaction(d[c(by, "param")], drop = TRUE, sep = "|")
  out <- lapply(split(d, key), function(x) {
    data.frame(x[1, c(by, "param")], n = nrow(x),
               bias = mean(x$mean - x$truth), rmse = sqrt(mean((x$mean - x$truth)^2)),
               coverage = mean(x$lower <= x$truth & x$truth <= x$upper),
               width = mean(x$upper - x$lower), row.names = NULL)
  })
  out <- do.call(rbind, out)
  out[do.call(order, out[c(by, "param")]), ]
}
fit_keys <- c("rep", "cell", "model")
fit_ok <- function(d) {
  bad <- !d$ok | (!is.na(d$rhat) & d$rhat >= 1.01) | (!is.na(d$divergences) & d$divergences > 0)
  !ave(bad, interaction(d[fit_keys], drop = TRUE), FUN = any)
}

res <- list()
if (!is.null(recovery)) {
  recovery$converged <- fit_ok(recovery)
  recovery$param <- sub("\\[sp\\d+\\]$", "[sp]", recovery$variable)
  # one row per fit, with its largest R-hat and smallest bulk ESS
  fits <- unique(recovery[c(fit_keys, "nspecies", "scenario", "converged", "seconds",
                            "divergences")])
  worst <- stats::aggregate(cbind(rhat = recovery$rhat, ess = -recovery$ess_bulk),
                            recovery[fit_keys], max)
  fits <- merge(fits, worst, by = fit_keys)
  res$fits <- do.call(rbind, lapply(split(fits, interaction(fits[c("nspecies", "scenario",
                                                                    "model")], drop = TRUE)),
    function(x) data.frame(x[1, c("nspecies", "scenario", "model")], n = nrow(x),
                           converged = mean(x$converged),
                           with_divergences = mean(x$divergences > 0, na.rm = TRUE),
                           max_rhat = max(x$rhat), min_ess = min(-x$ess),
                           median_hours = stats::median(x$seconds) / 3600, row.names = NULL)))
  res$fits <- res$fits[order(res$fits$nspecies, res$fits$scenario, res$fits$model), ]
  ok <- recovery[recovery$ok, ]
  # community and pooled species values, multispecies fits only
  res$recovery <- metrics(ok[ok$model == "ms", ], c("nspecies", "scenario"))
  # the first species of the rare and unverified scenarios, multispecies vs alone
  sp1 <- ok[ok$scenario %in% c("rare", "unverified") & grepl("\\[sp1\\]$", ok$variable), ]
  sp1$param <- sp1$variable
  res$sp1 <- metrics(sp1, c("nspecies", "scenario", "model"))
  # fs and mixed models on the fs and mixed datasets
  re <- ok[ok$scenario %in% c("fs", "mixed"), ]
  res$fs_mixed <- metrics(re, c("scenario", "model"))
}
if (!is.null(sbc)) {
  ok <- sbc[sbc$ok, ]
  # 100 possible ranks (0 to 99) in 10 bins of 10; chi-squared test of uniformity
  res$sbc_ranks <- ok[c("rep", "variable", "truth", "rank", "max_rhat", "divergences")]
  res$sbc <- do.call(rbind, lapply(split(ok, ok$variable), function(x) {
    counts <- tabulate(floor(x$rank / 10) + 1, nbins = 10)
    data.frame(variable = x$variable[1], n = nrow(x),
               chisq_p = stats::chisq.test(counts)$p.value,
               mean_rank = mean(x$rank), row.names = NULL)
  }))
  res$sbc_fits <- c(n = length(unique(sbc$rep)), failed = length(unique(sbc$rep[!sbc$ok])),
                    with_divergences = length(unique(ok$rep[ok$divergences > 0])))
}
saveRDS(res, file.path("vignettes", "multispecies-results.rds"))

md <- function(x, digits = 3) {
  x[] <- lapply(x, function(v) if (is.numeric(v)) signif(v, digits) else v)
  c(paste("|", paste(names(x), collapse = " | "), "|"),
    paste("|", paste(rep("---", ncol(x)), collapse = " | "), "|"),
    apply(x, 1, function(r) paste("|", paste(r, collapse = " | "), "|")))
}
out <- c("# Simulations of the multispecies model: results", "",
         "Written by `dev/ms-sim/summarise.R`.", "")
if (!is.null(res$fits)) {
  community <- !grepl("\\[sp\\]$", res$recovery$param)
  out <- c(out, "## Fits", "", md(res$fits), "",
           "## Recovery: community parameters", "", md(res$recovery[community, ]), "",
           "## Recovery: species values, pooled over species", "",
           md(res$recovery[!community, ]), "",
           "## First species: multispecies (ms) vs alone (single)", "", md(res$sp1), "",
           "## fs and mixed", "", md(res$fs_mixed), "")
}
if (!is.null(res$sbc)) {
  out <- c(out, "## Simulation-based calibration", "",
           paste0("Fits: ", res$sbc_fits[["n"]], ", failed: ", res$sbc_fits[["failed"]],
                  ", with divergences: ", res$sbc_fits[["with_divergences"]], "."), "",
           md(res$sbc), "")
}
writeLines(out, file.path("dev", "ms-sim", "results.md"))
