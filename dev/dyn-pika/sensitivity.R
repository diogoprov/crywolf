# =============================================================================
# dev/dyn-pika/sensitivity.R
# Sensitivity of the fit of run.R to the definition of the verified method.
# Goldman et al. (2023) class fresh scat as unambiguous when a fresh haypile was
# seen on the same survey, so in `pika` every survey with both has FS_U > 0 and
# the verified occasion depends on the automated haypile occasion. The notes
# behind the other unambiguous scat are consistent only from 2019 on.
#
#   base:     verified = sighting or FS_U > 0, as in run.R
#   no_rule:  verified = sighting, or FS_U > 0 on a survey without fresh haypile
#   sighting: verified = sighting; all fresh scat is automated
#
# Scat that is not verified stays in the automated scat occasion. The model and
# sampler settings are those of the default fit of run.R.
# Run from the package root: Rscript dev/dyn-pika/sensitivity.R
# =============================================================================

library(crywolf)
dir <- file.path("dev", "dyn-pika")
u <- ifelse(is.na(pika$FS_U), 0, pika$FS_U) > 0
year <- factor(pika$Year, levels = as.character(2010:2021))
plots <- unique(pika[c("plot", "Park")])
sites <- data.frame(site = plots$plot, park = plots$Park)

design <- function(scat_verified) {
  sign <- function(name, y) {
    data.frame(site = pika$plot, sign = name, y = as.integer(y), year = year)
  }
  auto <- rbind(sign("scat", pika$FS > 0 & !scat_verified), sign("haypile", pika$FH > 0),
                sign("call", pika$PC > 0))
  verified <- data.frame(site = pika$plot, y = as.integer(pika$PS > 0 | scat_verified),
                         year = year)
  fp_data(auto, sites, verified, season = "year")
}
designs <- list(base = design(u), no_rule = design(u & pika$FH == 0),
                sighting = design(rep(FALSE, nrow(pika))))

fit <- function(d) {
  fp_occu_dyn(d, psi1 = ~ park, col = ~ 1, ext = ~ 1, priors = fp_priors(), det = ~ sign,
              fp = ~ sign, random_recorder = FALSE, chains = 4, parallel_chains = 4,
              threads_per_chain = 2, iter_warmup = 1000, iter_sampling = 1000,
              adapt_delta = 0.95, seed = 2026, refresh = 0, show_messages = FALSE)
}
fits <- lapply(designs, fit)

keep <- c("col_(Intercept)", "ext_(Intercept)", "fp_(Intercept)", "fp_gap", "p11_ref", "p10_ref",
          "p10_bar", "p_verified", "psi_eq")
table <- do.call(rbind, lapply(names(fits), function(n) {
  s <- as.data.frame(summary(fits[[n]])$table)
  cbind(design = n, s[s$variable %in% keep, c("variable", "mean", "lower", "upper")])
}))
ppc <- function(f, method) fp_ppc(f, method = method, ndraws = 300, seed = 1)[c("p_value", "fit_obs", "fit_rep")]
res <- list(verified_detections = sapply(designs, function(d) sum(d$verified$y)),
            diagnostics = lapply(fits, `[[`, "diagnostics"), table = table,
            ppc_auto = lapply(fits, ppc, method = "automated"),
            ppc_verified = lapply(fits, ppc, method = "verified"))
saveRDS(res, file.path(dir, "sensitivity.rds"))
str(res[c("verified_detections", "diagnostics")], max.level = 3)
print(sapply(res$ppc_auto, `[[`, "p_value"))
print(sapply(res$ppc_verified, `[[`, "p_value"))
print(table, digits = 3)
