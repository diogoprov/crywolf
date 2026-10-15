# Full-length curiango fit (custom classifier, fp ~ 1, recorder effect, as in
# dev/precompute-curiango.R) with threads, compared with the saved fit of the
# previous model in dev/curiango-fits/custom_constant.rds.
# Usage, from the package root: Rscript dev/profile/full-fit-curiango.R [threads]
devtools::load_all(quiet = TRUE)
source("dev/profile/curiango-data.R")
args <- commandArgs(trailingOnly = TRUE)
threads <- if (length(args)) as.integer(args[1]) else 4
zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
auto <- data.frame(site_night = m$site_night,
                   y = as.integer(m$conf_custom >= 0.40),
                   bio = zs(m$bio), aci = zs(m$aci), recorder = m$site)
n <- curiango$nights
sites <- data.frame(site_night = n$site_night, temp = zs(n$temp),
                    rh = zs(n$rh), moon = zs(log(n$moon)))
ver <- data.frame(site_night = curiango$verified$site_night,
                  y = curiango$verified$y)
d <- fp_data(auto, sites, ver, site = "site_night", recorder = "recorder")
fp_compile()
t0 <- Sys.time()
fit <- fp_occu(d, occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1,
               random_recorder = TRUE, seed = 1, threads_per_chain = threads)
minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat("threads_per_chain =", threads, ":", round(minutes, 1), "min\n")
old <- readRDS("dev/curiango-fits/custom_constant.rds")
new_tab <- summary(fit)$table
cmp <- merge(old$table[, c("variable", "mean", "sd")],
             new_tab[, c("variable", "mean", "sd")], by = "variable",
             suffixes = c("_old", "_new"))
cmp$z <- (cmp$mean_new - cmp$mean_old) / cmp$sd_old
print(cmp, digits = 3)
print(rbind(old = unlist(old$diagnostics), new = unlist(fit$diagnostics)))
lo_old <- old$loo$estimates; lo_new <- loo(fit)$estimates
print(rbind(old = lo_old["elpd_loo", ], new = lo_new["elpd_loo", ]))
saveRDS(list(minutes = minutes, old_minutes = old$minutes, threads = threads,
             compare = cmp, diagnostics = list(old = old$diagnostics,
                                               new = fit$diagnostics),
             loo = list(old = lo_old, new = lo_new)),
        "dev/profile/full-fit-curiango.rds")
