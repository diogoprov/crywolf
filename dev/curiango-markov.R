# Fits the curiango data (custom classifier, MaxF1 threshold 0.40, fp ~ 1,
# det ~ bio + aci, recorder effect) with and without the Markov term
# (markov = "time"; recordings every 15 minutes, so a 30-minute interval is
# a gap that starts a new run: 1,747 of the 96,053 occasions that are not
# the first of their night).
# Result of 21/10/2026 (4 chains x 4 threads, 6 to 7 min each): det_y_prev
# 2.14 (2.01 to 2.28); LOO favours the Markov model by 450 (SE 55); PPC
# discrepancy observed/replicated 568/400 without and 429/350 with it, p = 0
# in both; occ_rh -0.78 -> -0.81, occ_temp -0.25 in both; p10_ref 0.33% ->
# 0.25%. Run from the package root: Rscript dev/curiango-markov.R
devtools::load_all(quiet = TRUE)
zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
auto <- data.frame(site_night = m$site_night, time = m$time,
                   y = as.integer(m$conf_custom >= 0.40),
                   bio = zs(m$bio), aci = zs(m$aci), recorder = m$site)
n <- curiango$nights
sites <- data.frame(site_night = n$site_night, temp = zs(n$temp), rh = zs(n$rh), moon = zs(log(n$moon)))
ver <- data.frame(site_night = curiango$verified$site_night, y = curiango$verified$y)
d <- suppressWarnings(fp_data(auto, sites, ver, site = "site_night", recorder = "recorder"))
fit <- function(markov) {
  t0 <- Sys.time()
  f <- fp_occu(d, occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1, markov = markov,
               random_recorder = TRUE, seed = 1, chains = 4, parallel_chains = 4,
               threads_per_chain = 4, refresh = 0, show_messages = FALSE)
  list(f = f, min = as.numeric(difftime(Sys.time(), t0, units = "mins")))
}
r0 <- fit(NULL); r1 <- fit("time")
p0 <- fp_ppc(r0$f, ndraws = 300, seed = 1); p1 <- fp_ppc(r1$f, ndraws = 300, seed = 1)
cat("minutes: no markov", round(r0$min, 1), " markov", round(r1$min, 1), "\n")
cat("PPC p-value: no markov", p0$p_value, " markov", p1$p_value, "\n")
cat("mean discrepancy obs/rep no markov:", mean(p0$fit_obs), mean(p0$fit_rep), " markov:", mean(p1$fit_obs), mean(p1$fit_rep), "\n")
tab <- function(f) { t <- summary(f)$table; t[t$variable %in% c("occ_temp","occ_rh","occ_moon","det_y_prev","p10_ref"), c("variable","mean","lower","upper")] }
print(tab(r0$f)); print(tab(r1$f))
print(r0$f$diagnostics); print(r1$f$diagnostics)
print(loo::loo_compare(list(no_markov = loo(r0$f), markov = loo(r1$f))))
saveRDS(list(p0 = p0[c("p_value","fit_obs","fit_rep")], p1 = p1[c("p_value","fit_obs","fit_rep")], t0 = tab(r0$f), t1 = tab(r1$f)), file.path("dev", "curiango-markov.rds"))
