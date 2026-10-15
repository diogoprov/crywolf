# Stan data of the curiango benchmark fit: custom classifier, MaxF1 threshold
# 0.40, occ ~ temp + rh + moon, det ~ bio + aci, fp ~ 1, recorder random
# intercept (the "constant" model of dev/precompute-curiango.R).
curiango_bench_data <- function(classifier = "custom", threshold = 0.40) {
  zs <- function(v) as.numeric(scale(v))
  m <- curiango$minutes
  auto <- data.frame(site_night = m$site_night,
                     y = as.integer(m[[paste0("conf_", classifier)]] >= threshold),
                     bio = zs(m$bio), aci = zs(m$aci), recorder = m$site)
  n <- curiango$nights
  sites <- data.frame(site_night = n$site_night, temp = zs(n$temp),
                      rh = zs(n$rh), moon = zs(log(n$moon)))
  ver <- data.frame(site_night = curiango$verified$site_night,
                    y = curiango$verified$y)
  d <- fp_data(auto, sites, ver, site = "site_night", recorder = "recorder")
  crywolf:::.build_stan_data(d, occ = ~ temp + rh + moon, det = ~ bio + aci,
                             fp = ~ 1, random_recorder = TRUE,
                             priors = fp_priors())$stan_data
}
