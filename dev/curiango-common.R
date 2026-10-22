# Curiango data in the form used by the precomputed vignette results:
# custom classifier, one-minute occasions, nights as sites, covariates
# standardized. Sourced by dev/precompute-markov.R and
# dev/precompute-threshold.R.
zs <- function(v) as.numeric(scale(v))
m <- curiango$minutes
nt <- curiango$nights
cur_occasions <- data.frame(file = paste(m$site_night, format(m$time, "%Y%m%d_%H%M")),
                            site_night = m$site_night, time = m$time,
                            bio = zs(m$bio), aci = zs(m$aci), recorder = m$site)
cur_detections <- data.frame(file = cur_occasions$file, confidence = m$conf_custom)
cur_sites <- data.frame(site_night = nt$site_night, temp = zs(nt$temp), rh = zs(nt$rh),
                        moon = zs(log(nt$moon)))
cur_verified <- data.frame(site_night = curiango$verified$site_night,
                           y = curiango$verified$y)
