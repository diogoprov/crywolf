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

# runs of consecutive detections within a night (15-minute steps)
runs <- function(a) {
  a <- a[order(a$site_night, a$time), ]
  same <- c(FALSE, a$site_night[-1] == a$site_night[-nrow(a)])
  step <- c(NA, diff(as.numeric(a$time)))
  adjacent <- same & step == 900
  start <- a$y == 1 & !(adjacent & c(0, a$y[-nrow(a)]) == 1)
  id <- cumsum(start)
  len <- tabulate(id[a$y == 1])
  table(factor(pmin(len, 5), levels = 1:5, labels = c(1:4, "5+")))
}
