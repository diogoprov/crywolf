# Pre-computes the fits shown in vignettes/curiango.Rmd. Each fit takes
# about 40 minutes (97,806 one-minute occasions with per-minute detection
# covariates), so the vignette reads saved summaries instead of refitting.
#
# Usage, from the package root (one process per classifier, run in parallel):
#   Rscript dev/precompute-curiango.R custom
#   Rscript dev/precompute-curiango.R birdnet
#   Rscript dev/precompute-curiango.R combine   # writes vignettes/curiango-fits.rds
#
# Rerun after any change to the Stan model or to fp_occu().

devtools::load_all(quiet = TRUE)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- file.path("dev", "curiango-fits")
dir.create(out_dir, showWarnings = FALSE)

# MaxF1 thresholds of Piatti et al. (2026)
thresholds <- c(custom = 0.40, birdnet = 0.10)

curiango_data <- function(classifier) {
  zs <- function(v) as.numeric(scale(v))
  m <- curiango$minutes
  auto <- data.frame(site_night = m$site_night,
                     y = as.integer(m[[paste0("conf_", classifier)]] >=
                                      thresholds[[classifier]]),
                     bio = zs(m$bio), aci = zs(m$aci), recorder = m$site)
  n <- curiango$nights
  sites <- data.frame(site_night = n$site_night, temp = zs(n$temp),
                      rh = zs(n$rh), moon = zs(log(n$moon)))
  ver <- data.frame(site_night = curiango$verified$site_night,
                    y = curiango$verified$y)
  fp_data(auto, sites, ver, site = "site_night", recorder = "recorder")
}

models <- list(standard = NULL, constant = ~ 1, covariates = ~ temp + rh + moon)

fit_one <- function(classifier, model) {
  d <- curiango_data(classifier)
  t0 <- Sys.time()
  m <- fp_occu(d, occ = ~ temp + rh + moon, det = ~ bio + aci,
               fp = models[[model]], random_recorder = TRUE, seed = 1)
  minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  grid <- seq(-2.5, 2.5, length.out = 41)
  nd <- function(var) {
    d0 <- data.frame(temp = 0, rh = 0, moon = 0)[rep(1, length(grid)), ]
    d0[[var]] <- grid
    d0
  }
  pred <- function(process) {
    purrr::map(c("temp", "rh", "moon"), function(v) {
      data.frame(covariate = v, value = grid,
                 predict(m, newdata = nd(v), process = process)[,
                   c("estimate", "lower", "upper")])
    }) |> purrr::list_rbind()
  }
  ppc <- fp_ppc(m, ndraws = 200, seed = 1)
  list(
    classifier = classifier, model = model, minutes = minutes,
    printed = utils::capture.output(print(m)),
    table = summary(m)$table, diagnostics = m$diagnostics,
    loo = loo(m),
    ppc = ppc[c("y", "yrep", "fit_obs", "fit_rep", "p_value", "fit_stat")],
    pred_occ = pred("occ"),
    pred_fp = if (model == "covariates") pred("fp") else NULL
  )
}

if (identical(args[1], "combine")) {
  files <- list.files(out_dir, pattern = "\\.rds$", full.names = TRUE)
  fits <- purrr::map(files, readRDS)
  names(fits) <- purrr::map_chr(fits, function(f) paste(f$classifier, f$model, sep = "_"))
  fits <- purrr::map(fits, function(f) {
    f$ppc$yrep <- f$ppc$yrep[seq_len(100), , drop = FALSE]   # keep the file small
    f
  })
  attr(fits, "info") <- list(
    date = Sys.Date(), crywolf = as.character(utils::packageVersion("crywolf")),
    cmdstan = as.character(cmdstanr::cmdstan_version()))
  saveRDS(crywolf:::.stamp(fits), file.path("vignettes", "curiango-fits.rds"), compress = "xz")
} else {
  classifier <- match.arg(args[1], names(thresholds))
  # an optional second argument points to a pre-compiled copy of the model,
  # so that edits to inst/stan/ during a long run do not affect it
  if (length(args) > 1) {
    fp_env <- get(".fp_env", envir = asNamespace("crywolf"))
    fp_env$model <- cmdstanr::cmdstan_model(args[2])
  }
  purrr::walk(names(models), function(model) {
    f <- fit_one(classifier, model)
    saveRDS(f, file.path(out_dir, paste0(classifier, "_", model, ".rds")))
    message(classifier, " / ", model, ": ", round(f$minutes, 1), " min")
  })
}
