# Comparison with a published analysis, run only with CRYWOLF_EXTENDED_TESTS=true
# (see tests/README.md). The data are downloaded from Zenodo, or read from the
# folder in CRYWOLF_CURIANGO_DIR when it holds the two files.

curiango_file <- function(f) {
  dir <- Sys.getenv("CRYWOLF_CURIANGO_DIR")
  if (nzchar(dir) && file.exists(file.path(dir, f))) return(file.path(dir, f))
  path <- file.path(tempdir(), f)
  if (file.exists(path)) return(path)
  url <- paste0("https://zenodo.org/api/records/21795898/files/", f, "/content")
  old <- options(timeout = max(3600, getOption("timeout")))
  on.exit(options(old))
  ok <- tryCatch(utils::download.file(url, path, mode = "wb", quiet = TRUE) == 0,
                 error = function(e) FALSE, warning = function(w) FALSE)
  if (!ok) {
    unlink(path)
    testthat::skip(paste("could not download", f, "from Zenodo"))
  }
  path
}

# The window-level data of the Custom classifier at its MaxF1 threshold (0.40),
# assembled as in the published analysis: one automated occasion per window,
# site-nights as sites, recorders as a random effect, and occupancy covariates
# averaged over each site-night.
curiango_windows <- function(thr = 0.40) {
  zc <- function(x) as.numeric(scale(x))
  df <- utils::read.table(curiango_file("df_night_clean.txt"), header = TRUE,
                          sep = "\t", stringsAsFactors = FALSE)
  df <- df[order(df$night_site_id, df$DataHora_Formatada), ]
  auto <- data.frame(site = df$night_site_id, y = as.integer(df$conf_div >= thr),
                     recorder = df$Site, bio = zc(df$BIO), aci = zc(df$ACI))
  key <- paste(df$Site, df$DataHora_Formatada)
  sites <- stats::aggregate(cbind(temp = Temp, rh = RH, moon = moon_illumination) ~
                              night_site_id, data = df, FUN = mean)
  sites <- data.frame(site = sites$night_site_id, temp = zc(sites$temp),
                      rh = zc(sites$rh), moon = zc(log(sites$moon)))
  cm <- utils::read.csv(curiango_file("CM-Custom.txt"), stringsAsFactors = FALSE)
  parts <- strsplit(cm$filename, "_")
  cm$key <- paste(cm$Site, vapply(parts, function(p) {
    paste(p[length(p) - 1], p[length(p)], sep = "_")
  }, character(1)))
  ver <- stats::aggregate(y ~ key, data.frame(key = cm$key,
                                              y = as.integer(cm$compau_annotation == 1)),
                          FUN = max)
  ver$site <- df$night_site_id[match(ver$key, key)]
  ver <- ver[!is.na(ver$site), c("site", "y")]
  # recording effort differs among site-nights by design, so the warning on
  # unequal effort is expected here
  suppressWarnings(fp_data(auto, sites, ver, recorder = "recorder"))
}

#' @srrstats {G5.4c} Comparison with the published estimates of Piatti et al. (2026), fitted with an independent Stan program on the same window-level data.
#' @srrstats {G5.11, G5.11a} The data are downloaded from Zenodo, and the test is skipped if the download fails.
test_that("estimates agree with a published analysis of the same data (extended)", {
  skip_if_not_extended()
  d <- curiango_windows()
  m <- suppressWarnings(fp_occu(
    d, occ = ~ temp + rh + moon, det = ~ bio + aci, fp = ~ 1,
    random_recorder = TRUE,
    priors = fp_priors(intercept_sd = 2.5, slope_sd = 2.5, recorder_sd = 2,
                       gap_mean = 0, gap_sd = 2.5),
    chains = 4, iter_warmup = 1000, iter_sampling = 1000, adapt_delta = 0.9,
    seed = 42
  ))
  # posterior means and sds of the published fit (Custom classifier, MaxF1)
  pub <- data.frame(
    par = c("occ_(Intercept)", "occ_temp", "occ_rh", "occ_moon"),
    mean = c(-1.55196847, -0.14004172, -0.57005637, 0.06017566),
    sd = c(0.06614236, 0.05892756, 0.05821186, 0.06818884)
  )
  post <- function(par) posterior::extract_variable(m$draws, par)
  est <- vapply(pub$par, function(p) mean(post(p)), numeric(1))
  # the two programs fit the same model, so the means differ only by Monte
  # Carlo error, far below a quarter of a posterior sd
  expect_true(all(abs(est - pub$mean) < 0.25 * pub$sd),
              label = paste(sprintf("%s %.3f vs %.3f", pub$par, est, pub$mean),
                            collapse = "; "))
  # the published p10: mean 3.081e-3, 90% interval 2.901e-3 to 3.267e-3
  p10 <- mean(post("p10_ref"))
  expect_lt(abs(p10 - 3.080786e-03), 0.25 * (3.266838e-03 - 2.900679e-03) / 3.29)
})
