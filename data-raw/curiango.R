# Builds the `curiango` dataset from the public deposit of Piatti et al. (2026),
# Zenodo record 21795898 (CC BY 4.0): https://doi.org/10.5281/zenodo.21795898
#
# The source table (df_night_clean.txt) has one to 20 rows (detection windows)
# per one-minute recording. Acoustic indices and weather are constant within a
# minute; a classifier's confidence score is non-zero only in the windows where
# it detected the species (the custom score usually repeats across all windows
# of a minute, the BirdNET score usually does not). We keep one row per minute,
# the occasion of the occupancy model, with the maximum score over its windows. The verified data are
# the independent ground-truth test set (CM-Custom.txt): about 1% of the
# recordings, chosen in proportion to the effort at each site and annotated by
# ear in 3-s segments; a minute is a verified detection if any of its segments
# holds the species.
#
# Run from the package root: Rscript data-raw/curiango.R

library(dplyr)
options(timeout = max(600, getOption("timeout")))

src <- file.path(tempdir(), "curiango-zenodo")
dir.create(src, showWarnings = FALSE)
get <- function(f) {
  dst <- file.path(src, f)
  if (!file.exists(dst)) {
    utils::download.file(
      paste0("https://zenodo.org/records/21795898/files/", f, "?download=1"),
      dst, mode = "wb", quiet = TRUE)
  }
  dst
}

raw <- utils::read.delim(get("df_night_clean.txt"))

# within a minute, a classifier's score is non-zero only in the windows where
# it detected the species, so the minute takes the maximum over its windows;
# the other columns are constant within a minute
minutes <- raw |>
  group_by(Site, DataHora_Formatada) |>
  summarise(
    site_night = first(night_site_id),
    time = as.POSIXct(first(datetime), tz = "America/Campo_Grande"),
    conf_custom = max(conf_div),
    conf_birdnet = max(conf_preb),
    bio = first(BIO),
    aci = first(ACI),
    .groups = "drop"
  ) |>
  transmute(site_night, site = Site, time, conf_custom, conf_birdnet, bio, aci,
            key = DataHora_Formatada) |>
  arrange(site_night, time)

nights <- raw |>
  distinct(Site, DataHora_Formatada, .keep_all = TRUE) |>
  group_by(site_night = night_site_id) |>
  summarise(site = first(Site), night = as.Date(first(night)),
            temp = mean(Temp), rh = mean(RH),
            moon = mean(moon_illumination), .groups = "drop") |>
  arrange(site_night)

cm <- utils::read.csv(get("CM-Custom.txt"))
parts <- strsplit(cm$filename, "_")
cm$key <- vapply(parts, function(p) paste(p[length(p) - 1], p[length(p)], sep = "_"),
                 character(1))
verified <- cm |>
  group_by(site = Site, key) |>
  summarise(y = as.integer(any(compau_annotation == 1)), .groups = "drop") |>
  inner_join(distinct(minutes, site, key, site_night, time), by = c("site", "key")) |>
  transmute(site_night, site, time, y) |>
  arrange(site_night, time)

minutes$key <- NULL
curiango <- list(minutes = as.data.frame(minutes),
                 nights = as.data.frame(nights),
                 verified = as.data.frame(verified))
message(sprintf("%d minutes, %d site-nights, %d verified minutes (%d dropped: no match)",
                nrow(minutes), nrow(nights), nrow(verified),
                n_distinct(paste(cm$Site, cm$key)) - nrow(verified)))
usethis::use_data(curiango, overwrite = TRUE, compress = "xz")
