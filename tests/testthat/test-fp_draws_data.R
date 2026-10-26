# a small export with calibration draws, as sheepdog writes it: two sites,
# grid of scores from 0.1 in steps of 0.1, two draws
write_cal_export <- function(dir) {
  dir.create(dir, showWarnings = FALSE)
  w <- function(d, name) utils::write.csv(d, file.path(dir, name), row.names = FALSE, na = "")
  w(data.frame(site = c("A", "A", "B", "B"), occasion_id = c("a1", "a2", "b1", "b2"),
               night = c("n1", "n2", "n1", "n1")), "occasions.csv")
  w(data.frame(site = c("A", "A", "B", "B"), occasion_id = c("a1", "a1", "b1", "b2"),
               species_code = "Sp one_sp-one", species = "sp_one",
               file_offset = c(0, 3, 0, 0), confidence = c(0.92, 0.18, 0.5, 0.31)),
    "detections.csv")
  w(data.frame(site = "B", occasion_id = "b2", species_code = "Sp one_sp-one",
               species = "sp_one", inclusion_prob = 0.1, verdict = "positive", present = 1),
    "recordings_verified.csv")
  w(data.frame(site = c("A", "B"), x = c(-1, 1)), "sites.csv")
  grid <- expand.grid(grid_index = 1:9, site = c("A", "B"), draw = c(3L, 8L),
                      stringsAsFactors = FALSE)
  # draw 3: every prediction true; draw 8: true only from a score of 0.9
  grid$p_tp <- ifelse(grid$draw == 3L, 1, as.numeric(grid$grid_index >= 9))
  grid$species <- "sp_one"
  grid$confidence <- 0.1 * grid$grid_index
  arrow::write_parquet(grid, file.path(dir, "calibration_draws.parquet"))
  jsonlite::write_json(list(schema = "sheepdog-occupancy", schema_version = 1L,
                            sheepdog_version = "3.2.1.9000", occasion_unit = "file",
                            min_confidence = 0.1,
                            calibration = list(engine = "brms", n_draws = 2L,
                                               grid_step = 0.1, grid_min = 0.1)),
                       file.path(dir, "schema.json"), auto_unbox = TRUE)
  dir
}

test_that("each calibration draw gives one version of the detections", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export")))
  v <- fp_draws_data(x, rule = "threshold", p_min = 0.5)
  expect_length(v, 2)
  expect_equal(attr(v, "draws"), c(3L, 8L))
  expect_s3_class(v[[1]], "fp_data")
  # draw 3: every occasion with a prediction; draw 8: only a1 (score 0.92)
  expect_equal(v[[1]]$auto$y, c(1L, 0L, 1L, 1L))
  expect_equal(v[[2]]$auto$y, c(1L, 0L, 0L, 0L))
  expect_equal(v[[1]]$verified$y, 1L)
  # probabilities of 0 and 1 make the sampled versions equal to these
  s <- fp_draws_data(x, seed = 1)
  expect_equal(s[[1]]$auto$y, v[[1]]$auto$y)
  expect_equal(s[[2]]$auto$y, v[[2]]$auto$y)
  expect_length(fp_draws_data(x, n_draws = 1), 1)
})

test_that("a calibration by month is matched by the month of recording", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-month"))
  occ <- utils::read.csv(file.path(dir, "occasions.csv"))
  occ$recording_date <- c("20240101", "20240201", "20240101", "20240201")
  utils::write.csv(occ, file.path(dir, "occasions.csv"), row.names = FALSE)
  cal <- as.data.frame(arrow::read_parquet(file.path(dir, "calibration_draws.parquet")))
  # January: every prediction true; February: none
  cal <- rbind(transform(cal, month = "202401", p_tp = 1),
               transform(cal, month = "202402", p_tp = 0))
  arrow::write_parquet(cal, file.path(dir, "calibration_draws.parquet"))
  sc <- jsonlite::read_json(file.path(dir, "schema.json"))
  sc$calibration$keys <- list("site", "species", "month", "grid_index")
  jsonlite::write_json(sc, file.path(dir, "schema.json"), auto_unbox = TRUE)
  x <- fp_read_sheepdog(dir)
  v <- fp_draws_data(x, rule = "threshold", p_min = 0.5)
  # a1 and b1 recorded in January, b2 in February (a2 has no prediction)
  expect_equal(v[[1]]$auto$y, c(1L, 0L, 1L, 0L))
})

test_that("predictions without calibration leave their occasions unsurveyed", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-nocal"))
  cal <- as.data.frame(arrow::read_parquet(file.path(dir, "calibration_draws.parquet")))
  arrow::write_parquet(cal[cal$site == "A", ], file.path(dir, "calibration_draws.parquet"))
  x <- fp_read_sheepdog(dir)
  expect_warning(
    v <- withCallingHandlers(fp_draws_data(x, rule = "threshold", p_min = 0.5),
                             crywolf_design_warning = function(w) invokeRestart("muffleWarning")),
    "2 occasion\\(s\\) are treated as not surveyed")
  # fp_data() keeps surveyed occasions only: site B drops out of auto
  expect_equal(v[[1]]$auto$y, c(1L, 0L))
  expect_equal(unique(v[[1]]$auto$site), "A")
  arrow::write_parquet(cal[cal$site == "A" & cal$species == "none", ], file.path(dir, "calibration_draws.parquet"))
  expect_error(suppressWarnings(fp_draws_data(fp_read_sheepdog(dir))),
               "No prediction could be matched")
})

test_that("keys of schema.json and the parquet must agree", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-keys"))
  cal <- as.data.frame(arrow::read_parquet(file.path(dir, "calibration_draws.parquet")))
  # a month column the schema does not announce
  arrow::write_parquet(rbind(transform(cal, month = "202401"),
                             transform(cal, month = "202402")),
                       file.path(dir, "calibration_draws.parquet"))
  expect_error(fp_draws_data(fp_read_sheepdog(dir)), "disagree")
  # a key column missing from both: two values per key and draw
  arrow::write_parquet(rbind(cal, cal), file.path(dir, "calibration_draws.parquet"))
  expect_error(fp_draws_data(fp_read_sheepdog(dir)), "more than one value")
})

test_that("the random numbers of the session are left untouched", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-rng")))
  set.seed(42)
  ref <- stats::runif(1)
  set.seed(42)
  invisible(fp_draws_data(x, seed = 1))
  expect_identical(stats::runif(1), ref)
})

test_that("duplicated occasions are refused", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-dup")))
  occ <- rbind(x$occasions, x$occasions[1, ])
  expect_error(fp_draws_data(x, occasions = occ), "one row per site")
})

test_that("the sites of the model can be other than the recorders", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  x <- fp_read_sheepdog(write_cal_export(file.path(tempdir(), "cal-export-nights")))
  occ <- x$occasions
  occ$site_night <- paste(occ$site, occ$night)
  nights <- data.frame(site_night = unique(occ$site_night))
  v <- fp_draws_data(x, rule = "threshold", occasions = occ, sites = nights,
                     site = "site_night")
  expect_equal(nrow(v[[1]]$sites), 3)
  expect_equal(v[[1]]$verified$site_night, "B n1")
})

test_that("an export without calibration draws is refused", {
  skip_if_not_installed("jsonlite")
  dir <- write_cal_export(file.path(tempdir(), "cal-export-none"))
  unlink(file.path(dir, "calibration_draws.parquet"))
  x <- fp_read_sheepdog(dir)
  expect_error(fp_draws_data(x), "no calibration draws")
})
