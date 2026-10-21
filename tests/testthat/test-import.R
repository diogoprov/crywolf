write_table <- function(path, rows) {
  header <- paste("Selection", "View", "Channel", "Begin Time (s)",
                  "End Time (s)", "Begin Path", "File Offset (s)",
                  "Species Code", "Common Name", "Confidence", sep = "\t")
  writeLines(c(header, rows), path)
}
row <- function(i, file, offset, code, conf) {
  paste(i, "Spectrogram 1", 1, offset, offset + 3, file, offset, code,
        "x", conf, sep = "\t")
}

test_that("fp_read_birdnet reads and standardizes selection tables", {
  dir <- tempfile()
  dir.create(file.path(dir, "output_SITE1"), recursive = TRUE)
  write_table(file.path(dir, "output_SITE1", "a.txt"), c(
    row(1, "N:\\\\rec\\\\SITE1_20240101_200000.wav", 0, "compau", 0.81),
    row(2, "/rec/SITE1_20240101_200000.wav", 63, "compau", 0.2)))
  write_table(file.path(dir, "b.txt"), row(1, "/rec/noStamp.wav", 3, "other", 0.5))
  d <- fp_read_birdnet(dir)
  expect_equal(nrow(d), 3)
  expect_true(all(c("begin_time", "confidence", "species_code", "file",
                    "start", "time", "table") %in% names(d)))
  s1 <- d[d$species_code == "compau", ]
  expect_equal(s1$file, rep("SITE1_20240101_200000.wav", 2))
  expect_equal(format(s1$time, "%H:%M:%S"), c("20:00:00", "20:01:03"))
  expect_true(is.na(d$time[d$file == "noStamp.wav"]))
  bad <- tempfile(fileext = ".txt")
  writeLines("Selection\tConfidence", bad)
  expect_error(fp_read_birdnet(bad), "lacks column")
})

test_that("fp_aggregate builds detections per occasion", {
  occ <- data.frame(site = c("A", "A", "B", "B"),
                    file = c("a1", "a2", "b1", "b2"), cov = 1:4)
  det <- data.frame(file = c("a1", "a1", "b1", "zz"),
                    species_code = c("sp", "sp", "sp", "sp"),
                    confidence = c(0.9, 0.6, 0.3, 0.9))
  expect_warning(out <- fp_aggregate(occ, det, threshold = 0.5), "not in")
  expect_equal(out$y, c(1L, 0L, 0L, 0L))
  expect_equal(out$n_detections, c(2L, 0L, 0L, 0L))
  expect_equal(out$score, c(0.9, 0, 0.3, 0))
  expect_equal(out$cov, 1:4)
  out2 <- suppressWarnings(fp_aggregate(occ, det, 0.2, species = "other"))
  expect_equal(sum(out2$y), 0)
  expect_error(fp_aggregate(rbind(occ, occ), det, 0.5), "one row per occasion")
})

test_that("verdicts become confirmation counts, with checks", {
  occ <- data.frame(site = c("A", "A", "B"), file = c("a1", "a2", "b1"))
  det <- data.frame(file = c("a1", "a2", "a2", "b1"), species_code = "sp",
                    confidence = c(0.9, 0.8, 0.7, 0.9))
  val <- data.frame(file = c("a1", "a2", "b1", "b1"), species_code = "sp",
                    confidence = c(0.9, 0.8, 0.9, 0.3),
                    verdict = c(TRUE, FALSE, NA, TRUE))
  expect_warning(expect_message(
    out <- fp_aggregate(occ, det, threshold = 0.5, validation = val),
    "below the threshold"), "more detections than")
  expect_equal(out$n_true, c(1L, 0L, NA))
  expect_equal(out$n_false, c(0L, 1L, NA))
  expect_warning(
    d <- fp_data(out, data.frame(site = c("A", "B")), confirmed = c("n_true", "n_false")),
    "Every automated occasion is a detection")
  expect_s3_class(d, "fp_data")
})

test_that("fp_read_validation maps the sampler's verdicts", {
  f <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(
    local = "S1", species_code = "sp", audio_origem = "D:\\rec\\S1_x.wav",
    inicio_clipe = c(1, 5), fim_clipe = c(4, 8), confidence = c(0.9, 0.7),
    veredito = c("positivo", "nao_validado")), f, row.names = FALSE)
  v <- fp_read_validation(f)
  expect_equal(v$verdict, c(TRUE, NA))
  expect_equal(v$file, c("S1_x.wav", "S1_x.wav"))
  utils::write.csv(data.frame(a = 1), f, row.names = FALSE)
  expect_error(fp_read_validation(f), "validacao_por_clipe")
})

test_that("fp_data warns about sites without occasions and constant data", {
  auto <- data.frame(site = c(1, 1, 2), y = c(0, 1, 0))
  expect_warning(fp_data(auto, data.frame(site = 1:3)), "no automated occasion")
  expect_warning(fp_data(data.frame(site = 1:2, y = 0), data.frame(site = 1:2)),
                 "non-detection")
})

test_that("empty tables, several species, NA scores and scaling are handled", {
  dir <- tempfile(); dir.create(dir)
  write_table(file.path(dir, "empty.txt"), character(0))
  write_table(file.path(dir, "full.txt"), c(
    row(1, "/r/S_20240101_200000.wav", 0, "sp", 0.9),
    row(2, "/r/S_20240101_200000.wav", 3, "other", 0.8)))
  d <- fp_read_birdnet(dir)
  expect_equal(nrow(d), 2)
  expect_type(d$confidence, "double")
  occ <- data.frame(site = "S", file = "S_20240101_200000.wav", bio = 5)
  expect_warning(fp_aggregate(occ, d, 0.5), "several species")
  d$confidence[1] <- NA
  expect_warning(out <- fp_aggregate(occ, d, 0.5, species = "sp"), "without a confidence")
  expect_equal(out$y, 0L)
  occ2 <- data.frame(site = c("S", "S"), file = c("a", "b"), bio = c(1, 3))
  out2 <- fp_aggregate(occ2, data.frame(file = "a", confidence = 0.9), 0.5,
                       standardize = "bio")
  expect_equal(out2$bio, c(-1, 1) / sqrt(2))
  expect_equal(attr(out2, "scaling")$center, 2)
})

test_that("fp_data warns about unequal effort and data without positives", {
  auto <- data.frame(site = c(1, rep(2, 12)), y = c(1, rep(0, 11), 1))
  expect_warning(fp_data(auto, data.frame(site = 1:2)), "very unequal")
  auto <- data.frame(site = c(1, 1, 2, 2), y = c(0, 1, 0, 1))
  expect_warning(fp_data(auto, data.frame(site = 1:2),
                         verified = data.frame(site = 1, y = 0)),
                 "No verified occasion is a detection")
  auto$t <- c(NA, 1, NA, 1); auto$f <- c(NA, 0, NA, 0)
  expect_warning(fp_data(auto, data.frame(site = 1:2), confirmed = c("t", "f")),
                 "No confirmed occasion holds a false observation")
})
