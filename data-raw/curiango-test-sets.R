# Writes data-raw/curiango-test-sets.csv: the test set (1 to 3) of each
# recording annotated by ear in the `curiango` dataset.
#
# The three test sets are folders kept by L. Piatti (Cjto-Teste_01 to _03);
# they are not part of the public deposit. Each holds a BirdNET run with a
# list of the files it analysed (BirdNET_SelectionTable.list.txt). Sets 1 and
# 2 are random draws; set 3 completes them with recordings chosen by BirdNET
# (see ?curiango). Two recordings are in set 3 and also in set 1 or 2; they
# keep the lower set, the random draw.
#
# Run from the package root: Rscript data-raw/curiango-test-sets.R <folder>
# where <folder> holds 03_output_Birdnet_ct01, 03_output_Birdnet-ct02 and
# 03_output_Birdnet-ct03.

base <- commandArgs(TRUE)[1]
dirs <- c("03_output_Birdnet_ct01", "03_output_Birdnet-ct02", "03_output_Birdnet-ct03")
sets <- do.call(rbind, lapply(seq_along(dirs), function(i) {
  f <- readLines(file.path(base, dirs[i], "BirdNET_SelectionTable.list.txt"))
  f <- sub("\\.flac$", "", sub(".*\\\\", "", f))  # Windows paths
  # DIVSPTS01_44p1K_<site>_<recorder>_<yyyymmdd>_<hhmmss>, in local time
  data.frame(site = sub("^DIVSPTS01_44p1K_([A-Z]+)_.*", "\\1", f),
             time = format(as.POSIXct(sub(".*_(\\d{8}_\\d{6})$", "\\1", f),
                                      format = "%Y%m%d_%H%M%S", tz = "UTC"),
                           "%Y-%m-%d %H:%M:%S"),
             test_set = i)
}))
sets <- sets[order(sets$test_set), ]
sets <- sets[!duplicated(sets[c("site", "time")]), ]
sets <- sets[order(sets$site, sets$time), ]
utils::write.csv(sets, "data-raw/curiango-test-sets.csv", row.names = FALSE)
message(nrow(sets), " recordings: ", paste(table(sets$test_set), collapse = ", "))
