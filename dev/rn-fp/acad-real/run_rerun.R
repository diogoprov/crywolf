# Refit, at adapt_delta 0.99, the fits of run_grid.R that had divergences
# (rerun-list.rds, built from the job files of the first run), with the same
# data and seeds. Job files go to dev/rn-fp/acad-real/jobs-ad99.
Sys.setenv(RERUN = "dev/rn-fp/acad-real/rerun-list.rds", ADAPT_DELTA = "0.99",
           JOBS = "dev/rn-fp/acad-real/jobs-ad99")
source("dev/rn-fp/acad-real/run_grid.R")
