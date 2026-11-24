# The grid of run_grid.R with the design of the FrogNet ACAD data (2023): 21
# sites, 17 occasions (mornings; 11-19 in the data, median 17), occupancy 0.8
# or 0.3 at the mean of x, and either every occasion verified or three. 24
# cells x 20 replicates, saved in dev/rn-fp/jobs-acad/. In the data the
# verified mornings are the automated ones; here, as in the model, they are
# separate occasions. Run from the package root:
#   WORKERS=16 Rscript dev/rn-fp/run_grid_acad.R
Sys.setenv(DESIGN = "acad")
source("dev/rn-fp/run_grid.R")
