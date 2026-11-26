#!/bin/sh
# Runs dev/jsdm-export/grid.R in two jobs (4 chains each), keeping the Mac awake.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NJOBS=2
JOB=1 Rscript dev/jsdm-export/grid.R > dev/jsdm-export/grid-job1.log 2>&1 &
p1=$!
JOB=2 Rscript dev/jsdm-export/grid.R > dev/jsdm-export/grid-job2.log 2>&1 &
p2=$!
caffeinate -w "$p1" &
caffeinate -w "$p2" &
wait
echo GRID_DONE
