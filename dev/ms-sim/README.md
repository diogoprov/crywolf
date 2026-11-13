# Simulations of the multispecies model (`fp_occu_ms()`)

`recovery.R` is the short recovery check of PR 2: one simulated community
(12 species, 80 sites, 8 automated occasions, 3 verified occasions at 30%
of the sites; occupancy `~ x` and detection `~ bio` by species, constant
false-positive rate with a gap by species), fitted once with 4 chains of
1000 + 1000 iterations. The simulator in the script is provisional;
`simulate_fp_ms()` comes in PR 4, and the simulations and SBC on AWS in
PR 5. Results in `recovery.rds`.

## Run of 02/11/2026

On the Mac, 2 threads per chain, 7.1 minutes. 43 of the 46 95% intervals
(community means, standard deviations among species, occupancy and
detection intercepts and `p10_ref` of each species) hold the true value;
the three misses are two occupancy intercepts and one `p10_ref`. Largest
R-hat 1.007, no divergences.

## PR 5: `run_sims.R` on AWS

`run_sims.R` holds the simulations of PR 5 (design in its header): recovery
with 5, 20 and 50 species, the "fs" and "mixed" scenarios of random effects
on logit(p10), and SBC. Timing pilot of 11/11/2026 on the Mac (ITER = 100,
one core per fit): 243 s for 5 species, 2,114 s for 20, 4,467 s for 50, 698 s
for the SBC fit. Run on spot with the recipe of `dev/dyn-sim/aws/` (commit
first):

```bash
export AWS_PROFILE=pcps
B=pcps-sim-560246495553
dev/dyn-sim/aws/pack_code.sh $B crywolf-ms
ID=$(JOBS=dev/ms-sim/jobs AZ=us-east-1d dev/dyn-sim/aws/launch_spot.sh $B crywolf-ms dev/ms-sim/run_sims.R c7a.48xlarge 190 pcps-sim-ec2 24)
JOBS=dev/ms-sim/jobs dev/dyn-sim/aws/supervise.sh $B crywolf-ms dev/ms-sim/run_sims.R 190 pcps-sim-ec2 24 $ID c7a.48xlarge c7i.48xlarge m7a.48xlarge
aws s3 sync s3://$B/fpcalib/runs/crywolf-ms/jobs dev/ms-sim/jobs
```
