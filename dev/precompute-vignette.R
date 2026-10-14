# Pre-computes the recorder model shown in vignettes/crywolf.Rmd
# ("Detection covariates and recorder effects"). It takes about a minute, so
# the vignette reads the saved summary instead of refitting it on every build.
#
# Run from the package root after any change to the Stan model or to
# fp_occu(), and commit the .rds:
#   Rscript dev/precompute-vignette.R

devtools::load_all(quiet = TRUE)

sim_rec <- simulate_fp(nsite = 300, nocc = 20, a_bio = 0.5, a_aci = -0.3,
                       sd_recorder = 0.6, nrec = 10,
                       conf_effect = 1.5, conf_cor = 0.7, seed = 2027)
m_rec <- fp_occu(sim_rec$data, occ = ~ x, det = ~ bio + aci, fp = ~ w,
                 random_recorder = TRUE, seed = 1)

rec_fit <- list(
  printed = utils::capture.output(print(m_rec)),
  table = summary(m_rec)$table,
  diagnostics = m_rec$diagnostics,
  truth = sim_rec$truth[c("b1", "a_bio", "a_aci", "sd_recorder",
                          "conf_effect")],
  info = list(date = Sys.Date(),
              crywolf = as.character(utils::packageVersion("crywolf")),
              cmdstan = as.character(cmdstanr::cmdstan_version()))
)
saveRDS(rec_fit, "vignettes/crywolf-recorder-fit.rds")
