# Hmsc setup for the export of occupancy draws (README.md). In the pilot,
# Hmsc fitted to the true states recovered the true residual pattern at only
# 0.55, so the draws could not be judged. Here, on the true states of the
# pilot dataset, three choices are compared:
#   - the scale of the comparison: residual correlation (as in the pilot)
#     against the true correlation of LL' + I, or the residual covariance
#     Omega = LL' against the true LL';
#   - the number of latent factors: Hmsc's default (shrinkage, up to many)
#     or fixed at two, as in the simulation;
#   - the chain length.
# Then the chosen setup is applied to the 20 state draws of the pilot.
#   Rscript dev/jsdm-export/hmsc_setup.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/jsdm-export/simulate_jsdm_fp.R")
sim <- simulate_jsdm_fp(K = 2, seed = 1001)   # the pilot dataset
d <- sim$data
S <- d$S; J <- d$J
X <- data.frame(x = d$x)
LL <- sim$truth$LL
off <- lower.tri(LL)
true_cor <- stats::cov2cor(LL + diag(J))   # latent residual correlation, probit-like

fit_hmsc <- function(y, nf = NULL, samples = 1000, transient = 2500, thin = 5) {
  colnames(y) <- paste0("sp", seq_len(J))
  units <- factor(seq_len(S))
  rl <- Hmsc::HmscRandomLevel(units = units)
  if (!is.null(nf)) rl <- Hmsc::setPriors(rl, nfMin = nf, nfMax = nf)
  m <- Hmsc::Hmsc(Y = y, XData = X, XFormula = ~ x, distr = "probit",
                  studyDesign = data.frame(site = units), ranLevels = list(site = rl))
  t0 <- Sys.time()
  m <- Hmsc::sampleMcmc(m, samples = samples, thin = thin, transient = transient,
                        nChains = 2, nParallel = 2, verbose = 0)
  mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  post <- Hmsc::convertToCodaObject(m)
  psrf <- coda::gelman.diag(post$Omega[[1]], multivariate = FALSE)$psrf[, 1]
  assoc <- Hmsc::computeAssociations(m)[[1]]
  omega <- Hmsc::getPostEstimate(m, "Omega")$mean
  list(cor = assoc$mean, omega = omega, mins = mins,
       max_psrf_omega = max(psrf, na.rm = TRUE))
}
score <- function(h) c(cor_scale = stats::cor(h$cor[off], true_cor[off]),
                       cov_scale = stats::cor(h$omega[off], LL[off]),
                       mean_abs_omega = mean(abs(h$omega[off])),
                       max_psrf = h$max_psrf_omega, minutes = h$mins)

set.seed(1)
setups <- list(default = list(), nf2 = list(nf = 2),
               nf2_long = list(nf = 2, samples = 2000, transient = 5000, thin = 10))
oracle <- lapply(setups, function(a) do.call(fit_hmsc, c(list(y = sim$truth$z), a)))
tab <- t(sapply(oracle, score))
print(round(tab, 3))

# the setup with the best covariance-scale recovery, on the pilot's 20 draws
best <- names(which.max(tab[, "cov_scale"]))
# fp_occu_ms() sorts the species as text (sp1, sp10, sp2, ...)
zd <- readRDS("dev/jsdm-export/z_draws.rds")$zd[1:20, , paste0("sp", seq_len(J)),
                                                 drop = FALSE]
draws <- lapply(1:20, function(k) do.call(fit_hmsc, c(list(y = zd[k, , ]), setups[[best]])))
mean_omega <- Reduce(`+`, lapply(draws, `[[`, "omega")) / 20
naive <- do.call(fit_hmsc, c(list(y = (d$y > 0) * 1L), setups[[best]]))
res <- rbind(oracle = c(cov_scale = stats::cor(oracle[[best]]$omega[off], LL[off]),
                        mean_abs = mean(abs(oracle[[best]]$omega[off]))),
             naive = c(stats::cor(naive$omega[off], LL[off]), mean(abs(naive$omega[off]))),
             draws_mean = c(stats::cor(mean_omega[off], LL[off]), mean(abs(mean_omega[off]))))
print(round(res, 3))
saveRDS(list(setups = tab, best = best, draws = res,
             true_mean_abs = mean(abs(LL[off]))), "dev/jsdm-export/hmsc_setup.rds")
