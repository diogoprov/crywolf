#' srr_stats
#'
#' Standards not yet met. They are tagged next to the code that addresses
#' them once done.
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstatsTODO {G5.4c} *Where applicable, stored values may be drawn from published paper outputs when applicable and where code from original implementations is not available.* Piatti et al. (2026) fitted the same model to the Pauraque data with an independent Stan program, at the level of detection windows; a test against their published estimates is pending.
#' @noRd
NULL

#' NA_standards
#'
#' @srrstatsNA {G2.14c} Missing values are not imputed: a missing detection is an unsurveyed occasion, which carries no information on detection, and missing covariates arise in ways the model cannot know, so the user imputes or drops them before fitting.
#' @srrstatsNA {G3.1, G3.1a} crywolf computes no covariance matrices; posterior summaries come from 'posterior'.
#' @srrstatsNA {G4.0} No function writes files named by the user; CmdStan output files and the compiled model are named and stored by the package.
#' @srrstatsNA {G5.11, G5.11a} The extended tests need no downloaded data: they simulate their data.
#' @srrstatsNA {BS1.0} The term "hyperparameter" is not used; the standard deviations of the random effects are named as such.
#' @srrstatsNA {BS1.3b} Only Stan's No-U-Turn sampler is offered: the discrete latent states are marginalised for Hamiltonian Monte Carlo, and other algorithms are not supported.
#' @srrstatsNA {BS1.5, BS5.4} R-hat, bulk and tail ESS, divergences, tree depth and E-BFMI are complementary diagnostics, always computed together, not alternative convergence checkers.
#' @srrstatsNA {BS2.4} Each prior parameter is one number shared by all coefficients of its class, so there are no vectors to match against the design.
#' @srrstatsNA {BS2.10} A single seed is accepted (checked by `.check_sampler()`), and CmdStan derives a distinct stream for each chain, so identical seeds cannot be given to different chains.
#' @srrstatsNA {BS2.11} crywolf has no starting-value argument of its own; starting values pass through `...` as cmdstanr's `init`.
#' @srrstatsNA {BS3.2} Perfectly collinear designs are diagnosed (BS3.1) but still sampled: the priors keep the posterior proper, and the posterior of a redundant coefficient follows its prior, as the warning says.
#' @srrstatsNA {BS4.4, BS4.6, BS4.7} Stan runs a fixed number of iterations and cannot stop on convergence; convergence is assessed after sampling with fixed thresholds.
#' @srrstatsNA {BS6.5} Optional standard ("may"); trace plots and posterior intervals are separate types of `plot()`.
#' @noRd
NULL
