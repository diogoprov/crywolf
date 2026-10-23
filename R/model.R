# Compilation and caching of the Stan model -----------------------------------

.fp_env <- new.env(parent = emptyenv())

#' Compile (or load) the Stan model
#'
#' Compiles the false-positive occupancy model the first time it is needed and
#' keeps the executable in the user cache directory
#' (`tools::R_user_dir("crywolf", "cache")`), so later sessions load it
#' without recompiling. [fp_occu()] calls this function itself; call it
#' directly to compile ahead of time, for example before launching parallel
#' workers in [fp_calibrate()].
#'
#' The model is compiled with threading enabled (`stan_threads`), so that
#' [fp_occu()] can split the likelihood of each chain among several threads
#' (argument `threads_per_chain`).
#'
#' @param force Logical. Recompile even if a current executable exists.
#' @param quiet Logical. Suppress compiler output.
#'
#' @return A `CmdStanModel` object (invisibly).
#' @export
#' @examples
#' \dontrun{
#' fp_compile()
#' }
fp_compile <- function(force = FALSE, quiet = TRUE) {
  if (!force && !is.null(.fp_env$model)) return(invisible(.fp_env$model))
  .check_cmdstan()
  src <- system.file("stan", "fp_occu.stan", package = "crywolf", mustWork = TRUE)
  dir <- tools::R_user_dir("crywolf", which = "cache")
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dst <- file.path(dir, "fp_occu.stan")
  # copy only when the source changed, so the executable stays up to date
  # without triggering a recompile in every session
  if (!file.exists(dst) || !identical(readLines(src), readLines(dst))) {
    file.copy(src, dst, overwrite = TRUE)
  }
  threads <- list(stan_threads = TRUE)
  mod <- suppressWarnings(
    cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = force,
                            cpp_options = threads)
  )
  # an executable built without threading (by an earlier version of the
  # package) is not rebuilt by cmdstanr unless asked
  if (!.has_threads(mod)) {
    mod <- cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = TRUE,
                                   cpp_options = threads)
  }
  .fp_env$model <- mod
  invisible(mod)
}

#' Path to the Stan program
#'
#' @return The path to the Stan file shipped with the package, for users who
#'   want to read the model or modify a copy of it.
#' @export
#' @examples
#' fp_stan_file()
fp_stan_file <- function() {
  system.file("stan", "fp_occu.stan", package = "crywolf", mustWork = TRUE)
}

.has_threads <- function(mod) {
  opts <- mod$cpp_options()
  v <- unlist(opts[tolower(names(opts)) == "stan_threads"])
  isTRUE(as.logical(v[1]))
}

.has_cmdstan <- function() {
  isTRUE(tryCatch(nzchar(cmdstanr::cmdstan_path()), error = function(e) FALSE))
}

.check_cmdstan <- function() {
  if (!.has_cmdstan()) {
    stop("CmdStan was not found. Install it once with ",
         "`cmdstanr::install_cmdstan()` and try again.", call. = FALSE)
  }
  v <- as.character(cmdstanr::cmdstan_version())
  if (utils::compareVersion(v, "2.40.0") < 0) {
    stop("crywolf needs CmdStan 2.40 or later (found ", v, "). Update it ",
         "with `cmdstanr::install_cmdstan()`.", call. = FALSE)
  }
  invisible(TRUE)
}
