# Compilation and caching of the Stan model -----------------------------------

.fp_env <- new.env(parent = emptyenv())

#' Compile (or load) the Stan model
#'
#' Loads the false-positive occupancy model compiled when the package was
#' installed from source. A binary package, or an installation without
#' CmdStan, comes without it: the model is then compiled the first time it
#' is needed and kept in the user cache directory
#' (`tools::R_user_dir("crywolf", "cache")`), so later sessions load it
#' without recompiling. The cache holds only the Stan program and its
#' executable, and a new version of the program replaces them. Under
#' `R CMD check` the cache is in [tempdir()], so nothing is written to the
#' user's filespace. [fp_occu()] calls this function itself; call it
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
#' @examplesIf instantiate::stan_cmdstan_exists()
#' \donttest{
#' fp_compile()
#' }
fp_compile <- function(force = FALSE, quiet = TRUE) {
  .compile_model("fp_occu", force = force, quiet = quiet)
}

#' where a model without its installed executable is compiled: the user cache,
#' or the session's temporary directory under R CMD check, which must not
#' write to the user's filespace
#' @noRd
.cache_dir <- function() {
  if (nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))) {
    return(file.path(tempdir(), "crywolf-cache"))
  }
  tools::R_user_dir("crywolf", which = "cache")
}

#' Compiles (or loads) one of the Stan programs of the package, `fp_occu` or
#' `fp_occu_dyn`; both #include fp_emissions.stanfunctions
#' @noRd
.compile_model <- function(name, force = FALSE, quiet = TRUE) {
  if (!force && !is.null(.fp_env[[name]])) return(invisible(.fp_env[[name]]))
  .check_cmdstan()
  # the executable compiled when the package was installed from source (if
  # it is missing, stan_package_model() would compile one without threads)
  src <- .stan_file(name)
  exe <- file.path(dirname(src),
                   if (.Platform$OS.type == "windows") paste0(name, ".exe") else name)
  mod <- if (!force && file.exists(exe)) {
    tryCatch(instantiate::stan_package_model(name = name, package = "crywolf"),
             error = function(e) NULL)
  }
  if (is.null(mod) || !.has_threads(mod)) {
    # binary packages (such as those on CRAN), or an installation without a
    # recent enough CmdStan, come without it: compile once into the user cache
    dir <- .cache_dir()
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    # the program and the functions it includes; a change to either
    # recompiles
    files <- c(src, file.path(dirname(src), "fp_emissions.stanfunctions"))
    changed <- FALSE
    for (f in files) {
      dst_f <- file.path(dir, basename(f))
      if (!file.exists(dst_f) || !identical(readLines(f), readLines(dst_f))) {
        file.copy(f, dst_f, overwrite = TRUE)
        changed <- TRUE
      }
    }
    dst <- file.path(dir, basename(src))
    threads <- list(stan_threads = TRUE)
    mod <- suppressWarnings(
      cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = force || changed,
                              cpp_options = threads)
    )
    if (!.has_threads(mod)) {
      mod <- cmdstanr::cmdstan_model(dst, quiet = quiet, force_recompile = TRUE,
                                     cpp_options = threads)
    }
  }
  .fp_env[[name]] <- .threads_option(mod)
  invisible(.fp_env[[name]])
}

#' cmdstanr 0.9.0 takes the options of an executable it loads from the
#' executable (as STAN_THREADS) but checks threads_per_chain against
#' `cpp_options()[["stan_threads"]]`, so it drops the threads with a warning;
#' add the lower-case name to a model compiled with threads
#' @noRd
.threads_option <- function(mod) {
  if (.has_threads(mod) && !isTRUE(mod$cpp_options()[["stan_threads"]])) {
    tryCatch(mod$.__enclos_env__$private$cpp_options_$stan_threads <- TRUE,
             error = function(e) NULL)
  }
  mod
}

#' Path to the Stan program
#'
#' @return The path to the Stan file shipped with the package, for users who
#'   want to read the model or modify a copy of it. The likelihood of the
#'   observations of each site is in `fp_emissions.stanfunctions`, in the
#'   same folder, which the program includes.
#' @export
#' @examples
#' fp_stan_file()
fp_stan_file <- function() {
  .stan_file("fp_occu")
}

#' Path of a Stan program of the package, installed or in the sources
#' @noRd
.stan_file <- function(name) {
  file <- paste0(name, ".stan")
  path <- system.file("bin", "stan", file, package = "crywolf")
  # under devtools::load_all() the program is still in the sources
  if (!nzchar(path)) {
    path <- system.file("src", "stan", file, package = "crywolf", mustWork = TRUE)
  }
  path
}

#' Whether a compiled model was built with threading
#' @noRd
.has_threads <- function(mod) {
  opts <- mod$cpp_options()
  v <- unlist(opts[tolower(names(opts)) == "stan_threads"])
  isTRUE(as.logical(v[1]))
}

.min_cmdstan <- "2.40.0"

#' the version of CmdStan that cmdstanr finds, or NULL
#' @noRd
.cmdstan_found <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(NULL)
  path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(path)) return(NULL)
  tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) NULL)
}

#' Whether a CmdStan version is recent enough for crywolf
#' @noRd
.cmdstan_recent <- function(v) utils::compareVersion(v, .min_cmdstan) >= 0

#' a CmdStan that crywolf can use: found, and recent enough (tests and
#' vignettes skip the fits otherwise)
#' @noRd
.has_cmdstan <- function() {
  v <- .cmdstan_found()
  !is.null(v) && .cmdstan_recent(v)
}

#' Stop unless a recent enough CmdStan is installed
#' @noRd
.check_cmdstan <- function() {
  v <- .cmdstan_found()
  if (is.null(v)) {
    stop("CmdStan was not found. Install it once with ",
         "`cmdstanr::install_cmdstan()` and try again.", call. = FALSE)
  }
  if (!.cmdstan_recent(v)) {
    stop("crywolf needs CmdStan ", sub("[.]0$", "", .min_cmdstan),
         " or later (found ", v, "). Update it ",
         "with `cmdstanr::install_cmdstan()`.", call. = FALSE)
  }
  invisible(TRUE)
}
