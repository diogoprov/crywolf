# Conditions raised by the methods and outputs of fits, fp_occu_mc(),
# fp_draws_data(), fp_read_sheepdog() and the CmdStan checks. Most use small
# fake objects with only the fields read before the condition, or mocked
# internals, so that no model is fitted.

# a small export with calibration draws (as in test-fp_draws_data.R): two
# sites, one species, a grid of scores from 0.1 in steps of 0.1, two draws
write_cond_export <- function(dir) {
  dir.create(dir, showWarnings = FALSE)
  w <- function(d, name) utils::write.csv(d, file.path(dir, name), row.names = FALSE, na = "")
  w(data.frame(site = c("A", "A", "B", "B"), occasion_id = c("a1", "a2", "b1", "b2")),
    "occasions.csv")
  w(data.frame(site = c("A", "B"), occasion_id = c("a1", "b1"),
               species_code = "Sp one_sp-one", species = "sp_one",
               file_offset = 0, confidence = c(0.92, 0.5)),
    "detections.csv")
  w(data.frame(site = c("A", "B"), x = c(-1, 1)), "sites.csv")
  grid <- expand.grid(grid_index = 1:9, site = c("A", "B"), draw = c(3L, 8L),
                      stringsAsFactors = FALSE)
  grid$p_tp <- 1
  grid$species <- "sp_one"
  grid$confidence <- 0.1 * grid$grid_index
  arrow::write_parquet(grid, file.path(dir, "calibration_draws.parquet"))
  jsonlite::write_json(list(schema = "sheepdog-occupancy", schema_version = 1L,
                            sheepdog_version = "3.2.1.9000", occasion_unit = "file",
                            min_confidence = 0.1,
                            calibration = list(engine = "brms", n_draws = 2L,
                                               grid_step = 0.1, grid_min = 0.1)),
                       file.path(dir, "schema.json"), auto_unbox = TRUE)
  dir
}

# a fake multispecies fit with only the fields given
fake_ms <- function(...) structure(list(...), class = "fp_occu_ms")

#' @srrstats {G5.2b} Conditions raised by the methods and outputs of fits, fp_occu_mc(), fp_draws_data(), fp_read_sheepdog() and the CmdStan checks.
test_that("fp_draws_data() checks the export and the occasions", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  expect_error(fp_draws_data(list(occasions = data.frame())),
               "`x` must be an export read by fp_read_sheepdog().", fixed = TRUE)

  x <- fp_read_sheepdog(write_cond_export(tempfile("cond-export")))
  expect_error(fp_draws_data(x, occasions = rbind(x$occasions, x$occasions[1, ])),
               "`occasions` must have one row per site and occasion_id.", fixed = TRUE)

  # the schema announces a key the parquet does not have
  x$schema$calibration$keys <- list("site", "species", "recorder", "grid_index")
  expect_error(fp_draws_data(x),
               "disagree on the keys of the calibration: the parquet lacks recorder.",
               fixed = TRUE)
})

test_that("fp_draws_data() refuses a species column with several species", {
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("arrow")
  dir <- test_path("fixtures", "sheepdog-export-1.3-ms")
  schema <- jsonlite::read_json(file.path(dir, "schema.json"), simplifyVector = TRUE)
  x <- suppressWarnings(fp_read_sheepdog(dir, species = schema$species))
  occ <- x$occasions
  occ$species <- "lep_fus"
  expect_error(fp_draws_data(x, occasions = occ),
               "`occasions` has a column `species`; with several species it is added by fp_draws_data().",
               fixed = TRUE)
})

test_that("fp_read_sheepdog() needs occasions.csv and warns on mixed species", {
  skip_if_not_installed("jsonlite")
  dir <- tempfile("cond-schema")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE))
  jsonlite::write_json(list(schema = "sheepdog-occupancy", schema_version = 1L),
                       file.path(dir, "schema.json"), auto_unbox = TRUE)
  expect_error(fp_read_sheepdog(dir), "The export has no occasions.csv.", fixed = TRUE)

  # the validation of this export holds several species (the clips drawn
  # for their high scores raise a warning of their own)
  ms <- test_path("fixtures", "sheepdog-export-1.3-ms")
  warnings_of <- function(expr) {
    w <- character()
    withCallingHandlers(expr, warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    })
    w
  }
  w <- warnings_of(fp_read_sheepdog(ms))
  expect_match(w, paste("The export holds several species, so `validation` mixes the verdicts",
                        "of all of them. Name the target species in `species`"),
               fixed = TRUE, all = FALSE)
  schema <- jsonlite::read_json(file.path(ms, "schema.json"), simplifyVector = TRUE)
  w <- warnings_of(fp_read_sheepdog(ms, species = schema$species[1]))
  expect_false(any(grepl("several species", w)))
})

test_that("fp_occu_mc() checks max_rhat before fitting", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  msg <- "`max_rhat` must be a single number of at least 1 (Inf keeps every fit)."
  expect_error(fp_occu_mc(list(s$data), max_rhat = NA), msg, fixed = TRUE)
  expect_error(fp_occu_mc(list(s$data), max_rhat = 0.9), msg, fixed = TRUE)
  expect_error(fp_occu_mc(list(s$data), max_rhat = c(1.01, 1.1)), msg, fixed = TRUE)
  expect_error(fp_occu_mc(list(s$data), max_rhat = "1.05"), msg, fixed = TRUE)
})

test_that("fp_occu_mc() warns when the fits have different numbers of draws", {
  s <- simulate_fp(nsite = 20, nocc = 3, seed = 1)
  n_chains <- c(1L, 2L)
  fit_no <- 0L
  local_mocked_bindings(
    fp_compile = function(...) invisible(NULL),
    # fit 1: one chain of 10 draws; fit 2: two chains (a failed chain in fit 1)
    fp_occu = function(...) {
      fit_no <<- fit_no + 1L
      nc <- n_chains[fit_no]
      draws <- posterior::as_draws_array(array(stats::rnorm(10 * nc), c(10, nc, 1),
                                               dimnames = list(NULL, NULL, "a")))
      structure(list(draws = draws,
                     diagnostics = list(divergences = 0, max_rhat = 1,
                                        min_ess_bulk = 100, min_ess_tail = 100)),
                class = "fp_occu")
    },
    summary.fp_occu = function(object, ...) list(table = data.frame(variable = "a"))
  )
  expect_warning(mc <- fp_occu_mc(list(s$data, s$data)),
                 "The fits have different numbers of draws (10 to 20)", fixed = TRUE)
  expect_equal(mc$diagnostics$n_draws, c(10, 20))
  expect_equal(posterior::nchains(mc$draws), 3)
})

test_that("multispecies outputs check their input", {
  expect_error(fp_richness(list()), "`object` must be an `fp_occu_ms` fit.", fixed = TRUE)
  expect_error(fp_richness(structure(list(), class = "fp_occu")),
               "`object` must be an `fp_occu_ms` fit.", fixed = TRUE)

  # no species are needed to reach the checks of the sites table
  local_mocked_bindings(.ms_species = function(object, ndraws) list())
  m <- fake_ms(data = list(sites = data.frame(site = "A", richness = 1, verified = 0)))
  expect_error(fp_richness(m),
               "The sites table already has column(s) richness, verified; rename them before calling fp_richness().",
               fixed = TRUE)

  m <- fake_ms(data = list(species = "species",
                           sites = data.frame(site = "A", species = "x")))
  expect_error(fp_sites(m),
               "The sites table already has a column `species`, the species column; rename it before calling fp_sites().",
               fixed = TRUE)

  m <- fake_ms(data = list(), stan_data = list(N1 = 0))
  expect_error(fp_ppc(m, method = "verified"),
               "The data of this multispecies fit have no verified occasions.", fixed = TRUE)

  m <- fake_ms(data = list(species = "species", species_names = c("a", "b")),
               formulas = list(occ = ~ 1, fp = NULL))
  expect_error(predict(m, process = "fp"),
               "This multispecies model has no false-positive process.", fixed = TRUE)
})

test_that("the parameters of a species must match the levels of the model", {
  # draws of one species with an intercept-only model
  dm <- matrix(0, 2, 4, dimnames = list(NULL, c("alpha[1]", "alpha_sp[1]",
                                                "beta_sp[1,1]", "delta0")))
  sd <- list(Kalpha = 1, Kfp = 0, Kbeta = 1, estimate_fp = 0, estimate_re = 0)
  object <- list(stan_data = sd,
                 design = list(species = "sp1", ver_names = c("(Intercept)", "a")))
  gen <- list(stan_data = list(Kalpha = 1, Kfp = 0, Kver_s = 1),
              design = list(ver_names = c("(Intercept)", "b")))
  expect_error(.ms_pars(dm, object, gen, 1),
               "Species sp1 has levels of `delta_s` the model does not have.", fixed = TRUE)
  gen$stan_data$Kalpha <- 2
  expect_error(.ms_pars(dm, object, gen, 1),
               "The rows of species sp1 do not have every level of the covariates of `det` or `fp`",
               fixed = TRUE)
})

test_that("predict() and fp_loo_weights() check their input", {
  m <- structure(list(formulas = list(occ = ~ 1, fp = NULL)), class = "fp_occu")
  expect_error(predict(m, process = "fp"),
               "This model has no false-positive process to predict.", fixed = TRUE)
  msg <- "Name every model, e.g. fp_loo_weights(constant = m0, fp_x = m1)."
  expect_error(fp_loo_weights(1, 2), msg, fixed = TRUE)
  expect_error(fp_loo_weights(a = 1, 2), msg, fixed = TRUE)
})

test_that("loo(reloo = TRUE) announces the refits and keeps sites it cannot refit", {
  set.seed(1)
  # log-likelihood of three sites: 50 iterations x 2 chains x 3 sites
  ll <- array(stats::rnorm(300, -2, 0.3), c(50, 2, 3),
              dimnames = list(NULL, NULL, paste0("log_lik[", 1:3, "]")))
  auto <- data.frame(site = c("A", "B", "C"), y = 0L, recorder = c("r1", "r1", "r2"))
  x <- structure(list(
    fit = list(draws = function(v) posterior::as_draws_array(ll)),
    data = list(auto = auto, sites = data.frame(site = c("A", "B", "C")),
                site = "site", y = "y", recorder = "recorder"),
    formulas = list(occ = ~ 1, fp = ~ 1),
    design = list(recorders = c("r1", "r2"))), class = "fp_occu")
  # every refit lacks recorder r2, as if each site had a recorder of its own
  local_mocked_bindings(
    fp_data = function(...) NULL,
    fp_occu = function(...) list(design = list(recorders = "r1"))
  )
  w <- character()
  expect_message(
    l <- withCallingHandlers(
      loo(x, reloo = TRUE, k_threshold = -Inf),
      warning = function(e) {
        w <<- c(w, conditionMessage(e))
        invokeRestart("muffleWarning")
      }),
    "Refitting the model without each of 3 site(s) with Pareto k > -Inf.", fixed = TRUE)
  expect_equal(w, sprintf(paste("Site %s is the only site with its recorder or with a level",
                                "of a factor covariate; its approximate value is kept."),
                          c("A", "B", "C")))
  # every site keeps its approximate value
  expect_identical(l, loo(x))
})

test_that("the CmdStan checks report a missing or old installation", {
  local_mocked_bindings(.cmdstan_found = function() NULL)
  expect_error(fp_compile(force = TRUE),
               "CmdStan was not found. Install it once with `cmdstanr::install_cmdstan()` and try again.",
               fixed = TRUE)
  expect_false(.has_cmdstan())
  local_mocked_bindings(.cmdstan_found = function() "2.33.1")
  expect_error(fp_compile(force = TRUE),
               "crywolf needs CmdStan 2.40 or later (found 2.33.1). Update it with `cmdstanr::install_cmdstan()`.",
               fixed = TRUE)
  expect_false(.has_cmdstan())
  local_mocked_bindings(.cmdstan_found = function() "2.40.0")
  expect_true(.check_cmdstan())
})
