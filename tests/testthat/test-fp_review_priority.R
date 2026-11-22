# .review_units() works on the table of fp_sites(), so most tests need no fit
fake_sites <- function(occupied, verified = rep(FALSE, length(occupied)), ...) {
  data.frame(plot = c(1, 2, 3, 4)[seq_along(occupied)], ..., occupied = occupied,
             verified = verified)
}

test_that("the size is the variance of the occupancy state, largest first", {
  st <- fake_sites(c(0.1, 0.5, 1, 0.9), verified = c(FALSE, FALSE, TRUE, FALSE))
  u <- crywolf:::.review_units(st, list(site = "plot"), "sp A", NULL)
  expect_named(u, c("site", "species_code", "size"))
  # the verified site is left out; the others ordered by q(1 - q)
  expect_identical(u$site, c("2", "1", "4"))
  expect_equal(u$size, c(0.25, 0.09, 0.09))
  expect_identical(unique(u$species_code), "sp A")
})

test_that("sites with a confirmed clip are left out", {
  st <- fake_sites(c(0.2, 0.5, 0.7))
  val <- data.frame(site = c("2", "3", "3"), species_code = "sp A",
                    verdict = c(FALSE, TRUE, NA))
  u <- crywolf:::.review_units(st, list(site = "plot"), "sp A", val)
  expect_identical(sort(u$site), c("1", "2"))
})

test_that("multispecies units are site x species, with labels mapped", {
  st <- fake_sites(rep(c(0.5, 0.2), 2), species = c("a", "b", "a", "b"))
  st$plot <- c(1, 1, 2, 2)
  data <- list(site = "plot", species = "species", species_names = c("a", "b"))
  u <- crywolf:::.review_units(st, data, NULL, NULL)
  expect_identical(u$species_code, c("a", "a", "b", "b"))
  map <- c(a = "Species a_A", b = "Species b_B")
  val <- data.frame(site = "1", species_code = "Species a_A", verdict = TRUE)
  u <- crywolf:::.review_units(st, data, map, val)
  expect_identical(paste(u$site, u$species_code),
                   c("2 Species a_A", "1 Species b_B", "2 Species b_B"))
})

#' @srrstats {G5.2, G5.2a, G5.2b} Every condition of fp_review_priority() is triggered and its text checked.
test_that("fp_review_priority() checks its arguments", {
  expect_error(fp_review_priority(list()),
               "`object` must be an fp_occu or fp_occu_ms model.", fixed = TRUE)
  fake <- structure(list(data = list(site = "plot")), class = "fp_occu")
  expect_error(fp_review_priority(fake),
               "`species_code` must be the label of the species in the predictions of sheepdog, one string.",
               fixed = TRUE)
  check <- crywolf:::.check_review_args
  data <- list(site = "plot", species = "species", species_names = c("a", "b"))
  expect_error(check(data, c(a = "A"), NULL),
               "`species_code` must be a character vector named by the species of the model: a, b.",
               fixed = TRUE)
  expect_error(check(data, NULL, data.frame(site = "1")),
               "`validation` lacks the column(s) verdict, species_code.", fixed = TRUE)
  val <- data.frame(site = "1", species_code = c("x", "y"), verdict = TRUE)
  expect_error(check(list(site = "plot"), "x", val),
               "`validation` holds clips of several species; keep those of the species of the model (fp_read_sheepdog(species = )).",
               fixed = TRUE)
  expect_warning(check(data, NULL, val),
                 "No `species_code` of `validation` is a species label of the units (a, b), so no unit is left out as confirmed; map the species of the model to their labels with `species_code`.",
                 fixed = TRUE)
  expect_silent(check(data, c(a = "x", b = "y"), val))
})

test_that("fp_review_priority() takes the occupancy of fp_sites()", {
  skip_if_no_cmdstan()
  s <- simulate_fp(nsite = 60, nocc = 6, seed = 5)
  m <- suppressWarnings(fp_occu(s$data, occ = ~ x, fp = ~ 1, random_recorder = FALSE,
                                chains = 2, iter_warmup = 300, iter_sampling = 300,
                                threads_per_chain = 1, seed = 5, refresh = 0,
                                show_messages = FALSE))
  st <- fp_sites(m, ndraws = 100, seed = 1)
  u <- fp_review_priority(m, species_code = "sp", ndraws = 100, seed = 1)
  expect_equal(nrow(u), sum(!st$verified))
  i <- match(u$site, as.character(st[[m$data$site]]))
  expect_equal(u$size, st$occupied[i] * (1 - st$occupied[i]))
  expect_true(all(u$size >= 0 & u$size <= 0.25))
  expect_false(is.unsorted(rev(u$size)))
})
