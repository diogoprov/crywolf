# Grid for the export of occupancy draws to gllvm (README.md): 3 numbers of
# sites x 2 strengths of the latent factors x 5 replicates = 30 datasets
# from the JSDM simulator (10 species, 10 occasions, two latent factors,
# verified data at 30% of the sites). For each, fp_occu_ms() is fitted, 20
# state draws are taken, and gllvm estimates LL' from the true states, the
# naive states and each draw (averaged over draws). One file per dataset in
# dev/jsdm-export/grid/, so the grid resumes where it stopped.
#   JOB=1 NJOBS=2 Rscript dev/jsdm-export/grid.R   # jobs 1..NJOBS split the grid
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/jsdm-export/simulate_jsdm_fp.R")
source("dev/jsdm-export/z_draws.R")
job <- as.integer(Sys.getenv("JOB", "1"))
njobs <- as.integer(Sys.getenv("NJOBS", "1"))
ndraws <- 20
smoke <- nzchar(Sys.getenv("SMOKE"))   # a short check of the script
out_dir <- if (smoke) file.path(tempdir(), "grid-smoke") else "dev/jsdm-export/grid"
dir.create(out_dir, showWarnings = FALSE)

grid <- expand.grid(rep = 1:5, load_sd = c(0.5, 1), S = c(50, 100, 200))
grid$seed <- 2000 + seq_len(nrow(grid))
# largest datasets first, so that the two jobs end together
grid <- grid[order(-grid$S, grid$seed), ]
mine <- grid[seq_len(nrow(grid)) %% njobs == job %% njobs, ]
if (smoke) mine <- data.frame(rep = 1, load_sd = 1, S = 50, seed = 1)

run_one <- function(S, load_sd, seed) {
  sim <- simulate_jsdm_fp(S = S, K = 2, load_sd = load_sd, seed = seed)
  d <- sim$data
  J <- d$J
  sp_names <- paste0("sp", seq_len(J))
  set.seed(seed)
  occ_rows <- function(n, y) {
    purrr::map_dfr(seq_len(S), function(i) purrr::map_dfr(seq_len(J), function(j) {
      if (n[i, j] == 0) return(NULL)
      data.frame(site = i, occasion = seq_len(n[i, j]),
                 y = sample(rep(c(1L, 0L), c(y[i, j], n[i, j] - y[i, j]))),
                 species = sp_names[j])
    }))
  }
  fd <- fp_data(auto = occ_rows(d$n, d$y),
                sites = data.frame(site = seq_len(S), x = d$x),
                verified = occ_rows(d$nv, d$v), species = "species")
  t0 <- Sys.time()
  fit <- fp_occu_ms(fd, occ = ~ x, det = ~ 1, fp = ~ 1, random_recorder = FALSE,
                    chains = 4, parallel_chains = 4, threads_per_chain = 1,
                    seed = seed, refresh = 0,
                    iter_warmup = if (smoke) 100 else 1000,
                    iter_sampling = if (smoke) 100 else 1000)
  mins_ms <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  diag <- fit$diagnostics
  zd <- z_draws(fit, ndraws = ndraws, seed = seed)
  # fp_occu_ms() sorts the species as text (sp1, sp10, sp2, ...)
  zd <- zd[, , sp_names, drop = FALSE]
  rm(fit)

  X <- data.frame(x = d$x)
  ll_gllvm <- function(y) {
    colnames(y) <- sp_names
    y <- y[, colSums(y) > 0 & colSums(y) < S, drop = FALSE]
    m <- tryCatch(gllvm::gllvm(y, X = X, formula = ~ x,
                               family = binomial(link = "logit"), num.lv = 2,
                               sd.errors = FALSE, seed = 1),
                  error = function(e) NULL)
    ll <- matrix(NA_real_, J, J, dimnames = list(sp_names, sp_names))
    if (is.null(m)) return(ll)
    lam <- m$params$theta %*% diag(m$params$sigma.lv, 2)
    ll[colnames(y), colnames(y)] <- lam %*% t(lam)
    ll
  }
  off <- lower.tri(sim$truth$LL)
  score <- function(ll) c(cor = suppressWarnings(stats::cor(ll[off], sim$truth$LL[off],
                                                            use = "complete.obs")),
                          mean_abs = mean(abs(ll[off]), na.rm = TRUE))
  ll_draws <- lapply(seq_len(ndraws), function(k) ll_gllvm(zd[k, , ]))
  ok <- !vapply(ll_draws, function(l) all(is.na(l)), logical(1))
  ll_mi <- Reduce(`+`, ll_draws[ok]) / sum(ok)
  scores <- rbind(oracle = score(ll_gllvm(sim$truth$z)),
                  naive = score(ll_gllvm((d$y > 0) * 1L)),
                  draws_mean = score(ll_mi),
                  draw_single_median = apply(sapply(ll_draws[ok], score), 1,
                                             stats::median))
  list(S = S, load_sd = load_sd, seed = seed, scores = scores,
       true_mean_abs = mean(abs(sim$truth$LL[off])), draws_ok = sum(ok),
       mins_fp_occu_ms = mins_ms, diagnostics = diag)
}

for (k in seq_len(nrow(mine))) {
  g <- mine[k, ]
  f <- file.path(out_dir, sprintf("S%d_load%.1f_seed%d.rds", g$S, g$load_sd, g$seed))
  if (file.exists(f)) next
  cat(format(Sys.time()), "start", basename(f), "\n")
  res <- run_one(g$S, g$load_sd, g$seed)
  saveRDS(res, f)
  cat(format(Sys.time()), "done", basename(f), sprintf("(%.1f min)", res$mins_fp_occu_ms),
      "\n")
}
