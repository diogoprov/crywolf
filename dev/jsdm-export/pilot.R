# Feasibility of exporting occupancy draws from fp_occu_ms() to gllvm and
# Hmsc (README.md). One dataset from the JSDM simulator (100 sites, 10
# species, 10 occasions, two latent factors, verified data at 30% of sites).
# The residual covariance LL' is estimated by gllvm (two latent variables,
# logit link) from: the true states (oracle), the naive states (at least one
# automated detection), and each of NDRAWS posterior draws of the states
# (multiple imputation, estimates averaged over draws). Hmsc is fitted to the
# true states and to ten draws.
#   NDRAWS=20 Rscript dev/jsdm-export/pilot.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/jsdm-export/simulate_jsdm_fp.R")   # copy of jsdm-proto:dev/jsdm/
source("dev/jsdm-export/z_draws.R")
ndraws <- as.integer(Sys.getenv("NDRAWS", "20"))
out_file <- Sys.getenv("OUT", "dev/jsdm-export/pilot.rds")
seed <- 1001

sim <- simulate_jsdm_fp(K = 2, seed = seed)
d <- sim$data
S <- d$S; J <- d$J
sp_names <- paste0("sp", seq_len(J))
# occasion-level rows from the aggregated counts (occasions are exchangeable)
set.seed(seed)
occ_rows <- function(n, y) {
  purrr::map_dfr(seq_len(S), function(i) purrr::map_dfr(seq_len(J), function(j) {
    if (n[i, j] == 0) return(NULL)
    data.frame(site = i, occasion = seq_len(n[i, j]),
               y = sample(rep(c(1L, 0L), c(y[i, j], n[i, j] - y[i, j]))),
               species = sp_names[j])
  }))
}
auto <- occ_rows(d$n, d$y)
verified <- occ_rows(d$nv, d$v)
sites <- data.frame(site = seq_len(S), x = d$x)
fd <- fp_data(auto = auto, sites = sites, verified = verified, species = "species")

# the fit takes about 15 min: its state draws are kept for reruns
z_file <- "dev/jsdm-export/z_draws.rds"
if (file.exists(z_file)) {
  zz <- readRDS(z_file)
  zd <- zz$zd
  mins_ms <- zz$mins
} else {
  t0 <- Sys.time()
  fit <- fp_occu_ms(fd, occ = ~ x, det = ~ 1, fp = ~ 1, random_recorder = FALSE,
                    chains = 4, parallel_chains = 4, threads_per_chain = 1,
                    seed = seed, refresh = 0)
  mins_ms <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  zd <- z_draws(fit, ndraws = 200, seed = seed)
  saveRDS(list(zd = zd, mins = mins_ms), z_file)
}
zd <- zd[seq_len(ndraws), , , drop = FALSE]

X <- data.frame(x = d$x)
ll_gllvm <- function(y) {
  y <- y[, colSums(y) > 0 & colSums(y) < S, drop = FALSE]
  m <- gllvm::gllvm(y, X = X, formula = ~ x, family = binomial(link = "logit"),
                    num.lv = 2, sd.errors = FALSE, seed = 1)
  lam <- m$params$theta %*% diag(m$params$sigma.lv, 2)
  ll <- matrix(NA, J, J, dimnames = list(sp_names, sp_names))
  ll[colnames(y), colnames(y)] <- lam %*% t(lam)
  ll
}
off <- lower.tri(sim$truth$LL)
score <- function(ll) c(cor = stats::cor(ll[off], sim$truth$LL[off], use = "complete.obs"),
                        mean_abs = mean(abs(ll[off]), na.rm = TRUE))
truth_z <- sim$truth$z
naive_z <- (d$y > 0) * 1L
colnames(truth_z) <- colnames(naive_z) <- sp_names
t0 <- Sys.time()
ll_draws <- lapply(seq_len(ndraws), function(k) ll_gllvm(zd[k, , ]))
mins_gllvm <- as.numeric(difftime(Sys.time(), t0, units = "mins")) / ndraws
ll_mi <- Reduce(`+`, ll_draws) / ndraws
gl <- rbind(oracle = score(ll_gllvm(truth_z)), naive = score(ll_gllvm(naive_z)),
            draws_mean = score(ll_mi),
            draw_single_median = apply(sapply(ll_draws, score), 1, stats::median))

# Hmsc: the true states and ten draws
hm_ll <- function(y) {
  m <- Hmsc::Hmsc(Y = y, XData = X, XFormula = ~ x, distr = "probit",
                  studyDesign = data.frame(site = factor(seq_len(S))),
                  ranLevels = list(site = Hmsc::HmscRandomLevel(units = factor(seq_len(S)))))
  m <- Hmsc::sampleMcmc(m, samples = 1000, thin = 5, transient = 2500, nChains = 2,
                        nParallel = 2, verbose = 0)
  Hmsc::computeAssociations(m)[[1]]$mean   # residual correlation (probit scale)
}
t0 <- Sys.time()
hm_true <- hm_ll(truth_z)
mins_hmsc <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
hm_draws <- lapply(1:10, function(k) hm_ll(zd[k, , ]))
true_cor <- stats::cov2cor(sim$truth$LL)   # residual correlation of the factors
cor_vs <- function(r) stats::cor(r[off], true_cor[off])
hm <- c(oracle = cor_vs(hm_true), draws_mean = cor_vs(Reduce(`+`, hm_draws) / 10))

res <- list(gllvm = gl, hmsc_cor = hm, minutes = c(fp_occu_ms = mins_ms,
            gllvm_per_draw = mins_gllvm, hmsc_per_fit = mins_hmsc),
            true_mean_abs = mean(abs(sim$truth$LL[off])))
saveRDS(res, out_file)
print(res, digits = 3)
