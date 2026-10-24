# Figures of the random_det scenarios (dev/random-det/simulate-scenarios.R):
# coverage and bias of the occupancy slope, recovery of the two site-effect
# standard deviations, and convergence, for each scenario and model.
# Run from the package root: Rscript dev/random-det/plot-scenarios.R
# Output: dev/random-det/scenarios-*.png
library(ggplot2)
res <- readRDS(file.path("dev", "random-det", "scenarios.rds"))
res <- res[!is.na(res$slope), ]
res$scenario <- factor(sprintf("sd det = %g, sd fp = %g", res$sigma_det, res$sigma_fp))
res$model <- factor(res$model, levels = c("base", "random_det", "random_fp", "both"))
out <- function(name) file.path("dev", "random-det", paste0("scenarios-", name, ".png"))

# per scenario and model: mean and Monte Carlo standard error
summ <- purrr::map(split(res, res[c("scenario", "model")], drop = TRUE), function(g) {
  n <- nrow(g)
  err <- g$slope + 0.8
  cover <- mean(g$covered)
  data.frame(scenario = g$scenario[1], model = g$model[1], n = n,
             coverage = cover, coverage_se = sqrt(cover * (1 - cover) / n),
             bias = mean(err), bias_se = stats::sd(err) / sqrt(n),
             rhat = mean(g$max_rhat > 1.01), divergent = mean(g$divergences > 0),
             minutes = stats::median(g$seconds) / 60)
}) |> purrr::list_rbind()

theme_set(theme_bw(base_size = 11))
dodge <- position_dodge(width = 0.5)

p <- ggplot(summ, aes(model, coverage, colour = model)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_pointrange(aes(ymin = coverage - 2 * coverage_se,
                      ymax = pmin(1, coverage + 2 * coverage_se))) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(x = NULL, y = "Coverage of the 95% interval of the occupancy slope") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(out("coverage"), p, width = 9, height = 3.5, dpi = 150)

p <- ggplot(summ, aes(model, bias, colour = model)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_pointrange(aes(ymin = bias - 2 * bias_se, ymax = bias + 2 * bias_se)) +
  facet_wrap(~ scenario, nrow = 1) +
  labs(x = NULL, y = "Bias of the posterior mean of the occupancy slope") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(out("bias"), p, width = 9, height = 3.5, dpi = 150)

# estimated standard deviations of the site effects against the truth
sds <- rbind(
  data.frame(res[c("scenario", "model")], effect = "detection (sd_det_site)",
             estimate = res$sd_det_est, truth = res$sigma_det),
  data.frame(res[c("scenario", "model")], effect = "false positives (sd_fp_site)",
             estimate = res$sd_fp_est, truth = res$sigma_fp))
# each standard deviation only for the models that estimate it (the
# prototype writes sd_det = 0 when the effect is off)
fits_it <- ifelse(grepl("^detection", sds$effect),
                  sds$model %in% c("random_det", "both"),
                  sds$model %in% c("random_fp", "both"))
sds <- sds[fits_it & !is.na(sds$estimate), ]
p <- ggplot(sds, aes(model, estimate, colour = model)) +
  geom_hline(aes(yintercept = truth), linetype = "dashed") +
  geom_boxplot(outlier.size = 0.6) +
  facet_grid(effect ~ scenario) +
  labs(x = NULL, y = "Posterior mean of the standard deviation") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(out("sd"), p, width = 9, height = 5, dpi = 150)

conv <- rbind(data.frame(summ[c("scenario", "model")], what = "R-hat > 1.01", share = summ$rhat),
              data.frame(summ[c("scenario", "model")], what = "any divergence", share = summ$divergent))
p <- ggplot(conv, aes(model, share, fill = model)) +
  geom_col() +
  facet_grid(what ~ scenario) +
  labs(x = NULL, y = "Share of fits") +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(out("convergence"), p, width = 9, height = 4.5, dpi = 150)

print(summ, digits = 2, row.names = FALSE)
