#!/usr/bin/env Rscript
# =============================================================================
# Publication-Ready Power Analysis Figures for PPI++
# =============================================================================
# Generates comprehensive power curves and sample size curves for:
# - Continuous outcomes (varying R²/correlation)
# - Binary outcomes (varying sensitivity/specificity)
#
# Output: High-quality figures suitable for academic publication

library(pppower)
library(ggplot2)
library(gridExtra)
library(grid)

# Palette objects are exported by pppower (loaded above).
# Fallback: source from package root if running in dev mode.
if (!exists("pppower_method_colors")) {
  tryCatch(source("R/palette.R", local = FALSE), error = function(e) {
    tryCatch(source("../R/palette.R", local = FALSE), error = function(e2) NULL)
  })
}

# Aliases for backward compatibility
theme_pub <- theme_pppower
colors_method <- pppower_method_colors
colors_quality <- pppower_quality_colors

cat("=============================================================================\n")
cat("Generating Publication-Ready Power Analysis Figures\n")
cat("=============================================================================\n\n")

# =============================================================================
# PART 1: CONTINUOUS OUTCOMES
# =============================================================================
cat("Part 1: Continuous Outcomes (varying R²)\n")
cat("----------------------------------------\n")

# Parameters
N_cont <- 10000          # unlabeled samples
sigma_y <- 1             # outcome SD
alpha <- 0.05
power_target <- 0.80

# Classifier quality levels (R² = correlation²)
r2_levels <- list(
  "Poor (R²=0.25)" = 0.25,
  "Fair (R²=0.50)" = 0.50,
  "Good (R²=0.70)" = 0.70,
  "Excellent (R²=0.90)" = 0.90
)

# --- Power Curves for Continuous ---
cat("  Generating power curves...\n")

n_values_cont <- c(50, 100, 200, 500)
delta_values_cont <- seq(0, 0.6, by = 0.02)

power_cont_results <- expand.grid(
  n = n_values_cont,
  delta = delta_values_cont,
  r2_label = names(r2_levels)
)
power_cont_results$power_ppi <- NA
power_cont_results$power_vanilla <- NA
power_cont_results$power_classical <- NA

for (i in seq_len(nrow(power_cont_results))) {
  n <- power_cont_results$n[i]
  delta <- power_cont_results$delta[i]
  r2_label <- power_cont_results$r2_label[i]
  r2 <- r2_levels[[r2_label]]

  rho <- sqrt(r2)
  var_f <- r2 * sigma_y^2
  var_res <- (1 - r2) * sigma_y^2
  cov_yf <- r2 * sigma_y^2

  # PPI++ power (oracle lambda = EIF-optimal)
  power_cont_results$power_ppi[i] <- power_ppi_mean(
    delta = delta, N = N_cont, n = n, alpha = alpha,
    var_f = var_f, var_res = var_res, cov_y_f = cov_yf
  )

  # Vanilla PPI power (lambda = 1)
  power_cont_results$power_vanilla[i] <- power_ppi_mean(
    delta = delta, N = N_cont, n = n, alpha = alpha,
    var_f = var_f, var_res = var_res, cov_y_f = cov_yf,
    lambda = 1, lambda_type = "user"
  )

  # Classical power
  se_classical <- sigma_y / sqrt(n)
  z_alpha <- qnorm(1 - alpha/2)
  if (delta == 0) {
    power_cont_results$power_classical[i] <- alpha
  } else {
    mu <- abs(delta) / se_classical
    power_cont_results$power_classical[i] <- 1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
  }
}

# Order factor levels
power_cont_results$r2_label <- factor(power_cont_results$r2_label,
                                       levels = names(r2_levels))
power_cont_results$n_label <- factor(paste0("n = ", power_cont_results$n),
                                      levels = paste0("n = ", sort(n_values_cont)))

# --- Sample Size Curves for Continuous ---
cat("  Generating sample size curves...\n")

delta_values_ss <- seq(0.1, 0.6, by = 0.02)
ss_cont_results <- expand.grid(
  delta = delta_values_ss,
  r2_label = names(r2_levels)
)
ss_cont_results$n_ppi <- NA
ss_cont_results$n_vanilla <- NA
ss_cont_results$n_classical <- NA

z_alpha <- qnorm(1 - alpha/2)
z_beta <- qnorm(power_target)

for (i in seq_len(nrow(ss_cont_results))) {
  delta <- ss_cont_results$delta[i]
  r2_label <- ss_cont_results$r2_label[i]
  r2 <- r2_levels[[r2_label]]

  rho <- sqrt(r2)
  sigma_y2 <- sigma_y^2
  sigma_f2 <- r2 * sigma_y^2
  cov_yf <- rho * sigma_y * sqrt(sigma_f2)

  # PPI++ required n (oracle lambda = EIF-optimal)
  ss_cont_results$n_ppi[i] <- tryCatch(
    power_ppi_mean(
      delta = delta, N = N_cont, n = NULL, power = power_target, alpha = alpha,
      sigma_y2 = sigma_y2, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
      lambda_mode = "oracle"
    ),
    error = function(e) NA
  )

  # Vanilla PPI required n (lambda = 1)
  ss_cont_results$n_vanilla[i] <- tryCatch(
    power_ppi_mean(
      delta = delta, N = N_cont, n = NULL, power = power_target, alpha = alpha,
      sigma_y2 = sigma_y2, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
      lambda_mode = "vanilla"
    ),
    error = function(e) NA
  )

  # Classical required n
  ss_cont_results$n_classical[i] <- sigma_y^2 * ((z_alpha + z_beta) / delta)^2
}

ss_cont_results$r2_label <- factor(ss_cont_results$r2_label, levels = names(r2_levels))

# =============================================================================
# PART 2: BINARY OUTCOMES
# =============================================================================
cat("\nPart 2: Binary Outcomes (varying Sensitivity/Specificity)\n")
cat("----------------------------------------------------------\n")

N_bin <- 10000
prevalence <- 0.30
alpha <- 0.05

# Classifier quality levels
classifier_levels <- list(
  "Poor (70%/70%)" = c(sens = 0.70, spec = 0.70),
  "Fair (80%/80%)" = c(sens = 0.80, spec = 0.80),
  "Good (85%/90%)" = c(sens = 0.85, spec = 0.90),
  "Excellent (95%/95%)" = c(sens = 0.95, spec = 0.95)
)

# --- Power Curves for Binary ---
cat("  Generating power curves...\n")

n_values_bin <- c(50, 100, 200, 500)
delta_values_bin <- seq(0, 0.15, by = 0.005)

power_bin_results <- expand.grid(
  n = n_values_bin,
  delta = delta_values_bin,
  clf_label = names(classifier_levels)
)
power_bin_results$power_ppi <- NA
power_bin_results$power_vanilla <- NA
power_bin_results$power_classical <- NA

var_y_bin <- prevalence * (1 - prevalence)

for (i in seq_len(nrow(power_bin_results))) {
  n <- power_bin_results$n[i]
  delta <- power_bin_results$delta[i]
  clf_label <- power_bin_results$clf_label[i]
  clf <- classifier_levels[[clf_label]]

  # PPI++ power (oracle lambda = EIF-optimal)
  power_bin_results$power_ppi[i] <- power_ppi_mean(
    delta = delta, N = N_bin, n = n, alpha = alpha,
    metrics = list(
      sensitivity = clf["sens"], specificity = clf["spec"],
      p_y = prevalence, m_obs = n
    ),
    metric_type = "classification"
  )

  # Vanilla PPI power (lambda = 1)
  power_bin_results$power_vanilla[i] <- power_ppi_mean(
    delta = delta, N = N_bin, n = n, alpha = alpha,
    metrics = list(
      sensitivity = clf["sens"], specificity = clf["spec"],
      p_y = prevalence, m_obs = n
    ),
    metric_type = "classification",
    lambda = 1, lambda_type = "user"
  )

  # Classical power
  se_classical <- sqrt(var_y_bin / n)
  z_alpha <- qnorm(1 - alpha/2)
  if (delta == 0) {
    power_bin_results$power_classical[i] <- alpha
  } else {
    mu <- abs(delta) / se_classical
    power_bin_results$power_classical[i] <- 1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
  }
}

power_bin_results$clf_label <- factor(power_bin_results$clf_label,
                                       levels = names(classifier_levels))
power_bin_results$n_label <- factor(paste0("n = ", power_bin_results$n),
                                     levels = paste0("n = ", sort(n_values_bin)))

# --- Sample Size Curves for Binary ---
cat("  Generating sample size curves...\n")

delta_values_ss_bin <- seq(0.02, 0.15, by = 0.005)
ss_bin_results <- expand.grid(
  delta = delta_values_ss_bin,
  clf_label = names(classifier_levels)
)
ss_bin_results$n_ppi <- NA
ss_bin_results$n_vanilla <- NA
ss_bin_results$n_classical <- NA

for (i in seq_len(nrow(ss_bin_results))) {
  delta <- ss_bin_results$delta[i]
  clf_label <- ss_bin_results$clf_label[i]
  clf <- classifier_levels[[clf_label]]

  # PPI++ required n (oracle lambda = EIF-optimal)
  ss_bin_results$n_ppi[i] <- tryCatch(
    power_ppi_mean(
      delta = delta, N = N_bin, n = NULL, power = power_target, alpha = alpha,
      metrics = list(
        sensitivity = clf["sens"], specificity = clf["spec"],
        p_y = prevalence, m_obs = 200
      ),
      metric_type = "classification",
      lambda_mode = "oracle"
    ),
    error = function(e) NA
  )

  # Vanilla PPI required n (lambda = 1)
  ss_bin_results$n_vanilla[i] <- tryCatch(
    power_ppi_mean(
      delta = delta, N = N_bin, n = NULL, power = power_target, alpha = alpha,
      metrics = list(
        sensitivity = clf["sens"], specificity = clf["spec"],
        p_y = prevalence, m_obs = 200
      ),
      metric_type = "classification",
      lambda_mode = "vanilla"
    ),
    error = function(e) NA
  )

  # Classical required n
  ss_bin_results$n_classical[i] <- var_y_bin * ((z_alpha + z_beta) / delta)^2
}

ss_bin_results$clf_label <- factor(ss_bin_results$clf_label,
                                    levels = names(classifier_levels))

# =============================================================================
# PART 3: CREATE PUBLICATION FIGURES
# =============================================================================
cat("\nPart 3: Creating Publication Figures\n")
cat("-------------------------------------\n")

# --- Figure 1: Continuous Power Curves (by classifier quality) ---
cat("  Figure 1: Continuous power curves...\n")

# Reshape for plotting
power_cont_long <- rbind(
  data.frame(power_cont_results[, c("n", "delta", "r2_label", "n_label")],
             power = power_cont_results$power_ppi, method = "PPI++"),
  data.frame(power_cont_results[, c("n", "delta", "r2_label", "n_label")],
             power = power_cont_results$power_vanilla, method = "Vanilla PPI"),
  data.frame(power_cont_results[, c("n", "delta", "r2_label", "n_label")],
             power = power_cont_results$power_classical, method = "Classical")
)

fig1 <- ggplot(power_cont_long[power_cont_long$n == 200, ],
               aes(x = delta, y = power, color = method, linetype = method)) +
  geom_line(linewidth = 1.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", alpha = 0.5) +
  facet_wrap(~r2_label, nrow = 1) +
  scale_color_manual(values = colors_method) +
  scale_linetype_manual(values = pppower_method_linetypes) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A. Continuous Outcomes: Power Curves by Prediction Quality",
    subtitle = expression(paste("n = 200 labeled, N = 10,000 unlabeled, ", sigma[Y], " = 1")),
    x = "Effect Size (Delta)",
    y = "Power",
    color = "Method",
    linetype = "Method"
  ) +
  theme_pub()

# --- Figure 2: Binary Power Curves (by classifier quality) ---
cat("  Figure 2: Binary power curves...\n")

power_bin_long <- rbind(
  data.frame(power_bin_results[, c("n", "delta", "clf_label", "n_label")],
             power = power_bin_results$power_ppi, method = "PPI++"),
  data.frame(power_bin_results[, c("n", "delta", "clf_label", "n_label")],
             power = power_bin_results$power_vanilla, method = "Vanilla PPI"),
  data.frame(power_bin_results[, c("n", "delta", "clf_label", "n_label")],
             power = power_bin_results$power_classical, method = "Classical")
)

fig2 <- ggplot(power_bin_long[power_bin_long$n == 200, ],
               aes(x = delta, y = power, color = method, linetype = method)) +
  geom_line(linewidth = 1.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", alpha = 0.5) +
  facet_wrap(~clf_label, nrow = 1) +
  scale_color_manual(values = colors_method) +
  scale_linetype_manual(values = pppower_method_linetypes) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "B. Binary Outcomes: Power Curves by Classifier Quality",
    subtitle = "n = 200 labeled, N = 10,000 unlabeled, prevalence = 30%",
    x = "Effect Size (Delta)",
    y = "Power",
    color = "Method",
    linetype = "Method"
  ) +
  theme_pub()

# --- Figure 3: Sample Size Curves (Continuous) ---
cat("  Figure 3: Sample size curves (continuous)...\n")

ss_cont_long <- rbind(
  data.frame(ss_cont_results[, c("delta", "r2_label")],
             n = ss_cont_results$n_ppi, method = "PPI++"),
  data.frame(ss_cont_results[, c("delta", "r2_label")],
             n = ss_cont_results$n_vanilla, method = "Vanilla PPI"),
  data.frame(ss_cont_results[, c("delta", "r2_label")],
             n = ss_cont_results$n_classical, method = "Classical")
)
ss_cont_long <- ss_cont_long[!is.na(ss_cont_long$n), ]

fig3 <- ggplot(ss_cont_long, aes(x = delta, y = n, color = method, linetype = method)) +
  geom_line(linewidth = 1.8) +
  facet_wrap(~r2_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = colors_method) +
  scale_linetype_manual(values = pppower_method_linetypes) +
  labs(
    title = "C. Continuous Outcomes: Required Sample Size for 80% Power",
    subtitle = expression(paste("N = 10,000 unlabeled, ", sigma[Y], " = 1, ", alpha, " = 0.05")),
    x = "Effect Size (Delta)",
    y = "Required Labeled n",
    color = "Method",
    linetype = "Method"
  ) +
  theme_pub()

# --- Figure 4: Sample Size Curves (Binary) ---
cat("  Figure 4: Sample size curves (binary)...\n")

ss_bin_long <- rbind(
  data.frame(ss_bin_results[, c("delta", "clf_label")],
             n = ss_bin_results$n_ppi, method = "PPI++"),
  data.frame(ss_bin_results[, c("delta", "clf_label")],
             n = ss_bin_results$n_vanilla, method = "Vanilla PPI"),
  data.frame(ss_bin_results[, c("delta", "clf_label")],
             n = ss_bin_results$n_classical, method = "Classical")
)
ss_bin_long <- ss_bin_long[!is.na(ss_bin_long$n), ]

fig4 <- ggplot(ss_bin_long, aes(x = delta, y = n, color = method, linetype = method)) +
  geom_line(linewidth = 1.8) +
  facet_wrap(~clf_label, nrow = 1, scales = "free_y") +
  scale_color_manual(values = colors_method) +
  scale_linetype_manual(values = pppower_method_linetypes) +
  labs(
    title = "D. Binary Outcomes: Required Sample Size for 80% Power",
    subtitle = "N = 10,000 unlabeled, prevalence = 30%, alpha = 0.05",
    x = "Effect Size (Delta)",
    y = "Required Labeled n",
    color = "Method",
    linetype = "Method"
  ) +
  theme_pub()

# --- Figure 5: Power by Sample Size (heatmap style) ---
cat("  Figure 5: Power gain summary...\n")

# Compute power gain for each scenario
gain_summary <- data.frame(
  outcome = c(rep("Continuous", 4), rep("Binary", 4)),
  quality = c(names(r2_levels), names(classifier_levels)),
  power_ppi = NA,
  power_classical = NA
)

# Fixed parameters for comparison
n_fixed <- 200
delta_cont_fixed <- 0.3
delta_bin_fixed <- 0.05

for (i in 1:4) {
  r2 <- r2_levels[[i]]
  var_f <- r2 * sigma_y^2
  var_res <- (1 - r2) * sigma_y^2
  cov_yf <- r2 * sigma_y^2

  gain_summary$power_ppi[i] <- power_ppi_mean(
    delta = delta_cont_fixed, N = N_cont, n = n_fixed, alpha = alpha,
    var_f = var_f, var_res = var_res, cov_y_f = cov_yf
  )
  gain_summary$power_classical[i] <- {
    se <- sigma_y / sqrt(n_fixed)
    mu <- abs(delta_cont_fixed) / se
    1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
  }
}

for (i in 1:4) {
  clf <- classifier_levels[[i]]

  gain_summary$power_ppi[i + 4] <- power_ppi_mean(
    delta = delta_bin_fixed, N = N_bin, n = n_fixed, alpha = alpha,
    metrics = list(sensitivity = clf["sens"], specificity = clf["spec"],
                   p_y = prevalence, m_obs = n_fixed),
    metric_type = "classification"
  )
  gain_summary$power_classical[i + 4] <- {
    se <- sqrt(var_y_bin / n_fixed)
    mu <- abs(delta_bin_fixed) / se
    1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
  }
}

gain_summary$gain <- gain_summary$power_ppi - gain_summary$power_classical
gain_summary$quality <- factor(gain_summary$quality,
                                levels = c(names(r2_levels), names(classifier_levels)))

fig5 <- ggplot(gain_summary, aes(x = quality, y = gain, fill = outcome)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  scale_fill_manual(values = pppower_outcome_colors) +
  labs(
    title = "E. Power Gain from PPI++ vs Classical",
    subtitle = sprintf("n = %d, Continuous: Delta = %.1f, Binary: Delta = %.2f",
                       n_fixed, delta_cont_fixed, delta_bin_fixed),
    x = "Prediction Quality",
    y = "Power Gain (PPI++ - Classical)",
    fill = "Outcome Type"
  ) +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# =============================================================================
# PART 4: SAVE FIGURES
# =============================================================================
cat("\nPart 4: Saving Figures\n")
cat("----------------------\n")

# Helper: save PDF + PNG at 400 DPI
save_fig <- function(plot, path, width = 12, height = 4.5) {
  ggsave(paste0(path, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(path, ".png"), plot, width = width, height = height,
         dpi = 400, bg = "white")
}

# Individual figures
save_fig(fig1, "simulation_studies/fig_power_continuous")
save_fig(fig2, "simulation_studies/fig_power_binary")
save_fig(fig3, "simulation_studies/fig_samplesize_continuous")
save_fig(fig4, "simulation_studies/fig_samplesize_binary")
save_fig(fig5, "simulation_studies/fig_power_gain_summary", width = 10, height = 5)

cat("  Saved individual figures (PDF + PNG)\n")

# Combined figure (main figure for paper)
combined <- grid.arrange(
  fig1, fig2, fig3, fig4,
  ncol = 1,
  heights = c(1, 1, 1, 1)
)

save_fig(combined, "simulation_studies/fig_main_power_analysis", width = 12, height = 15)

cat("  Saved combined figure (PDF + PNG)\n")

# =============================================================================
# PART 5: SUMMARY STATISTICS
# =============================================================================
cat("\n=============================================================================\n")
cat("SUMMARY: Sample Size Reduction from PPI++\n")
cat("=============================================================================\n\n")

cat("CONTINUOUS OUTCOMES (Delta = 0.3, 80% power):\n")
cat("---------------------------------------\n")
for (label in names(r2_levels)) {
  r2 <- r2_levels[[label]]
  rho <- sqrt(r2)
  sigma_y2 <- sigma_y^2
  sigma_f2 <- r2 * sigma_y^2
  cov_yf <- rho * sigma_y * sqrt(sigma_f2)

  n_ppi <- tryCatch(
    power_ppi_mean(delta = 0.3, N = N_cont, n = NULL, power = 0.8, alpha = 0.05,
                   sigma_y2 = sigma_y2, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
                   lambda_mode = "oracle"),
    error = function(e) NA
  )
  n_classical <- sigma_y^2 * ((z_alpha + z_beta) / 0.3)^2
  reduction <- (1 - n_ppi / n_classical) * 100

  cat(sprintf("  %s: n_PPI++ = %d vs n_Classical = %d (%.0f%% reduction)\n",
              label, ceiling(n_ppi), ceiling(n_classical), reduction))
}

cat("\nBINARY OUTCOMES (Delta = 0.05, 80% power):\n")
cat("--------------------------------------\n")
for (label in names(classifier_levels)) {
  clf <- classifier_levels[[label]]

  n_ppi <- tryCatch(
    power_ppi_mean(delta = 0.05, N = N_bin, n = NULL, power = 0.8, alpha = 0.05,
                   metrics = list(sensitivity = clf["sens"], specificity = clf["spec"],
                                  p_y = prevalence, m_obs = 200),
                   metric_type = "classification",
                   lambda_mode = "oracle"),
    error = function(e) NA
  )
  n_classical <- var_y_bin * ((z_alpha + z_beta) / 0.05)^2
  reduction <- (1 - n_ppi / n_classical) * 100

  cat(sprintf("  %s: n_PPI++ = %d vs n_Classical = %d (%.0f%% reduction)\n",
              label, ceiling(n_ppi), ceiling(n_classical), reduction))
}

cat("\n=============================================================================\n")
cat("Figures saved to simulation_studies/\n")
cat("=============================================================================\n")
