#!/usr/bin/env Rscript
# =============================================================================
# Power Curves for Binary Outcomes using Sensitivity/Specificity
# =============================================================================
# This script demonstrates the recommended user-friendly interface for
# computing PPI++ (= EIF-optimal) power with binary classification metrics.

library(pppower)

# -----------------------------------------------------------------------------
# Example: AI-assisted disease screening
# -----------------------------------------------------------------------------
# You have an AI model that predicts disease status from images.
# You know: sensitivity = 0.85, specificity = 0.90, prevalence = 0.30
# Question: How much power do you have to detect a 5% difference in prevalence?

cat("=== Power Curves for Binary Classification (PPI++ / EIF) ===\n\n")

# Fixed parameters
sensitivity <- 0.85
specificity <- 0.90
prevalence <- 0.30
N <- 5000           # unlabeled samples (cheap to get AI predictions)
alpha <- 0.05

# Vary: labeled sample size and effect size
n_values <- c(50, 100, 200, 500)
delta_values <- seq(0, 0.15, by = 0.01)

# Compute power grid
results <- expand.grid(n = n_values, delta = delta_values)
results$power_ppiplus <- NA
results$power_classical <- NA

for (i in seq_len(nrow(results))) {
  n <- results$n[i]
  delta <- results$delta[i]

  # PPI++ power (uses AI predictions)
  results$power_ppiplus[i] <- power_ppi_mean(
    delta = delta,
    N = N,
    n = n,
    alpha = alpha,
    metrics = list(
      sensitivity = sensitivity,
      specificity = specificity,
      p_y = prevalence,
      m_obs = n
    ),
    metric_type = "classification"
  )

  # Classical power (no AI predictions) - just uses labeled data
  # Var(Y) for Bernoulli = p(1-p)
  var_y <- prevalence * (1 - prevalence)
  se_classical <- sqrt(var_y / n)
  z_alpha <- qnorm(1 - alpha / 2)
  mu <- abs(delta) / se_classical
  results$power_classical[i] <- 1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
}

# Print table
cat("Labeled n | Delta | PPI++ Power | Classical Power | Gain\n")
cat(strrep("-", 60), "\n")
for (i in seq_len(nrow(results))) {
  if (results$delta[i] %in% c(0, 0.05, 0.10, 0.15)) {
    gain <- results$power_ppiplus[i] - results$power_classical[i]
    cat(sprintf("%9d | %5.2f | %11.3f | %15.3f | %+.3f\n",
                results$n[i], results$delta[i],
                results$power_ppiplus[i], results$power_classical[i], gain))
  }
}

# -----------------------------------------------------------------------------
# Plot power curves
# -----------------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Reshape for plotting
  plot_data <- rbind(
    data.frame(n = results$n, delta = results$delta,
               power = results$power_ppiplus, method = "PPI++ (with AI)"),
    data.frame(n = results$n, delta = results$delta,
               power = results$power_classical, method = "Classical (no AI)")
  )
  plot_data$n_label <- paste0("n = ", plot_data$n)
  plot_data$n_label <- factor(plot_data$n_label,
                               levels = paste0("n = ", sort(unique(plot_data$n))))

  p <- ggplot(plot_data, aes(x = delta, y = power, color = method, linetype = method)) +
    geom_line(linewidth = 1.8) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50", alpha = 0.7) +
    geom_hline(yintercept = alpha, linetype = "dotted", color = "gray50", alpha = 0.7) +
    facet_wrap(~n_label, nrow = 1) +
    scale_color_manual(values = pppower_method_colors) +
    scale_linetype_manual(values = pppower_method_linetypes) +
    labs(
      title = "Power Curves: PPI++ vs Classical Estimation",
      subtitle = sprintf("Sensitivity=%.0f%%, Specificity=%.0f%%, Prevalence=%.0f%%, N=%d unlabeled",
                         sensitivity*100, specificity*100, prevalence*100, N),
      x = "Effect Size (Delta)",
      y = "Power",
      color = "Method",
      linetype = "Method"
    ) +
    theme_pppower() +
    annotate("text", x = 0.14, y = 0.82, label = "80% power",
             size = 3.5, color = "gray40", hjust = 1)

  ggsave("simulation_studies/power_curve_binary.pdf", p,
         width = 12, height = 4.5, device = cairo_pdf)
  ggsave("simulation_studies/power_curve_binary.png", p,
         width = 12, height = 4.5, dpi = 400, bg = "white")
  cat("\nPlot saved to: simulation_studies/power_curve_binary.{pdf,png}\n")
}

# -----------------------------------------------------------------------------
# Efficiency summary
# -----------------------------------------------------------------------------
cat("\n=== Efficiency Summary ===\n")
cat(sprintf("AI classifier: sensitivity=%.0f%%, specificity=%.0f%%\n",
            sensitivity*100, specificity*100))

# Compute effective sample size multiplier
# For binary: error_rate determines variance reduction
error_rate <- (1 - sensitivity) * prevalence + (1 - specificity) * (1 - prevalence)
cat(sprintf("Error rate: %.1f%%\n", error_rate * 100))

# The variance reduction from PPI++ depends on the correlation between Y and f
# For classification: Cov(Y,f) = TP_rate - p_y * p_hat
tp_rate <- sensitivity * prevalence
fp_rate <- (1 - specificity) * (1 - prevalence)
p_hat <- tp_rate + fp_rate
cov_yf <- tp_rate - prevalence * p_hat
var_y <- prevalence * (1 - prevalence)
var_f <- p_hat * (1 - p_hat)
rho <- cov_yf / sqrt(var_y * var_f)

cat(sprintf("Correlation(Y, f): %.3f\n", rho))
cat(sprintf("Variance reduction factor: %.1f%%\n", (1 - rho^2) * 100))
cat(sprintf("Effective sample size multiplier: %.1fx\n", 1 / (1 - rho^2)))

cat("\n=== Interpretation ===\n")
cat("With this AI classifier, PPI++ reduces the required labeled sample size\n")
cat(sprintf("by approximately %.0f%% compared to classical estimation.\n", rho^2 * 100))
