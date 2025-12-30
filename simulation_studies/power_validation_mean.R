#' Power Validation for PPI Mean Estimation
#'
#' This script validates that the analytical power formulas in pppower
#' match empirical power from Monte Carlo simulation.
#'
#' Key question: Does the closed-form power formula accurately predict
#' the actual rejection rate when we run the test many times?

source("simulation_studies/sim_utils.R")

# =============================================================================
# Configuration
# =============================================================================

# Sample sizes
n_values <- c(50, 100, 200)   # labeled
N <- 1000                      # unlabeled

# Effect sizes (under H1)
delta_values <- seq(0, 0.5, by = 0.1)

# DGP parameters
sigma_y <- 1      # SD of outcome
sigma_f <- 0.8    # SD of predictor (determines predictor quality)
rho <- 0.7        # correlation between Y and f

# MC settings
R <- 500          # replicates (increase for more precision)
alpha <- 0.05     # significance level

# =============================================================================
# Run Validation: PPI
# =============================================================================

cat("\n=== Validating Vanilla PPI Power Formula ===\n\n")

results_ppi <- power_curve_grid(
  n_values = n_values,
  delta_values = delta_values,
  N = N,
  sigma_y = sigma_y,
  sigma_f = sigma_f,
  rho = rho,
  R = R,
  use_ppi_plus = FALSE,
  seed = 42
)

cat("Vanilla PPI Results:\n")
cat(sprintf("%-5s %-8s %-12s %-12s %-10s\n",
            "n", "delta", "Emp Power", "Ana Power", "Diff"))
cat(rep("-", 55), "\n", sep = "")

for (i in 1:nrow(results_ppi)) {
  row <- results_ppi[i, ]
  diff <- row$empirical_power - row$analytical_power
  cat(sprintf("%-5d %-8.2f %-12.3f %-12.3f %-10.3f\n",
              row$n, row$delta, row$empirical_power, row$analytical_power, diff))
}

# =============================================================================
# Run Validation: PPI++
# =============================================================================

cat("\n\n=== Validating PPI++ Power Formula ===\n\n")

results_ppiplus <- power_curve_grid(
  n_values = n_values,
  delta_values = delta_values,
  N = N,
  sigma_y = sigma_y,
  sigma_f = sigma_f,
  rho = rho,
  R = R,
  use_ppi_plus = TRUE,
  seed = 42
)

cat("PPI++ Results:\n")
cat(sprintf("%-5s %-8s %-12s %-12s %-10s\n",
            "n", "delta", "Emp Power", "Ana Power", "Diff"))
cat(rep("-", 55), "\n", sep = "")

for (i in 1:nrow(results_ppiplus)) {
  row <- results_ppiplus[i, ]
  diff <- row$empirical_power - row$analytical_power
  cat(sprintf("%-5d %-8.2f %-12.3f %-12.3f %-10.3f\n",
              row$n, row$delta, row$empirical_power, row$analytical_power, diff))
}

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n\n=== Summary ===\n\n")

# Mean absolute difference
mad_ppi <- mean(abs(results_ppi$empirical_power - results_ppi$analytical_power))
mad_ppiplus <- mean(abs(results_ppiplus$empirical_power - results_ppiplus$analytical_power))

cat(sprintf("Mean Absolute Difference (PPI):   %.4f\n", mad_ppi))
cat(sprintf("Mean Absolute Difference (PPI++): %.4f\n", mad_ppiplus))

# Interpretation
cat("\n")
if (mad_ppi < 0.05 && mad_ppiplus < 0.05) {
  cat("SUCCESS: Analytical formulas closely match empirical power.\n")
} else {
  cat("WARNING: Larger than expected differences. Check assumptions.\n")
}

# =============================================================================
# Optional: Simple Plot (if ggplot2 available)
# =============================================================================

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  # Combine results
  results_ppi$method <- "PPI"
  results_ppiplus$method <- "PPI++"
  all_results <- rbind(results_ppi, results_ppiplus)

  p <- ggplot(all_results, aes(x = analytical_power, y = empirical_power, color = method)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = 2, alpha = 0.7) +
    facet_wrap(~ n, labeller = labeller(n = function(x) paste("n =", x))) +
    labs(
      title = "Analytical vs Empirical Power",
      subtitle = sprintf("N = %d, rho = %.1f, R = %d MC replicates", N, rho, R),
      x = "Analytical Power",
      y = "Empirical Power",
      color = "Method"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

  print(p)

  # Save plot
  ggsave("simulation_studies/power_validation_plot.png", p, width = 8, height = 5, dpi = 150)
  cat("\nPlot saved to: simulation_studies/power_validation_plot.png\n")
}
