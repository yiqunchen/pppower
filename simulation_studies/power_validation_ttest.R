#' Power Validation for PPI Two-Sample t-Test (Difference in Means)
#'
#' This script validates the analytical power formula for testing
#' H0: theta_A - theta_B = 0 under PPI.

source("simulation_studies/sim_utils.R")

#' Single MC replicate for two-sample PPI test
#'
#' @param n_A, n_B Labeled sample sizes per group
#' @param N_A, N_B Unlabeled sample sizes per group
#' @param delta True difference in means (theta_A - theta_B)
#' @param sigma_y SD of Y (same for both groups)
#' @param sigma_f SD of f (same for both groups)
#' @param rho Correlation between Y and f (same for both groups)
#' @param alpha Significance level
#' @return List with delta_hat, se, z_stat, reject
ppi_ttest_one_rep <- function(n_A, n_B, N_A, N_B, delta, sigma_y = 1, sigma_f = 1,
                               rho = 0.5, alpha = 0.05) {
  # Group A: mean = delta/2
  dat_L_A <- generate_bivariate_normal(n_A, mu_y = delta/2, sigma_y, sigma_f, rho)
  dat_U_A <- generate_bivariate_normal(N_A, mu_y = delta/2, sigma_y, sigma_f, rho)

  # Group B: mean = -delta/2
  dat_L_B <- generate_bivariate_normal(n_B, mu_y = -delta/2, sigma_y, sigma_f, rho)
  dat_U_B <- generate_bivariate_normal(N_B, mu_y = -delta/2, sigma_y, sigma_f, rho)

  # PPI estimators for each group
  theta_A_hat <- mean(dat_U_A$f) + mean(dat_L_A$y - dat_L_A$f)
  theta_B_hat <- mean(dat_U_B$f) + mean(dat_L_B$y - dat_L_B$f)

  # Difference estimator
  delta_hat <- theta_A_hat - theta_B_hat

  # Variance components
  var_eps <- sigma_y^2 - 2 * rho * sigma_y * sigma_f + sigma_f^2  # Var(Y - f)

  var_A <- var(dat_U_A$f) / N_A + var(dat_L_A$y - dat_L_A$f) / n_A
  var_B <- var(dat_U_B$f) / N_B + var(dat_L_B$y - dat_L_B$f) / n_B

  se <- sqrt(var_A + var_B)

  # Test
  z_stat <- delta_hat / se
  z_crit <- qnorm(1 - alpha / 2)
  reject <- abs(z_stat) > z_crit

  list(
    delta_hat = delta_hat,
    se = se,
    z_stat = z_stat,
    reject = reject
  )
}

#' Compute empirical power for two-sample test
compute_ttest_power <- function(R = 500, n_A, n_B, N_A, N_B, delta,
                                 sigma_y = 1, sigma_f = 1, rho = 0.5,
                                 alpha = 0.05, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  rejects <- replicate(R, {
    res <- ppi_ttest_one_rep(n_A, n_B, N_A, N_B, delta, sigma_y, sigma_f, rho, alpha)
    res$reject
  })

  emp_power <- mean(rejects)
  mc_se <- sqrt(emp_power * (1 - emp_power) / R)

  # Analytical power
  var_eps <- sigma_y^2 - 2 * rho * sigma_y * sigma_f + sigma_f^2

  # Variance of delta_hat under PPI
  var_delta <- (sigma_f^2 / N_A + var_eps / n_A) + (sigma_f^2 / N_B + var_eps / n_B)
  se_analytical <- sqrt(var_delta)

  z_alpha <- qnorm(1 - alpha / 2)
  analytical_power <- 1 - pnorm(z_alpha - delta / se_analytical) +
                      pnorm(-z_alpha - delta / se_analytical)

  list(
    empirical_power = emp_power,
    mc_se = mc_se,
    analytical_power = analytical_power,
    delta = delta,
    n_A = n_A,
    n_B = n_B
  )
}

# =============================================================================
# Run Validation
# =============================================================================

cat("\n=== Validating PPI Two-Sample t-Test Power ===\n\n")

# Configuration
n_values <- c(50, 100, 150)  # labeled per group
N <- 1000                     # unlabeled per group
delta_values <- seq(0, 0.6, by = 0.15)
sigma_y <- 1
sigma_f <- 0.8
rho <- 0.6
R <- 500
alpha <- 0.05

# Run validation
results <- list()
k <- 1

for (n in n_values) {
  for (delta in delta_values) {
    res <- compute_ttest_power(
      R = R, n_A = n, n_B = n, N_A = N, N_B = N, delta = delta,
      sigma_y = sigma_y, sigma_f = sigma_f, rho = rho, alpha = alpha,
      seed = 42 + k
    )
    results[[k]] <- as.data.frame(res)
    k <- k + 1
  }
}

results_df <- do.call(rbind, results)

# Print results
cat(sprintf("%-5s %-8s %-12s %-12s %-10s\n",
            "n", "delta", "Emp Power", "Ana Power", "Diff"))
cat(rep("-", 55), "\n", sep = "")

for (i in 1:nrow(results_df)) {
  row <- results_df[i, ]
  diff <- row$empirical_power - row$analytical_power
  cat(sprintf("%-5d %-8.2f %-12.3f %-12.3f %-10.3f\n",
              row$n_A, row$delta, row$empirical_power, row$analytical_power, diff))
}

# Summary
mad <- mean(abs(results_df$empirical_power - results_df$analytical_power))
cat(sprintf("\nMean Absolute Difference: %.4f\n", mad))

if (mad < 0.05) {
  cat("SUCCESS: t-test power formula validated.\n")
} else {
  cat("WARNING: Check formula or increase R.\n")
}
