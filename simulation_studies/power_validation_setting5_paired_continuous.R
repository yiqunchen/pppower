# =============================================================================
# SETTING 5: PAIRED T-TEST - CONTINUOUS OUTCOME
# =============================================================================
cat("Setting 5: Paired t-Test - Continuous Outcome\n")
cat("----------------------------------------------\n")
setting5_start <- Sys.time()

# Parameters
n_values_paired <- c(20, 40, 60, 80, 100)
N_paired_values <- c(200, 500)
delta_paired <- 0.3
sigma_D2 <- 1  # Variance of differences
rho_D_values <- c(0.5, 0.7, 0.9)
alpha <- 0.05
R <- 1000  # Monte Carlo replications

# Calculate total iterations for progress
total_configs_5 <- length(rho_D_values) * length(n_values_paired) * length(N_paired_values)
current_config_5 <- 0

results_paired <- data.frame()

for (rho_D in rho_D_values) {
  # Construct covariance matrix for (D, f^D)
  sigma_fD2 <- 1
  cov_D_fD <- rho_D * sqrt(sigma_D2 * sigma_fD2)
  Sigma_paired <- matrix(c(sigma_D2, cov_D_fD, cov_D_fD, sigma_fD2), 2, 2)
  L_paired <- chol(Sigma_paired)

  for (N_paired in N_paired_values) {
    for (n in n_values_paired) {
      current_config_5 <- current_config_5 + 1
      progress_bar(current_config_5, total_configs_5, "Setting 5",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho_D, N_paired, n),
                   start_time = setting5_start)

      # Oracle lambda (EIF-optimal based on true population parameters)
      lambda_oracle <- cov_D_fD / ((1 + n / N_paired) * sigma_fD2)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_paired(delta_paired, n, N_paired, sigma_D2, rho_D, alpha)
      power_theo_classical <- theo_power_classical_onesample(delta_paired, n, sigma_D2, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Generate paired (D, f^D) under H1: E[D] = delta
        Z_L <- matrix(rnorm(n * 2), n, 2) %*% L_paired
        D_L <- Z_L[, 1] + delta_paired  # Differences
        fD_L <- Z_L[, 2] + delta_paired  # Predicted differences

        Z_U <- matrix(rnorm(N_paired * 2), N_paired, 2) %*% L_paired
        fD_U <- Z_U[, 2] + delta_paired

        # PPI++ test with oracle lambda
        test_ppi <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0, alpha = alpha,
                                      lambda_oracle = lambda_oracle)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical paired t-test
        test_classical <- classical_mean_test(D_L, theta0 = 0, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R

      results_paired <- rbind(results_paired, data.frame(
        setting = "Paired t-Test (Continuous)",
        rho_D = rho_D,
        n = n,
        N = N_paired,
        delta = delta_paired,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting5_start, units = "secs")))
