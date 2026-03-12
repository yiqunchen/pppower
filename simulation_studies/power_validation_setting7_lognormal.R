# =============================================================================
# SETTING 7: MEAN ESTIMATION - LOG-NORMAL OUTCOME (RIGHT-SKEWED)
# =============================================================================
cat("Setting 7: Mean Estimation - Log-Normal Outcome (Skewed)\n")
cat("---------------------------------------------------------\n")
setting7_start <- Sys.time()

# Parameters for log-normal: if X ~ N(mu, sigma^2), then exp(X) ~ LogNormal
# We want E[Y] = exp(mu + sigma^2/2), Var(Y) = (exp(sigma^2) - 1) * exp(2*mu + sigma^2)
n_values_lognorm <- c(20, 40, 60, 80, 100)
N_lognorm_values <- c(200, 500)
mu_log <- 0  # log-scale mean
sigma_log <- 0.5  # log-scale sd

# True mean and variance of log-normal
true_mean_lognorm <- exp(mu_log + sigma_log^2 / 2)
sigma_Y2_lognorm <- (exp(sigma_log^2) - 1) * exp(2 * mu_log + sigma_log^2)

# Effect size (shift in log-scale mean)
delta_log <- 0.15  # shift in log-scale -> multiplicative effect on original scale
delta_lognorm <- exp(mu_log + delta_log + sigma_log^2 / 2) - true_mean_lognorm

rho_values_lognorm <- c(0.5, 0.7, 0.9)
alpha <- 0.05
R <- 1000  # Monte Carlo replications

# Calculate total iterations for progress
total_configs_7 <- length(rho_values_lognorm) * length(n_values_lognorm) * length(N_lognorm_values)
current_config_7 <- 0

results_lognorm <- data.frame()

for (rho in rho_values_lognorm) {
  # For log-normal, we generate correlated normals then exponentiate
  # Correlation in log-scale translates approximately to correlation in original scale

  for (N_lognorm in N_lognorm_values) {
    for (n in n_values_lognorm) {
      current_config_7 <- current_config_7 + 1
      progress_bar(current_config_7, total_configs_7, "Setting 7",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho, N_lognorm, n),
                   start_time = setting7_start)

      # Oracle lambda
      sigma_f2_lognorm <- sigma_Y2_lognorm  # assume same variance for predictions
      cov_Yf_lognorm <- rho * sqrt(sigma_Y2_lognorm * sigma_f2_lognorm)
      lambda_oracle <- cov_Yf_lognorm / ((1 + n / N_lognorm) * sigma_f2_lognorm)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_mean(delta_lognorm, n, N_lognorm,
                                            sigma_Y2_lognorm, sigma_f2_lognorm, cov_Yf_lognorm, alpha)
      power_theo_classical <- theo_power_classical_onesample(delta_lognorm, n, sigma_Y2_lognorm, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      # Correlation matrix for log-scale
      Sigma_log <- matrix(c(1, rho, rho, 1), 2, 2) * sigma_log^2
      L_log <- chol(Sigma_log)

      for (r in 1:R) {
        # Generate correlated log-normals under H1
        Z_L <- matrix(rnorm(n * 2), n, 2) %*% L_log
        Y_L <- exp(Z_L[, 1] + mu_log + delta_log)  # Shifted mean
        f_L <- exp(Z_L[, 2] + mu_log + delta_log)

        Z_U <- matrix(rnorm(N_lognorm * 2), N_lognorm, 2) %*% L_log
        f_U <- exp(Z_U[, 2] + mu_log + delta_log)

        # PPI++ test with oracle lambda
        test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = true_mean_lognorm, alpha = alpha,
                                      lambda_oracle = lambda_oracle)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical test
        test_classical <- classical_mean_test(Y_L, theta0 = true_mean_lognorm, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R

      results_lognorm <- rbind(results_lognorm, data.frame(
        setting = "Mean (Log-Normal)",
        rho = rho,
        n = n,
        N = N_lognorm,
        delta = delta_lognorm,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting7_start, units = "secs")))
