# =============================================================================
# SETTING 1: MEAN ESTIMATION - CONTINUOUS OUTCOME
# =============================================================================
cat("Setting 1: Mean Estimation - Continuous Outcome\n")
cat("------------------------------------------------\n")
setting1_start <- Sys.time()

# Parameters
n_values <- c(20, 40, 60, 80, 100)
N_values <- c(200, 500)
rho_values <- c(0.5, 0.7, 0.9)
delta <- 0.2
sigma_Y2 <- 1
alpha <- 0.05
R <- 1000  # Monte Carlo replications

settings_cont_mean <- list(
  n_values = n_values,
  N_values = N_values,
  rho_values = rho_values,
  sigma_Y2 = sigma_Y2,
  sigma_f2 = 1,
  alpha = alpha,
  R = R
)

# Calculate total iterations for progress
total_configs_1 <- length(rho_values) * length(n_values) * length(N_values)
current_config_1 <- 0

results_cont_mean <- data.frame()

for (rho in rho_values) {
  # Correlation structure for bivariate normal
  sigma_f2 <- 1
  cov_Yf <- rho * sqrt(sigma_Y2 * sigma_f2)
  Sigma <- matrix(c(sigma_Y2, cov_Yf, cov_Yf, sigma_f2), 2, 2)
  L <- chol(Sigma)

  for (N in N_values) {
    for (n in n_values) {
      current_config_1 <- current_config_1 + 1
      progress_bar(current_config_1, total_configs_1, "Setting 1",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho, N, n),
                   start_time = setting1_start)

      # Oracle lambda (EIF-optimal based on true population parameters)
      lambda_oracle <- cov_Yf / ((1 + n / N) * sigma_f2)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_mean(delta, n, N, sigma_Y2, sigma_f2, cov_Yf, alpha)
      power_theo_classical <- theo_power_classical_onesample(delta, n, sigma_Y2, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Generate bivariate normal (Y, f) under H1
        Z_L <- matrix(rnorm(n * 2), n, 2) %*% L
        Y_L <- Z_L[, 1] + delta  # Shift mean by delta
        f_L <- Z_L[, 2] + delta

        Z_U <- matrix(rnorm(N * 2), N, 2) %*% L
        f_U <- Z_U[, 2] + delta

        # PPI++ test with oracle lambda
        test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha,
                                      lambda_oracle = lambda_oracle)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical test
        test_classical <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R

      results_cont_mean <- rbind(results_cont_mean, data.frame(
        setting = "Mean (Continuous)",
        rho = rho,
        n = n,
        N = N,
        delta = delta,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting1_start, units = "secs")))
