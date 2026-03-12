# =============================================================================
# SETTING 3: TWO-SAMPLE T-TEST - CONTINUOUS OUTCOME
# =============================================================================
cat("Setting 3: Two-Sample t-Test - Continuous Outcome\n")
cat("--------------------------------------------------\n")
setting3_start <- Sys.time()

# Parameters
n_values_ttest <- c(20, 40, 60, 80, 100)  # Per group
N_ttest_values <- c(200, 500)  # Per group
delta_ttest <- 0.3
sigma_Y2 <- 1
rho_values <- c(0.5, 0.7, 0.9)
alpha <- 0.05
R <- 1000  # Monte Carlo replications

# Calculate total iterations for progress
total_configs_3 <- length(rho_values) * length(n_values_ttest) * length(N_ttest_values)
current_config_3 <- 0

results_cont_ttest <- data.frame()

for (rho in rho_values) {
  sigma_f2 <- 1
  cov_Yf <- rho * sqrt(sigma_Y2 * sigma_f2)
  Sigma <- matrix(c(sigma_Y2, cov_Yf, cov_Yf, sigma_f2), 2, 2)
  L <- chol(Sigma)

  for (N_ttest in N_ttest_values) {
    for (n in n_values_ttest) {
      current_config_3 <- current_config_3 + 1
      progress_bar(current_config_3, total_configs_3, "Setting 3",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho, N_ttest, n),
                   start_time = setting3_start)

      # Oracle lambda (EIF-optimal, same for both groups due to same rho)
      lambda_oracle <- cov_Yf / ((1 + n / N_ttest) * sigma_f2)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_twosample(
        delta_ttest, n, n, N_ttest, N_ttest,
        sigma_Y2, sigma_f2, cov_Yf,
        sigma_Y2, sigma_f2, cov_Yf, alpha
      )
      power_theo_classical <- theo_power_classical_twosample(delta_ttest, n, sigma_Y2, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Group A (treatment): mean = delta/2
        Z_A_L <- matrix(rnorm(n * 2), n, 2) %*% L
        Y_A <- Z_A_L[, 1] + delta_ttest / 2
        f_A_L <- Z_A_L[, 2] + delta_ttest / 2
        Z_A_U <- matrix(rnorm(N_ttest * 2), N_ttest, 2) %*% L
        f_A_U <- Z_A_U[, 2] + delta_ttest / 2

        # Group B (control): mean = -delta/2
        Z_B_L <- matrix(rnorm(n * 2), n, 2) %*% L
        Y_B <- Z_B_L[, 1] - delta_ttest / 2
        f_B_L <- Z_B_L[, 2] - delta_ttest / 2
        Z_B_U <- matrix(rnorm(N_ttest * 2), N_ttest, 2) %*% L
        f_B_U <- Z_B_U[, 2] - delta_ttest / 2

        # PPI++ test with oracle lambda
        test_ppi <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                                  delta0 = 0, alpha = alpha,
                                  lambda_A_oracle = lambda_oracle,
                                  lambda_B_oracle = lambda_oracle)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical test
        test_classical <- classical_ttest(Y_A, Y_B, delta0 = 0, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R

      results_cont_ttest <- rbind(results_cont_ttest, data.frame(
        setting = "t-Test (Continuous)",
        rho = rho,
        n = n,
        N = N_ttest,
        delta = delta_ttest,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting3_start, units = "secs")))
