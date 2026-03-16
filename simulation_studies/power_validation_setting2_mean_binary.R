# =============================================================================
# SETTING 2: MEAN ESTIMATION - BINARY OUTCOME
# =============================================================================
cat("Setting 2: Mean Estimation - Binary Outcome\n")
cat("--------------------------------------------\n")
setting2_start <- Sys.time()

# Parameters
n_values_bin <- c(20, 40, 60, 80, 100)
N_bin_values <- c(200, 500)
p0 <- 0.3  # Base prevalence
delta_bin <- 0.05  # Effect size
classifier_settings <- list(
  c(sens = 0.70, spec = 0.70),
  c(sens = 0.85, spec = 0.85),
  c(sens = 0.95, spec = 0.95)
)
alpha <- 0.05
R <- 5000  # Monte Carlo replications

# Calculate total iterations for progress
total_configs_2 <- length(classifier_settings) * length(n_values_bin) * length(N_bin_values)
current_config_2 <- 0

results_bin_mean <- data.frame()

for (classifier in classifier_settings) {
  sens <- classifier["sens"]
  spec <- classifier["spec"]

  p <- p0 + delta_bin  # True prevalence under H1

  # Derived quantities
  p_f <- sens * p + (1 - spec) * (1 - p)  # P(f = 1)
  sigma_Y2 <- p * (1 - p)
  sigma_f2 <- p_f * (1 - p_f)
  cov_Yf <- sens * p - p * p_f

  # Correlation
  rho <- cov_Yf / sqrt(sigma_Y2 * sigma_f2)

  for (N_bin in N_bin_values) {
    for (n in n_values_bin) {
      current_config_2 <- current_config_2 + 1
      progress_bar(current_config_2, total_configs_2, "Setting 2",
                   sprintf("sens=%.0f%%, N=%d, n=%d: ", sens * 100, N_bin, n),
                   start_time = setting2_start)

      # Oracle lambda (EIF-optimal based on true population parameters)
      lambda_oracle <- cov_Yf / ((1 + n / N_bin) * sigma_f2)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_mean(delta_bin, n, N_bin, sigma_Y2, sigma_f2, cov_Yf, alpha)
      power_theo_classical <- theo_power_classical_onesample(delta_bin, n, sigma_Y2, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Generate binary outcome
        Y_L <- rbinom(n, 1, p)

        # Generate classifier predictions conditional on Y
        f_L <- ifelse(Y_L == 1,
                      rbinom(n, 1, sens),      # TP or FN
                      rbinom(n, 1, 1 - spec))  # FP or TN

        # Unlabeled
        Y_U_latent <- rbinom(N_bin, 1, p)
        f_U <- ifelse(Y_U_latent == 1,
                      rbinom(N_bin, 1, sens),
                      rbinom(N_bin, 1, 1 - spec))

        # PPI++ test with oracle lambda
        test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = p0, alpha = alpha,
                                      lambda_oracle = lambda_oracle)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical test
        test_classical <- classical_mean_test(Y_L, theta0 = p0, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R
      mc_se_emp_ppi <- sqrt(power_emp_ppi * (1 - power_emp_ppi) / R)
      mc_se_emp_classical <- sqrt(power_emp_classical * (1 - power_emp_classical) / R)

      results_bin_mean <- rbind(results_bin_mean, data.frame(
        setting = "Mean (Binary)",
        sens = sens,
        spec = spec,
        n = n,
        N = N_bin,
        delta = delta_bin,
        rho = rho,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        mc_se_emp_ppi = mc_se_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical,
        mc_se_emp_classical = mc_se_emp_classical,
        mc_reps = R
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting2_start, units = "secs")))
