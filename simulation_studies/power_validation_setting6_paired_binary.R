# =============================================================================
# SETTING 6: PAIRED PROPORTION TEST - BINARY OUTCOME (McNemar-style)
# =============================================================================
cat("Setting 6: Paired Proportion Test - Binary Outcome\n")
cat("---------------------------------------------------\n")
setting6_start <- Sys.time()

# For paired binary data: D = Y_A - Y_B can be -1, 0, or 1
# Var(D) depends on marginal and joint probabilities
# Simpler approach: treat as difference in matched proportions

n_values_paired_bin <- c(20, 40, 60, 80, 100)
N_paired_bin_values <- c(200, 500)
p0_paired <- 0.3  # Base probability for each condition
delta_paired_bin <- 0.08  # Difference in paired proportions
classifier_settings <- list(
  c(sens = 0.70, spec = 0.70),
  c(sens = 0.85, spec = 0.85),
  c(sens = 0.95, spec = 0.95)
)
alpha <- 0.05
R <- 1000  # Monte Carlo replications

# Calculate total iterations for progress
total_configs_6 <- length(classifier_settings) * length(n_values_paired_bin) * length(N_paired_bin_values)
current_config_6 <- 0

results_paired_bin <- data.frame()

for (classifier in classifier_settings) {
  sens <- classifier["sens"]
  spec <- classifier["spec"]

  # Under H1: p_A = p0 + delta/2, p_B = p0 - delta/2
  p_A <- p0_paired + delta_paired_bin / 2
  p_B <- p0_paired - delta_paired_bin / 2

  # Paired Bernoulli with within-subject Pearson correlation
  rho_within <- 0.3  # Within-subject correlation
  p_11 <- p_A * p_B + rho_within * sqrt(p_A * (1 - p_A) * p_B * (1 - p_B))
  p_11 <- min(max(p_11, max(0, p_A + p_B - 1)), min(p_A, p_B))
  p_10 <- p_A - p_11
  p_01 <- p_B - p_11
  p_00 <- 1 - p_11 - p_10 - p_01
  p_vec <- c(p_11, p_10, p_01, p_00)
  p_vec <- pmax(p_vec, 0)
  p_vec <- p_vec / sum(p_vec)
  p_11 <- p_vec[1]; p_10 <- p_vec[2]; p_01 <- p_vec[3]; p_00 <- p_vec[4]

  # Variance of D = Y_A - Y_B
  # E[D] = p_A - p_B = delta
  # Var(D) = Var(Y_A) + Var(Y_B) - 2*Cov(Y_A, Y_B)
  cov_AB <- p_11 - p_A * p_B
  sigma_D2_bin <- p_A * (1 - p_A) + p_B * (1 - p_B) - 2 * cov_AB

  for (N_paired_bin in N_paired_bin_values) {
    for (n in n_values_paired_bin) {
      current_config_6 <- current_config_6 + 1
      progress_bar(current_config_6, total_configs_6, "Setting 6",
                   sprintf("sens=%.0f%%, N=%d, n=%d: ", sens * 100, N_paired_bin, n),
                   start_time = setting6_start)

      # Theoretical power using package functions (binary paired CLT)
      power_theo_ppi <- theo_power_ppi_paired_binary(
        delta_paired_bin, n, N_paired_bin,
        p_A, p_B, rho_within, sens, spec, alpha
      )
      power_theo_classical <- theo_power_classical_onesample(delta_paired_bin, n, sigma_D2_bin, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Generate paired binary (Y_A, Y_B) via joint Bernoulli
        pairs <- rmultinom(1, n, prob = p_vec)
        Y_A <- c(rep(1, pairs[1] + pairs[2]), rep(0, pairs[3] + pairs[4]))
        Y_B <- c(rep(1, pairs[1]), rep(0, pairs[2]), rep(1, pairs[3]), rep(0, pairs[4]))
        perm <- sample.int(n)
        Y_A <- Y_A[perm]
        Y_B <- Y_B[perm]
        D_L <- Y_A - Y_B

        # Classifier predictions
        f_A <- ifelse(Y_A == 1, rbinom(n, 1, sens), rbinom(n, 1, 1 - spec))
        f_B <- ifelse(Y_B == 1, rbinom(n, 1, sens), rbinom(n, 1, 1 - spec))
        fD_L <- f_A - f_B

        # Unlabeled (joint Bernoulli)
        pairs_U <- rmultinom(1, N_paired_bin, prob = p_vec)
        Y_A_U <- c(rep(1, pairs_U[1] + pairs_U[2]), rep(0, pairs_U[3] + pairs_U[4]))
        Y_B_U <- c(rep(1, pairs_U[1]), rep(0, pairs_U[2]), rep(1, pairs_U[3]), rep(0, pairs_U[4]))
        perm_U <- sample.int(N_paired_bin)
        Y_A_U <- Y_A_U[perm_U]
        Y_B_U <- Y_B_U[perm_U]
        f_A_U <- ifelse(Y_A_U == 1, rbinom(N_paired_bin, 1, sens), rbinom(N_paired_bin, 1, 1 - spec))
        f_B_U <- ifelse(Y_B_U == 1, rbinom(N_paired_bin, 1, sens), rbinom(N_paired_bin, 1, 1 - spec))
        fD_U <- f_A_U - f_B_U

        # PPI++ test (use plug-in lambda for paired binary since oracle is complex)
        test_ppi <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0, alpha = alpha)
        rejections_ppi <- rejections_ppi + test_ppi$reject

        # Classical test
        test_classical <- classical_mean_test(D_L, theta0 = 0, alpha = alpha)
        rejections_classical <- rejections_classical + test_classical$reject
      }

      power_emp_ppi <- rejections_ppi / R
      power_emp_classical <- rejections_classical / R

      results_paired_bin <- rbind(results_paired_bin, data.frame(
        setting = "Paired Proportion Test (Binary)",
        sens = sens,
        spec = spec,
        n = n,
        N = N_paired_bin,
        delta = delta_paired_bin,
        rho_within = rho_within,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting6_start, units = "secs")))
