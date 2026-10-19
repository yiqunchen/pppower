# =============================================================================
# SETTING 8: MEAN ESTIMATION - HEAVY-TAILED (t-DISTRIBUTION)
# =============================================================================
cat("Setting 8: Mean Estimation - Heavy-Tailed (t-distribution, df=5)\n")
cat("-----------------------------------------------------------------\n")
setting8_start <- Sys.time()

# Parameters for t-distribution with df=5
# Var(t_df) = df/(df-2) for df > 2
df_t <- 5
n_values_t <- c(20, 40, 60, 80, 100)
N_t_values <- c(200, 500)
sigma_Y2_t <- df_t / (df_t - 2)  # Variance of t_5
delta_t <- 0.3

rho_values_t <- c(0.3, 0.45, 0.5, 0.7, 0.9)  # Gaussian copula correlation; R1: added low values

# R1 fix: the t5 marginals have a slightly lower correlation than the Gaussian
# copula parameter, so compute the actual correlation once per rho.
rng_state_8 <- .Random.seed
set.seed(8)
z1_8 <- rnorm(2e6); z2_8 <- rnorm(2e6)
rho_t_actual <- sapply(rho_values_t, function(r)
  cor(qt(pnorm(z1_8), df_t), qt(pnorm(r * z1_8 + sqrt(1 - r^2) * z2_8), df_t)))
names(rho_t_actual) <- rho_values_t
rm(z1_8, z2_8)
assign(".Random.seed", rng_state_8, envir = .GlobalEnv)
alpha <- 0.05
R <- 2000  # Monte Carlo replications (R1 revision: was 1000)

# Calculate total iterations for progress
total_configs_8 <- length(rho_values_t) * length(n_values_t) * length(N_t_values)
current_config_8 <- 0

results_tdist <- data.frame()

for (rho in rho_values_t) {
  for (N_t in N_t_values) {
    for (n in n_values_t) {
      current_config_8 <- current_config_8 + 1
      progress_bar(current_config_8, total_configs_8, "Setting 8",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho, N_t, n),
                   start_time = setting8_start)

      # Oracle lambda
      sigma_f2_t <- sigma_Y2_t
      rho_actual <- rho_t_actual[[as.character(rho)]]
      cov_Yf_t <- rho_actual * sqrt(sigma_Y2_t * sigma_f2_t)
      lambda_oracle <- cov_Yf_t / ((1 + n / N_t) * sigma_f2_t)

      # Theoretical power using package functions
      power_theo_ppi <- theo_power_ppi_mean(delta_t, n, N_t, sigma_Y2_t, sigma_f2_t, cov_Yf_t, alpha)
      power_theo_classical <- theo_power_classical_onesample(delta_t, n, sigma_Y2_t, alpha)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        # Generate correlated t-distributed variables using copula approach
        # Step 1: Generate correlated normals
        Sigma_norm <- matrix(c(1, rho, rho, 1), 2, 2)
        L_norm <- chol(Sigma_norm)
        Z_L <- matrix(rnorm(n * 2), n, 2) %*% L_norm

        # Step 2: Transform to uniform via normal CDF
        U_L <- pnorm(Z_L)

        # Step 3: Transform to t via t inverse CDF, then scale
        Y_L <- qt(U_L[, 1], df = df_t) + delta_t
        f_L <- qt(U_L[, 2], df = df_t) + delta_t

        # Same for unlabeled
        Z_U <- matrix(rnorm(N_t * 2), N_t, 2) %*% L_norm
        U_U <- pnorm(Z_U)
        f_U <- qt(U_U[, 2], df = df_t) + delta_t

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

      results_tdist <- rbind(results_tdist, data.frame(
        setting = "Mean (t-dist, df=5)",
        rho = rho,
        rho_actual = rho_actual,
        n = n,
        N = N_t,
        delta = delta_t,
        power_theo_ppi = power_theo_ppi,
        power_emp_ppi = power_emp_ppi,
        power_theo_classical = power_theo_classical,
        power_emp_classical = power_emp_classical,
        mc_reps = R
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting8_start, units = "secs")))
