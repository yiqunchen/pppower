#!/usr/bin/env Rscript
# =============================================================================
# Setting 17: N/n Ratio Sensitivity
# =============================================================================
# Current sims always have N >> n (ratio 10-200x).
# This setting varies the N/n ratio from 1 (N = n, no unlabeled advantage)
# to large (N >> n), showing the PPI advantage vanishing as N/n -> 1.
#
# Fixed n = 100, varying N = n, 2n, 5n, 10n, 50n, 100n.

cat("Setting 17: N/n Ratio Sensitivity\n")
cat("-----------------------------------\n")
setting17_start <- Sys.time()

library(pppower)

# Parameters
n_17 <- 50
Nn_ratios <- c(1, 1.5, 2, 3, 5, 10, 20, 50, 100)
N_values_17 <- n_17 * Nn_ratios
rho_values_17 <- c(0.5, 0.7, 0.9)
delta_17 <- 0.2
sigma_Y2_17 <- 1
sigma_f2_17 <- 1
alpha_17 <- 0.05
R_17 <- 1000

grid_17 <- expand.grid(
  rho = rho_values_17,
  N = N_values_17,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_17 <- nrow(grid_17)

results_Nn_ratio <- vector("list", total_17)

for (i in seq_len(total_17)) {
  row <- grid_17[i, ]
  Nn_ratio <- row$N / n_17
  progress_bar(i, total_17, "Setting 17",
               sprintf("rho=%.1f, N/n=%.1f: ", row$rho, Nn_ratio),
               start_time = setting17_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_17 * sigma_f2_17)
  Sigma <- matrix(c(sigma_Y2_17, cov_Yf, cov_Yf, sigma_f2_17), 2, 2)
  L <- chol(Sigma)

  lambda_oracle <- cov_Yf / ((1 + n_17 / row$N) * sigma_f2_17)

  power_theo_ppi <- theo_power_ppi_mean(delta_17, n_17, row$N,
                                         sigma_Y2_17, sigma_f2_17, cov_Yf, alpha_17)
  power_theo_classical <- theo_power_classical_onesample(delta_17, n_17,
                                                          sigma_Y2_17, alpha_17)

  rej_ppi <- 0
  rej_classical <- 0

  for (r in seq_len(R_17)) {
    Z_L <- matrix(rnorm(n_17 * 2), n_17, 2) %*% L
    Y_L <- Z_L[, 1] + delta_17
    f_L <- Z_L[, 2] + delta_17

    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    f_U <- Z_U[, 2] + delta_17

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha_17,
                                  lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha_17)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_Nn_ratio[[i]] <- data.frame(
    setting = "N/n Ratio",
    rho = row$rho,
    n = n_17,
    N = row$N,
    Nn_ratio = Nn_ratio,
    delta = delta_17,
    power_theo_ppi = power_theo_ppi,
    power_emp_ppi = rej_ppi / R_17,
    power_theo_classical = power_theo_classical,
    power_emp_classical = rej_classical / R_17,
    stringsAsFactors = FALSE
  )
}

results_Nn_ratio <- do.call(rbind, results_Nn_ratio)

cat(sprintf("Setting 17 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting17_start, units = "secs")))
cat("  Max |theo - emp| PPI++:",
    max(abs(results_Nn_ratio$power_theo_ppi - results_Nn_ratio$power_emp_ppi)), "\n")
# Show convergence: PPI -> classical as N/n -> 1
r1 <- results_Nn_ratio[results_Nn_ratio$Nn_ratio == 1, ]
cat(sprintf("  At N/n=1:  PPI++ power ≈ classical (max diff = %.4f)\n\n",
            max(abs(r1$power_emp_ppi - r1$power_emp_classical))))
