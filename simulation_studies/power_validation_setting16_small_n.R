#!/usr/bin/env Rscript
# =============================================================================
# Setting 16: Small n Regime (CLT Breakdown)
# =============================================================================
# Pushes to very small labeled samples (n = 15, 20, 25, 30, 50) to test
# where the normal approximation underlying the power formula breaks down.
# Same DGP as Setting 1 (Gaussian mean) — if the formula fails here, it's
# not the DGP's fault but the asymptotic approximation.

cat("Setting 16: Small n Regime (CLT Breakdown)\n")
cat("--------------------------------------------\n")
setting16_start <- Sys.time()

library(pppower)

# Parameters — focus on very small n
n_values_16 <- c(15, 20, 25, 30, 50, 100)
N_values_16 <- c(200, 500)
rho_values_16 <- c(0.5, 0.7, 0.9)
delta_16 <- 0.3   # slightly larger delta so power is non-trivial even at small n
sigma_Y2_16 <- 1
sigma_f2_16 <- 1
alpha_16 <- 0.05
R_16 <- 2000  # more reps for tighter MC estimates at small n

grid_16 <- expand.grid(
  rho = rho_values_16,
  N = N_values_16,
  n = n_values_16,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_16 <- nrow(grid_16)

results_small_n <- vector("list", total_16)

for (i in seq_len(total_16)) {
  row <- grid_16[i, ]
  progress_bar(i, total_16, "Setting 16",
               sprintf("rho=%.1f, N=%d, n=%d: ", row$rho, row$N, row$n),
               start_time = setting16_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_16 * sigma_f2_16)
  Sigma <- matrix(c(sigma_Y2_16, cov_Yf, cov_Yf, sigma_f2_16), 2, 2)
  L <- chol(Sigma)

  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2_16)

  power_theo_ppi <- theo_power_ppi_mean(delta_16, row$n, row$N,
                                         sigma_Y2_16, sigma_f2_16, cov_Yf, alpha_16)
  power_theo_classical <- theo_power_classical_onesample(delta_16, row$n,
                                                          sigma_Y2_16, alpha_16)

  rej_ppi <- 0
  rej_classical <- 0

  for (r in seq_len(R_16)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1] + delta_16
    f_L <- Z_L[, 2] + delta_16

    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    f_U <- Z_U[, 2] + delta_16

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha_16,
                                  lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha_16)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_small_n[[i]] <- data.frame(
    setting = "Small n",
    rho = row$rho,
    n = row$n,
    N = row$N,
    delta = delta_16,
    power_theo_ppi = power_theo_ppi,
    power_emp_ppi = rej_ppi / R_16,
    power_theo_classical = power_theo_classical,
    power_emp_classical = rej_classical / R_16,
    stringsAsFactors = FALSE
  )
}

results_small_n <- do.call(rbind, results_small_n)

cat(sprintf("Setting 16 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting16_start, units = "secs")))
# Show discrepancy at small n
small <- results_small_n[results_small_n$n <= 25, ]
cat(sprintf("  Max |theo - emp| at n <= 25: %.4f\n",
            max(abs(small$power_theo_ppi - small$power_emp_ppi))))
large <- results_small_n[results_small_n$n >= 50, ]
cat(sprintf("  Max |theo - emp| at n >= 50: %.4f\n\n",
            max(abs(large$power_theo_ppi - large$power_emp_ppi))))
