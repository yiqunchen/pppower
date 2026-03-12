#!/usr/bin/env Rscript
# =============================================================================
# Setting 15: Power vs Effect Size (Delta)
# =============================================================================
# Standard practice: show power as a function of delta at fixed n.
# Sweeps delta from 0 to 0.5 at several fixed n values.

cat("Setting 15: Power vs Effect Size (Delta)\n")
cat("------------------------------------------\n")
setting15_start <- Sys.time()

library(pppower)

# Parameters
delta_values_15 <- seq(0.0, 0.5, by = 0.025)
n_values_15 <- c(20, 40, 60, 80, 100)
N_15 <- 500
rho_values_15 <- c(0.5, 0.7, 0.9)
sigma_Y2_15 <- 1
sigma_f2_15 <- 1
alpha_15 <- 0.05
R_15 <- 1000

grid_15 <- expand.grid(
  rho = rho_values_15,
  n = n_values_15,
  delta = delta_values_15,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_15 <- nrow(grid_15)

results_power_vs_delta <- vector("list", total_15)

for (i in seq_len(total_15)) {
  row <- grid_15[i, ]
  progress_bar(i, total_15, "Setting 15",
               sprintf("rho=%.1f, n=%d, delta=%.3f: ", row$rho, row$n, row$delta),
               start_time = setting15_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_15 * sigma_f2_15)
  Sigma <- matrix(c(sigma_Y2_15, cov_Yf, cov_Yf, sigma_f2_15), 2, 2)
  L <- chol(Sigma)

  lambda_oracle <- cov_Yf / ((1 + row$n / N_15) * sigma_f2_15)

  # Theoretical power
  power_theo_ppi <- theo_power_ppi_mean(row$delta, row$n, N_15,
                                         sigma_Y2_15, sigma_f2_15, cov_Yf, alpha_15)
  power_theo_classical <- theo_power_classical_onesample(row$delta, row$n,
                                                          sigma_Y2_15, alpha_15)

  # Monte Carlo
  rej_ppi <- 0
  rej_classical <- 0

  for (r in seq_len(R_15)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1] + row$delta
    f_L <- Z_L[, 2] + row$delta

    Z_U <- matrix(rnorm(N_15 * 2), N_15, 2) %*% L
    f_U <- Z_U[, 2] + row$delta

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha_15,
                                  lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha_15)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_power_vs_delta[[i]] <- data.frame(
    setting = "Power vs Delta",
    rho = row$rho,
    n = row$n,
    N = N_15,
    delta = row$delta,
    power_theo_ppi = power_theo_ppi,
    power_emp_ppi = rej_ppi / R_15,
    power_theo_classical = power_theo_classical,
    power_emp_classical = rej_classical / R_15,
    stringsAsFactors = FALSE
  )
}

results_power_vs_delta <- do.call(rbind, results_power_vs_delta)

cat(sprintf("Setting 15 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting15_start, units = "secs")))
cat("  Max |theo - emp| PPI++:",
    max(abs(results_power_vs_delta$power_theo_ppi - results_power_vs_delta$power_emp_ppi)), "\n\n")
