#!/usr/bin/env Rscript
# =============================================================================
# Setting 14: Plugin Lambda Power Validation
# =============================================================================
# All prior settings use oracle lambda* (true population moments).
# This setting validates that the analytical power formula still holds when
# lambda is estimated from data (plugin), i.e., the practical scenario.
#
# Design: For each (rho, n, N), run R replications:
#   (a) Oracle lambda:  lambda* = cov_Yf / ((1 + n/N) * sigma_f2)
#   (b) Plugin lambda:  lambda_hat = cov_hat(Y,f) / ((1 + n/N) * var_hat(f))
# Record empirical power for both and compare to the same theoretical formula.

cat("Setting 14: Plugin Lambda Power Validation\n")
cat("--------------------------------------------\n")
setting14_start <- Sys.time()

library(pppower)

# Parameters
n_values_14 <- c(20, 40, 60, 80, 100)
N_values_14 <- c(200, 500)
rho_values_14 <- c(0.5, 0.7, 0.9)
delta_14 <- 0.2
sigma_Y2_14 <- 1
sigma_f2_14 <- 1
alpha_14 <- 0.05
R_14 <- 1000

grid_14 <- expand.grid(
  rho = rho_values_14,
  N = N_values_14,
  n = n_values_14,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_14 <- nrow(grid_14)

results_plugin_lambda <- vector("list", total_14)

for (i in seq_len(total_14)) {
  row <- grid_14[i, ]
  progress_bar(i, total_14, "Setting 14",
               sprintf("rho=%.1f, N=%d, n=%d: ", row$rho, row$N, row$n),
               start_time = setting14_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_14 * sigma_f2_14)
  Sigma <- matrix(c(sigma_Y2_14, cov_Yf, cov_Yf, sigma_f2_14), 2, 2)
  L <- chol(Sigma)

  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2_14)

  # Theoretical power (same formula for both — the question is whether
  # plugin achieves the same empirical power as oracle)
  power_theo <- theo_power_ppi_mean(delta_14, row$n, row$N,
                                     sigma_Y2_14, sigma_f2_14, cov_Yf, alpha_14)

  rej_oracle <- 0
  rej_plugin <- 0

  for (r in seq_len(R_14)) {
    # Generate bivariate normal (Y, f) under H1
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1] + delta_14
    f_L <- Z_L[, 2] + delta_14

    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    f_U <- Z_U[, 2] + delta_14

    # (a) Oracle lambda
    test_oracle <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0,
                                     alpha = alpha_14,
                                     lambda_oracle = lambda_oracle)
    rej_oracle <- rej_oracle + test_oracle$reject

    # (b) Plugin lambda (estimated from labeled data)
    cov_hat <- cov(Y_L, f_L)
    varf_hat <- var(f_L)
    if (varf_hat > 0) {
      lambda_plugin <- cov_hat / ((1 + row$n / row$N) * varf_hat)
    } else {
      lambda_plugin <- 0
    }
    test_plugin <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0,
                                     alpha = alpha_14,
                                     lambda_oracle = lambda_plugin)
    rej_plugin <- rej_plugin + test_plugin$reject
  }

  results_plugin_lambda[[i]] <- data.frame(
    setting = "Plugin Lambda",
    rho = row$rho,
    n = row$n,
    N = row$N,
    delta = delta_14,
    power_theo = power_theo,
    power_emp_oracle = rej_oracle / R_14,
    power_emp_plugin = rej_plugin / R_14,
    stringsAsFactors = FALSE
  )
}

results_plugin_lambda <- do.call(rbind, results_plugin_lambda)

cat(sprintf("Setting 14 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting14_start, units = "secs")))
# Summary
cat("  Max |theo - oracle|:  ",
    max(abs(results_plugin_lambda$power_theo - results_plugin_lambda$power_emp_oracle)), "\n")
cat("  Max |theo - plugin|:  ",
    max(abs(results_plugin_lambda$power_theo - results_plugin_lambda$power_emp_plugin)), "\n")
cat("  Max |oracle - plugin|:",
    max(abs(results_plugin_lambda$power_emp_oracle - results_plugin_lambda$power_emp_plugin)), "\n\n")
