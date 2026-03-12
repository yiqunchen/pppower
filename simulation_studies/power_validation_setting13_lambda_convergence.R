#!/usr/bin/env Rscript
# =============================================================================
# Setting 13: Plugin Lambda Convergence
# =============================================================================
# Validates that the plugin lambda estimator converges to the oracle lambda*
# as the labeled sample size n increases.

cat("Setting 13: Plugin Lambda Convergence\n")
cat("---------------------------------------\n")
setting13_start <- Sys.time()

library(pppower)

rho_vals_13 <- c(0.5, 0.7, 0.9)
N_13        <- 1000
n_vals_13   <- c(20, 50, 100, 200, 500)
R_13        <- 1000
sigma_Y2_13 <- 1.0
sigma_f2_13 <- 1.0
delta_13    <- 0.2

grid_13 <- expand.grid(
  rho = rho_vals_13,
  n   = n_vals_13,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_13 <- nrow(grid_13)

results_lambda_convergence <- vector("list", total_13)

for (i in seq_len(total_13)) {
  row <- grid_13[i, ]
  progress_bar(i, total_13, "Setting 13",
               sprintf("rho=%.1f, n=%d: ", row$rho, row$n),
               start_time = setting13_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_13 * sigma_f2_13)
  Sigma  <- matrix(c(sigma_Y2_13, cov_Yf, cov_Yf, sigma_f2_13), 2, 2)
  L      <- chol(Sigma)

  lambda_oracle <- cov_Yf / ((1 + row$n / N_13) * sigma_f2_13)

  lambda_hat_vec <- numeric(R_13)

  for (r in seq_len(R_13)) {
    set.seed(13000 + (i - 1) * R_13 + r)

    # Labeled data
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1] + delta_13
    f_L <- Z_L[, 2] + delta_13

    # Unlabeled data
    Z_U <- matrix(rnorm(N_13 * 2), N_13, 2) %*% L
    f_U <- Z_U[, 2] + delta_13

    # Plugin lambda estimate: Cov(Y,f)_hat / ((1 + n/N) * Var(f)_hat)
    cov_hat  <- cov(Y_L, f_L)
    varf_hat <- var(f_L)

    if (varf_hat > 0) {
      lambda_hat_vec[r] <- cov_hat / ((1 + row$n / N_13) * varf_hat)
    } else {
      lambda_hat_vec[r] <- 0
    }
  }

  results_lambda_convergence[[i]] <- data.frame(
    rho = row$rho,
    n = row$n,
    N = N_13,
    lambda_oracle = lambda_oracle,
    lambda_hat_mean = mean(lambda_hat_vec),
    lambda_hat_sd = sd(lambda_hat_vec),
    lambda_hat_median = median(lambda_hat_vec),
    bias = mean(lambda_hat_vec) - lambda_oracle,
    rmse = sqrt(mean((lambda_hat_vec - lambda_oracle)^2)),
    stringsAsFactors = FALSE
  )
}

results_lambda_convergence <- do.call(rbind, results_lambda_convergence)

cat(sprintf("Setting 13 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting13_start, units = "secs")))
cat("  Lambda convergence summary (n = 500):\n")
subset_2k <- results_lambda_convergence[results_lambda_convergence$n == 500, ]
for (j in seq_len(nrow(subset_2k))) {
  r <- subset_2k[j, ]
  cat(sprintf("    rho=%.1f: oracle=%.4f, mean_hat=%.4f, SD=%.4f, RMSE=%.4f\n",
              r$rho, r$lambda_oracle, r$lambda_hat_mean, r$lambda_hat_sd, r$rmse))
}
