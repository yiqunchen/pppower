#!/usr/bin/env Rscript
# =============================================================================
# Setting 12: Type I Error Calibration
# =============================================================================
# Re-runs Settings 1 and 2 DGPs with Delta = 0 (null hypothesis true) to
# verify that rejection rates are calibrated at the nominal level alpha = 0.05.

cat("Setting 12: Type I Error Calibration\n")
cat("--------------------------------------\n")
setting12_start <- Sys.time()

library(pppower)

alpha_12 <- 0.05
R_12     <- 2000

# ---- Setting 1 DGP (Gaussian mean) under H0 ----
rho_vals_12 <- c(0.5, 0.7, 0.9)
n_vals_12   <- c(50, 100, 200, 500)
N_vals_12   <- c(1000, 5000)
sigma_Y2_12 <- 1.0
sigma_f2_12 <- 1.0

grid_12a <- expand.grid(
  rho = rho_vals_12,
  n   = n_vals_12,
  N   = N_vals_12,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_12a <- nrow(grid_12a)

results_type1_cont <- vector("list", total_12a)

for (i in seq_len(total_12a)) {
  row <- grid_12a[i, ]
  progress_bar(i, total_12a, "Setting 12a",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N))

  cov_Yf <- row$rho * sqrt(sigma_Y2_12 * sigma_f2_12)
  Sigma  <- matrix(c(sigma_Y2_12, cov_Yf, cov_Yf, sigma_f2_12), 2, 2)
  L      <- chol(Sigma)
  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2_12)

  rej_ppi      <- 0
  rej_classical <- 0

  for (r in seq_len(R_12)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1]   # Delta = 0
    f_L <- Z_L[, 2]

    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    f_U <- Z_U[, 2]

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0,
                                  alpha = alpha_12, lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha_12)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_type1_cont[[i]] <- data.frame(
    setting = "Gaussian mean (H0)",
    rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12,
    rejection_classical = rej_classical / R_12,
    stringsAsFactors = FALSE
  )
}

results_type1_cont <- do.call(rbind, results_type1_cont)

# ---- Setting 2 DGP (Binary mean) under H0 ----
p0_12     <- 0.3
sens_spec_12 <- list(c(0.70, 0.70), c(0.85, 0.85), c(0.95, 0.95))
n_bin_12  <- c(100, 200, 500, 1000)
N_bin_12  <- c(5000, 10000)

grid_12b <- expand.grid(
  ss_ix = seq_along(sens_spec_12),
  n     = n_bin_12,
  N     = N_bin_12,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_12b <- nrow(grid_12b)

results_type1_bin <- vector("list", total_12b)

for (i in seq_len(total_12b)) {
  row  <- grid_12b[i, ]
  ss   <- sens_spec_12[[row$ss_ix]]
  sens_i <- ss[1]; spec_i <- ss[2]

  progress_bar(i, total_12b, "Setting 12b",
               sprintf("sens=%.2f, n=%d, N=%d: ", sens_i, row$n, row$N))

  moments <- binary_moments_from_sens_spec(p = p0_12, sens = sens_i, spec = spec_i)
  lambda_oracle <- moments$cov_y_f / ((1 + row$n / row$N) * moments$sigma_f2)

  rej_ppi <- 0
  rej_classical <- 0

  for (r in seq_len(R_12)) {
    Y_L <- rbinom(row$n, 1, p0_12)
    f_L <- ifelse(Y_L == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))

    Y_U_hidden <- rbinom(row$N, 1, p0_12)
    f_U <- ifelse(Y_U_hidden == 1, rbinom(row$N, 1, sens_i),
                  rbinom(row$N, 1, 1 - spec_i))

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = p0_12,
                                  alpha = alpha_12, lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_mean_test(Y_L, theta0 = p0_12, alpha = alpha_12)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_type1_bin[[i]] <- data.frame(
    setting = "Binary mean (H0)",
    sens = sens_i, spec = spec_i,
    n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12,
    rejection_classical = rej_classical / R_12,
    stringsAsFactors = FALSE
  )
}

results_type1_bin <- do.call(rbind, results_type1_bin)

# Combine
results_type1_error <- list(
  continuous = results_type1_cont,
  binary     = results_type1_bin
)

# Check calibration: 95% CI for alpha=0.05 with R=2000 is approx [0.040, 0.060]
ci_lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / R_12)
ci_hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / R_12)

calibrated_cont <- all(
  results_type1_cont$rejection_ppi >= ci_lo &
  results_type1_cont$rejection_ppi <= ci_hi
)
calibrated_bin <- all(
  results_type1_bin$rejection_ppi >= ci_lo &
  results_type1_bin$rejection_ppi <= ci_hi
)

cat(sprintf("Setting 12 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting12_start, units = "secs")))
cat(sprintf("  95%% CI for nominal alpha=0.05: [%.3f, %.3f]\n", ci_lo, ci_hi))
cat(sprintf("  Continuous: all PPI++ Type I errors within CI? %s\n",
            ifelse(calibrated_cont, "YES", "NO")))
cat(sprintf("  Binary:     all PPI++ Type I errors within CI? %s\n",
            ifelse(calibrated_bin, "YES", "NO")))
cat(sprintf("  Cont range:  [%.4f, %.4f]\n",
            min(results_type1_cont$rejection_ppi),
            max(results_type1_cont$rejection_ppi)))
cat(sprintf("  Bin range:   [%.4f, %.4f]\n",
            min(results_type1_bin$rejection_ppi),
            max(results_type1_bin$rejection_ppi)))
