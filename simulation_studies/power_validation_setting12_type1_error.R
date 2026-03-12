#!/usr/bin/env Rscript
# =============================================================================
# Setting 12: Type I Error Calibration
# =============================================================================
# Re-runs Settings 1-8 DGPs with Delta = 0 (null hypothesis true) to verify
# that rejection rates are calibrated at the nominal level alpha = 0.05.
# Uses R = 2000 for tighter MC estimates.

cat("Setting 12: Type I Error Calibration\n")
cat("--------------------------------------\n")
setting12_start <- Sys.time()

library(pppower)

alpha_12 <- 0.05
R_12     <- 2000

# =============================================================================
# 12a: Setting 1 DGP — Gaussian mean under H0
# =============================================================================
cat("  12a: Gaussian mean (Setting 1 DGP)...\n")
rho_vals_12 <- c(0.5, 0.7, 0.9)
n_vals_12   <- c(20, 40, 60, 80, 100)
N_vals_12   <- c(200, 500)
sigma_Y2_12 <- 1.0
sigma_f2_12 <- 1.0

grid_12a <- expand.grid(
  rho = rho_vals_12, n = n_vals_12, N = N_vals_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12a <- nrow(grid_12a)
results_type1_cont <- vector("list", total_12a)

for (i in seq_len(total_12a)) {
  row <- grid_12a[i, ]
  progress_bar(i, total_12a, "Setting 12a",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N),
               start_time = setting12_start)
  cov_Yf <- row$rho * sqrt(sigma_Y2_12 * sigma_f2_12)
  Sigma  <- matrix(c(sigma_Y2_12, cov_Yf, cov_Yf, sigma_f2_12), 2, 2)
  L      <- chol(Sigma)
  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2_12)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    test_ppi <- run_ppi_mean_test(Z_L[,1], Z_L[,2], Z_U[,2], theta0 = 0,
                                  alpha = alpha_12, lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(Z_L[,1], theta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_cont[[i]] <- data.frame(
    setting = "1_Gaussian_mean", rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_cont <- do.call(rbind, results_type1_cont)

# =============================================================================
# 12b: Setting 2 DGP — Binary mean under H0
# =============================================================================
cat("  12b: Binary mean (Setting 2 DGP)...\n")
p0_12        <- 0.3
sens_spec_12 <- list(c(0.70, 0.70), c(0.85, 0.85), c(0.95, 0.95))
n_bin_12     <- c(20, 40, 60, 80, 100)
N_bin_12     <- c(200, 500)

grid_12b <- expand.grid(
  ss_ix = seq_along(sens_spec_12), n = n_bin_12, N = N_bin_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12b <- nrow(grid_12b)
results_type1_bin <- vector("list", total_12b)

for (i in seq_len(total_12b)) {
  row <- grid_12b[i, ]
  ss <- sens_spec_12[[row$ss_ix]]; sens_i <- ss[1]; spec_i <- ss[2]
  progress_bar(i, total_12b, "Setting 12b",
               sprintf("sens=%.2f, n=%d, N=%d: ", sens_i, row$n, row$N),
               start_time = setting12_start)
  moments <- binary_moments_from_sens_spec(p = p0_12, sens = sens_i, spec = spec_i)
  lambda_oracle <- moments$cov_y_f / ((1 + row$n / row$N) * moments$sigma_f2)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    Y_L <- rbinom(row$n, 1, p0_12)
    f_L <- ifelse(Y_L == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))
    Y_U_h <- rbinom(row$N, 1, p0_12)
    f_U <- ifelse(Y_U_h == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))
    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = p0_12,
                                  alpha = alpha_12, lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(Y_L, theta0 = p0_12, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_bin[[i]] <- data.frame(
    setting = "2_Binary_mean", sens = sens_i, spec = spec_i,
    n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_bin <- do.call(rbind, results_type1_bin)

# =============================================================================
# 12c: Setting 3 DGP — Two-sample t-test (continuous) under H0
# =============================================================================
cat("  12c: Two-sample t-test continuous (Setting 3 DGP)...\n")
# Use a subset of the grid to keep runtime manageable
n_tt_12  <- c(20, 40, 60, 80, 100)
N_tt_12  <- c(200, 500)
rho_tt_12 <- c(0.5, 0.7, 0.9)

grid_12c <- expand.grid(
  rho = rho_tt_12, n = n_tt_12, N = N_tt_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12c <- nrow(grid_12c)
results_type1_ttest <- vector("list", total_12c)

for (i in seq_len(total_12c)) {
  row <- grid_12c[i, ]
  progress_bar(i, total_12c, "Setting 12c",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N),
               start_time = setting12_start)
  cov_Yf <- row$rho * sqrt(sigma_Y2_12 * sigma_f2_12)
  Sigma  <- matrix(c(sigma_Y2_12, cov_Yf, cov_Yf, sigma_f2_12), 2, 2)
  L      <- chol(Sigma)
  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2_12)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    # Group A: mean 0 under H0
    Z_A <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Z_AU <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    # Group B: mean 0 under H0
    Z_B <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Z_BU <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    test_ppi <- run_ppi_ttest(
      Z_A[,1], Z_A[,2], Z_AU[,2],
      Z_B[,1], Z_B[,2], Z_BU[,2],
      delta0 = 0, alpha = alpha_12,
      lambda_A_oracle = lambda_oracle, lambda_B_oracle = lambda_oracle
    )
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_ttest(Z_A[,1], Z_B[,1], delta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_ttest[[i]] <- data.frame(
    setting = "3_Ttest_continuous", rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_ttest <- do.call(rbind, results_type1_ttest)

# =============================================================================
# 12d: Setting 4 DGP — Two-sample proportion test (binary) under H0
# =============================================================================
cat("  12d: Two-sample proportion test binary (Setting 4 DGP)...\n")
n_prop_12 <- c(20, 40, 60, 80, 100)
N_prop_12 <- c(200, 500)

grid_12d <- expand.grid(
  ss_ix = seq_along(sens_spec_12), n = n_prop_12, N = N_prop_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12d <- nrow(grid_12d)
results_type1_prop <- vector("list", total_12d)

for (i in seq_len(total_12d)) {
  row <- grid_12d[i, ]
  ss <- sens_spec_12[[row$ss_ix]]; sens_i <- ss[1]; spec_i <- ss[2]
  progress_bar(i, total_12d, "Setting 12d",
               sprintf("sens=%.2f, n=%d, N=%d: ", sens_i, row$n, row$N),
               start_time = setting12_start)
  moments_p <- binary_moments_from_sens_spec(p = p0_12, sens = sens_i, spec = spec_i)
  lambda_oracle <- moments_p$cov_y_f / ((1 + row$n / row$N) * moments_p$sigma_f2)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    # Both groups have same prevalence under H0
    Y_A <- rbinom(row$n, 1, p0_12)
    f_A_L <- ifelse(Y_A == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))
    Y_A_U <- rbinom(row$N, 1, p0_12)
    f_A_U <- ifelse(Y_A_U == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))
    Y_B <- rbinom(row$n, 1, p0_12)
    f_B_L <- ifelse(Y_B == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))
    Y_B_U <- rbinom(row$N, 1, p0_12)
    f_B_U <- ifelse(Y_B_U == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))
    test_ppi <- run_ppi_ttest(
      Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
      delta0 = 0, alpha = alpha_12,
      lambda_A_oracle = lambda_oracle, lambda_B_oracle = lambda_oracle
    )
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_ttest(Y_A, Y_B, delta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_prop[[i]] <- data.frame(
    setting = "4_Ttest_binary", sens = sens_i, spec = spec_i,
    n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_prop <- do.call(rbind, results_type1_prop)

# =============================================================================
# 12e: Setting 5 DGP — Paired t-test (continuous) under H0
# =============================================================================
cat("  12e: Paired t-test continuous (Setting 5 DGP)...\n")
n_paired_12 <- c(20, 40, 60, 80, 100)
N_paired_12 <- c(200, 500)
sigma_D2_12 <- 1.0
sigma_fD2_12 <- 1.0

grid_12e <- expand.grid(
  rho = rho_vals_12, n = n_paired_12, N = N_paired_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12e <- nrow(grid_12e)
results_type1_paired <- vector("list", total_12e)

for (i in seq_len(total_12e)) {
  row <- grid_12e[i, ]
  progress_bar(i, total_12e, "Setting 12e",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N),
               start_time = setting12_start)
  cov_DfD <- row$rho * sqrt(sigma_D2_12 * sigma_fD2_12)
  Sigma_p <- matrix(c(sigma_D2_12, cov_DfD, cov_DfD, sigma_fD2_12), 2, 2)
  L_p     <- chol(Sigma_p)
  lambda_oracle <- cov_DfD / ((1 + row$n / row$N) * sigma_fD2_12)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L_p
    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L_p
    # Delta = 0 under H0
    test_ppi <- run_ppi_mean_test(Z_L[,1], Z_L[,2], Z_U[,2], theta0 = 0,
                                  alpha = alpha_12, lambda_oracle = lambda_oracle)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(Z_L[,1], theta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_paired[[i]] <- data.frame(
    setting = "5_Paired_continuous", rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_paired <- do.call(rbind, results_type1_paired)

# =============================================================================
# 12f: Setting 6 DGP — Paired proportion test (binary) under H0
# =============================================================================
cat("  12f: Paired proportion test binary (Setting 6 DGP)...\n")
rho_within_12 <- 0.3
n_pairedbin_12 <- c(20, 40, 60, 80, 100)
N_pairedbin_12 <- c(200, 500)

grid_12f <- expand.grid(
  ss_ix = seq_along(sens_spec_12), n = n_pairedbin_12, N = N_pairedbin_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12f <- nrow(grid_12f)
results_type1_paired_bin <- vector("list", total_12f)

for (i in seq_len(total_12f)) {
  row <- grid_12f[i, ]
  ss <- sens_spec_12[[row$ss_ix]]; sens_i <- ss[1]; spec_i <- ss[2]
  progress_bar(i, total_12f, "Setting 12f",
               sprintf("sens=%.2f, n=%d, N=%d: ", sens_i, row$n, row$N),
               start_time = setting12_start)

  # Under H0: p_A = p_B = p0
  # Joint Bernoulli with correlation rho_within
  p_both <- p0_12^2 + rho_within_12 * sqrt(p0_12 * (1 - p0_12) *
                                             p0_12 * (1 - p0_12))

  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    # Generate correlated pairs under H0
    U <- runif(row$n)
    Y_A <- as.integer(U < p0_12)
    Y_B <- ifelse(Y_A == 1,
                  rbinom(row$n, 1, p_both / p0_12),
                  rbinom(row$n, 1, (p0_12 - p_both) / (1 - p0_12)))
    D_L <- Y_A - Y_B
    f_A_L <- ifelse(Y_A == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))
    f_B_L <- ifelse(Y_B == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))
    fD_L <- f_A_L - f_B_L

    U_u <- runif(row$N)
    Y_A_U <- as.integer(U_u < p0_12)
    Y_B_U <- ifelse(Y_A_U == 1,
                    rbinom(row$N, 1, p_both / p0_12),
                    rbinom(row$N, 1, (p0_12 - p_both) / (1 - p0_12)))
    f_A_U <- ifelse(Y_A_U == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))
    f_B_U <- ifelse(Y_B_U == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))
    fD_U <- f_A_U - f_B_U

    # Plugin lambda from labeled data
    if (var(fD_L) > 0) {
      lambda_hat <- cov(D_L, fD_L) / ((1 + row$n / row$N) * var(fD_L))
    } else {
      lambda_hat <- 0
    }
    test_ppi <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0,
                                  alpha = alpha_12, lambda_oracle = lambda_hat)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(D_L, theta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_paired_bin[[i]] <- data.frame(
    setting = "6_Paired_binary", sens = sens_i, spec = spec_i,
    n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_paired_bin <- do.call(rbind, results_type1_paired_bin)

# =============================================================================
# 12g: Setting 7 DGP — Log-normal under H0
# =============================================================================
cat("  12g: Log-normal (Setting 7 DGP)...\n")
mu_log_12  <- 0
sigma_log_12 <- 0.5
# Under H0 the true mean is exp(mu + sigma^2/2)
theta0_lognorm <- exp(mu_log_12 + sigma_log_12^2 / 2)
n_ln_12 <- c(20, 40, 60, 80, 100)
N_ln_12 <- c(200, 500)

grid_12g <- expand.grid(
  rho = rho_vals_12, n = n_ln_12, N = N_ln_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12g <- nrow(grid_12g)
results_type1_lognorm <- vector("list", total_12g)

for (i in seq_len(total_12g)) {
  row <- grid_12g[i, ]
  progress_bar(i, total_12g, "Setting 12g",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N),
               start_time = setting12_start)
  Sigma_ln <- matrix(c(1, row$rho, row$rho, 1), 2, 2) * sigma_log_12^2
  L_ln <- chol(Sigma_ln)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L_ln + mu_log_12
    Y_L <- exp(Z_L[, 1])
    f_L <- exp(Z_L[, 2])
    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L_ln + mu_log_12
    f_U <- exp(Z_U[, 2])
    # Plugin lambda
    if (var(f_L) > 0) {
      lam <- cov(Y_L, f_L) / ((1 + row$n / row$N) * var(f_L))
    } else {
      lam <- 0
    }
    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = theta0_lognorm,
                                  alpha = alpha_12, lambda_oracle = lam)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(Y_L, theta0 = theta0_lognorm, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_lognorm[[i]] <- data.frame(
    setting = "7_LogNormal", rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_lognorm <- do.call(rbind, results_type1_lognorm)

# =============================================================================
# 12h: Setting 8 DGP — t-distributed (df=5) under H0
# =============================================================================
cat("  12h: t-distribution df=5 (Setting 8 DGP)...\n")
df_t <- 5
sigma_Y2_t <- df_t / (df_t - 2)   # 5/3
n_t_12 <- c(20, 40, 60, 80, 100)
N_t_12 <- c(200, 500)

grid_12h <- expand.grid(
  rho = rho_vals_12, n = n_t_12, N = N_t_12,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
)
total_12h <- nrow(grid_12h)
results_type1_tdist <- vector("list", total_12h)

for (i in seq_len(total_12h)) {
  row <- grid_12h[i, ]
  progress_bar(i, total_12h, "Setting 12h",
               sprintf("rho=%.1f, n=%d, N=%d: ", row$rho, row$n, row$N),
               start_time = setting12_start)
  Sigma_t <- matrix(c(1, row$rho, row$rho, 1), 2, 2)
  L_t <- chol(Sigma_t)
  rej_ppi <- 0; rej_cl <- 0
  for (r in seq_len(R_12)) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L_t
    U_L <- pnorm(Z_L)
    Y_L <- qt(U_L[, 1], df = df_t)   # mean 0 under H0
    f_L <- qt(U_L[, 2], df = df_t)
    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L_t
    U_U <- pnorm(Z_U)
    f_U <- qt(U_U[, 2], df = df_t)
    # Plugin lambda
    if (var(f_L) > 0) {
      lam <- cov(Y_L, f_L) / ((1 + row$n / row$N) * var(f_L))
    } else {
      lam <- 0
    }
    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0,
                                  alpha = alpha_12, lambda_oracle = lam)
    rej_ppi <- rej_ppi + test_ppi$reject
    test_cl <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha_12)
    rej_cl <- rej_cl + test_cl$reject
  }
  results_type1_tdist[[i]] <- data.frame(
    setting = "8_t_dist", rho = row$rho, n = row$n, N = row$N,
    rejection_ppi = rej_ppi / R_12, rejection_classical = rej_cl / R_12,
    stringsAsFactors = FALSE
  )
}
results_type1_tdist <- do.call(rbind, results_type1_tdist)

# =============================================================================
# Combine all results
# =============================================================================
results_type1_error <- list(
  continuous   = results_type1_cont,
  binary       = results_type1_bin,
  ttest_cont   = results_type1_ttest,
  ttest_bin    = results_type1_prop,
  paired_cont  = results_type1_paired,
  paired_bin   = results_type1_paired_bin,
  lognormal    = results_type1_lognorm,
  tdist        = results_type1_tdist
)

# Check calibration: 95% CI for alpha=0.05 with R=2000 is approx [0.040, 0.060]
ci_lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / R_12)
ci_hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / R_12)

cat(sprintf("\nSetting 12 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting12_start, units = "secs")))
cat(sprintf("  95%% CI for nominal alpha=0.05: [%.3f, %.3f]\n", ci_lo, ci_hi))

for (nm in names(results_type1_error)) {
  df <- results_type1_error[[nm]]
  rng <- range(df$rejection_ppi)
  all_ok <- all(df$rejection_ppi >= ci_lo & df$rejection_ppi <= ci_hi)
  cat(sprintf("  %-15s: PPI++ range [%.4f, %.4f] — %s\n",
              nm, rng[1], rng[2], ifelse(all_ok, "ALL OK", "SOME OUTSIDE CI")))
}
