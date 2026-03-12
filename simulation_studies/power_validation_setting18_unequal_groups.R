#!/usr/bin/env Rscript
# =============================================================================
# Setting 18: Unequal Group Sizes (Two-Sample)
# =============================================================================
# Settings 3-4 always use n_A = n_B, N_A = N_B. In practice designs are
# often unbalanced. This setting validates the analytical power formula for
# asymmetric two-sample designs.
#
# Design: fix total labeled budget n_A + n_B = 200, vary allocation ratio.
# Similarly for unlabeled budget.

cat("Setting 18: Unequal Group Sizes (Two-Sample)\n")
cat("----------------------------------------------\n")
setting18_start <- Sys.time()

library(pppower)

# Parameters
# Allocation ratios: n_A / n_B
alloc_ratios <- c(1, 1.5, 2, 3, 4)  # 1 = balanced
n_total <- 100   # total labeled
N_total <- 600   # total unlabeled
rho_values_18 <- c(0.5, 0.7, 0.9)
delta_18 <- 0.3
sigma_Y2_18 <- 1
sigma_f2_18 <- 1
alpha_18 <- 0.05
R_18 <- 1000

grid_18 <- expand.grid(
  rho = rho_values_18,
  alloc = alloc_ratios,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_18 <- nrow(grid_18)

results_unequal <- vector("list", total_18)

for (i in seq_len(total_18)) {
  row <- grid_18[i, ]
  # Split budgets by allocation ratio
  n_A <- round(n_total * row$alloc / (1 + row$alloc))
  n_B <- n_total - n_A
  N_A <- round(N_total * row$alloc / (1 + row$alloc))
  N_B <- N_total - N_A

  progress_bar(i, total_18, "Setting 18",
               sprintf("rho=%.1f, n_A=%d, n_B=%d: ", row$rho, n_A, n_B),
               start_time = setting18_start)

  cov_Yf <- row$rho * sqrt(sigma_Y2_18 * sigma_f2_18)
  Sigma <- matrix(c(sigma_Y2_18, cov_Yf, cov_Yf, sigma_f2_18), 2, 2)
  L <- chol(Sigma)

  lambda_oracle_A <- cov_Yf / ((1 + n_A / N_A) * sigma_f2_18)
  lambda_oracle_B <- cov_Yf / ((1 + n_B / N_B) * sigma_f2_18)

  # Theoretical power
  power_theo_ppi <- theo_power_ppi_twosample(
    delta_18, n_A, n_B, N_A, N_B,
    sigma_Y2_18, sigma_f2_18, cov_Yf,
    sigma_Y2_18, sigma_f2_18, cov_Yf, alpha_18
  )
  # Classical: Welch t-test approximation uses harmonic-mean-like df
  # pwr uses equal-n formula, so we use the smaller group for conservative estimate
  power_theo_classical <- theo_power_classical_twosample(delta_18, min(n_A, n_B),
                                                          sigma_Y2_18, alpha_18)

  rej_ppi <- 0
  rej_classical <- 0

  for (r in seq_len(R_18)) {
    # Group A: mean = delta/2
    Z_A_L <- matrix(rnorm(n_A * 2), n_A, 2) %*% L
    Y_A <- Z_A_L[, 1] + delta_18 / 2
    f_A_L <- Z_A_L[, 2] + delta_18 / 2
    Z_A_U <- matrix(rnorm(N_A * 2), N_A, 2) %*% L
    f_A_U <- Z_A_U[, 2] + delta_18 / 2

    # Group B: mean = -delta/2
    Z_B_L <- matrix(rnorm(n_B * 2), n_B, 2) %*% L
    Y_B <- Z_B_L[, 1] - delta_18 / 2
    f_B_L <- Z_B_L[, 2] - delta_18 / 2
    Z_B_U <- matrix(rnorm(N_B * 2), N_B, 2) %*% L
    f_B_U <- Z_B_U[, 2] - delta_18 / 2

    test_ppi <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                              delta0 = 0, alpha = alpha_18,
                              lambda_A_oracle = lambda_oracle_A,
                              lambda_B_oracle = lambda_oracle_B)
    rej_ppi <- rej_ppi + test_ppi$reject

    test_cl <- classical_ttest(Y_A, Y_B, delta0 = 0, alpha = alpha_18)
    rej_classical <- rej_classical + test_cl$reject
  }

  results_unequal[[i]] <- data.frame(
    setting = "Unequal Groups",
    rho = row$rho,
    alloc_ratio = row$alloc,
    n_A = n_A,
    n_B = n_B,
    N_A = N_A,
    N_B = N_B,
    delta = delta_18,
    power_theo_ppi = power_theo_ppi,
    power_emp_ppi = rej_ppi / R_18,
    power_emp_classical = rej_classical / R_18,
    stringsAsFactors = FALSE
  )
}

results_unequal <- do.call(rbind, results_unequal)

cat(sprintf("Setting 18 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting18_start, units = "secs")))
cat("  Max |theo - emp| PPI++:",
    max(abs(results_unequal$power_theo_ppi - results_unequal$power_emp_ppi)), "\n\n")
