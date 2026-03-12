#!/usr/bin/env Rscript
# =============================================================================
# Setting 19: Misspecified Rho (Planning vs True)
# =============================================================================
# Practical robustness: the user plans a study using rho_planning (from
# a pilot or published report) but the true rho differs.  We compute:
# (a) The required n at rho_planning  ->  n_planned
# (b) The actual power at n_planned under rho_true (analytical + MC)
#
# This quantifies how sensitive the sample size formula is to
# misspecification of the prediction quality parameter.

cat("Setting 19: Misspecified Rho (Planning vs True)\n")
cat("------------------------------------------------\n")
setting19_start <- Sys.time()

library(pppower)

# Parameters
rho_planning_values <- c(0.5, 0.7, 0.9)
rho_true_offsets <- c(-0.20, -0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15, 0.20)
N_19 <- 500
sigma_Y2_19 <- 1
sigma_f2_19 <- 1
delta_19 <- 0.2
target_power <- 0.80
alpha_19 <- 0.05
R_19 <- 2000  # more reps for robustness setting

# Build grid: rho_planning x rho_true_offset, filtering infeasible combos
grid_rows <- list()
for (rho_plan in rho_planning_values) {
  for (offset in rho_true_offsets) {
    rho_true <- rho_plan + offset
    if (rho_true > 0.01 && rho_true < 0.99) {
      grid_rows[[length(grid_rows) + 1]] <- data.frame(
        rho_planning = rho_plan,
        rho_true = rho_true,
        offset = offset,
        stringsAsFactors = FALSE
      )
    }
  }
}
grid_19 <- do.call(rbind, grid_rows)
total_19 <- nrow(grid_19)

results_misspecified_rho <- vector("list", total_19)

for (i in seq_len(total_19)) {
  row <- grid_19[i, ]
  progress_bar(i, total_19, "Setting 19",
               sprintf("plan=%.2f, true=%.2f: ", row$rho_planning, row$rho_true),
               start_time = setting19_start)

  # Step 1: Compute n_planned using rho_planning
  cov_Yf_plan <- row$rho_planning * sqrt(sigma_Y2_19 * sigma_f2_19)
  n_planned <- power_ppi_mean(
    delta = delta_19, N = N_19, power = target_power,
    sigma_y2 = sigma_Y2_19, sigma_f2 = sigma_f2_19,
    cov_y_f = cov_Yf_plan, alpha = alpha_19
  )
  n_planned <- ceiling(n_planned)

  # Step 2: Compute analytical power at n_planned under rho_true
  cov_Yf_true <- row$rho_true * sqrt(sigma_Y2_19 * sigma_f2_19)
  power_theo_true <- theo_power_ppi_mean(
    delta_19, n_planned, N_19,
    sigma_Y2_19, sigma_f2_19, cov_Yf_true, alpha_19
  )

  # Also compute classical n for reference
  n_classical <- ceiling(sigma_Y2_19 / (delta_19 / (qnorm(1 - alpha_19/2) + qnorm(target_power)))^2)
  power_classical <- theo_power_classical_onesample(delta_19, n_planned, sigma_Y2_19, alpha_19)

  # Step 3: Monte Carlo power at n_planned under rho_true
  Sigma_true <- matrix(c(sigma_Y2_19, cov_Yf_true, cov_Yf_true, sigma_f2_19), 2, 2)
  L_true <- chol(Sigma_true)
  lambda_oracle_true <- cov_Yf_true / ((1 + n_planned / N_19) * sigma_f2_19)

  rej_ppi <- 0
  for (r in seq_len(R_19)) {
    Z_L <- matrix(rnorm(n_planned * 2), n_planned, 2) %*% L_true
    Y_L <- Z_L[, 1] + delta_19
    f_L <- Z_L[, 2] + delta_19

    Z_U <- matrix(rnorm(N_19 * 2), N_19, 2) %*% L_true
    f_U <- Z_U[, 2] + delta_19

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha_19,
                                  lambda_oracle = lambda_oracle_true)
    rej_ppi <- rej_ppi + test_ppi$reject
  }

  results_misspecified_rho[[i]] <- data.frame(
    setting = "Misspecified Rho",
    rho_planning = row$rho_planning,
    rho_true = row$rho_true,
    offset = row$offset,
    n_planned = n_planned,
    n_classical = n_classical,
    power_target = target_power,
    power_theo_true = power_theo_true,
    power_emp_true = rej_ppi / R_19,
    power_classical = power_classical,
    N = N_19,
    delta = delta_19,
    R = R_19,
    stringsAsFactors = FALSE
  )
}

results_setting19 <- do.call(rbind, results_misspecified_rho)

# Summary
cat("\n\nSetting 19 Summary:\n")
cat("===================\n")
for (rp in rho_planning_values) {
  sub <- results_setting19[results_setting19$rho_planning == rp, ]
  cat(sprintf("\nrho_planning = %.2f (n_planned = %d):\n", rp, sub$n_planned[1]))
  for (j in seq_len(nrow(sub))) {
    r <- sub[j, ]
    flag <- if (r$power_emp_true < target_power - 0.05) " ** UNDERPOWERED" else ""
    cat(sprintf("  rho_true = %.2f: analytical = %.3f, MC = %.3f%s\n",
                r$rho_true, r$power_theo_true, r$power_emp_true, flag))
  }
}

elapsed_19 <- difftime(Sys.time(), setting19_start, units = "secs")
cat(sprintf("\nSetting 19 completed in %s\n", format_seconds(as.numeric(elapsed_19))))
