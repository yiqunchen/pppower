#!/usr/bin/env Rscript
# =============================================================================
# Setting 10: Sample Size Inversion Validation
# =============================================================================
# Verifies that n_required_*() functions return sample sizes that achieve
# (or exceed) the target power, both analytically and via Monte Carlo.
#
# Covers: mean (continuous), mean (binary), t-test (via mean per-group),
#         paired (continuous).

cat("Setting 10: Sample Size Inversion Validation\n")
cat("-----------------------------------------------\n")
setting10_start <- Sys.time()

library(pppower)

target_powers <- c(0.60, 0.70, 0.80, 0.90)
alpha_10      <- 0.05
R_10          <- 2000

results_10 <- list()
counter    <- 0

# Count total configs for progress bar
total_10 <- length(target_powers) * (
  3 +   # mean_continuous: 3 rho values
  3 +   # mean_binary: 3 sens/spec combos
  3 +   # ttest_continuous: 3 rho values
  1     # paired_continuous: 1 rho_D
)

# =============================================================================
# Mean (Continuous)
# =============================================================================
rho_values_10 <- c(0.5, 0.7, 0.9)
N_10          <- 500
delta_10      <- 0.2
sigma_y2_10   <- 1.0

for (rho in rho_values_10) {
  sigma_f2 <- rho^2 * sigma_y2_10
  cov_yf   <- rho * sqrt(sigma_y2_10 * sigma_f2)

  for (tp in target_powers) {
    counter <- counter + 1
    progress_bar(counter, total_10, "Setting 10",
                 sprintf("mean-cont rho=%.1f, target=%.2f: ", rho, tp),
                 start_time = setting10_start)

    n_star <- tryCatch(
      power_ppi_mean(
        delta = delta_10, N = N_10, n = NULL, power = tp, alpha = alpha_10,
        sigma_y2 = sigma_y2_10, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
        lambda_mode = "oracle"
      ),
      error = function(e) NA_integer_
    )

    if (is.na(n_star)) {
      results_10[[counter]] <- data.frame(
        design = "mean_continuous", rho = rho, target_power = tp,
        n_star = NA, analytical_achieved = NA, empirical_achieved = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    analytical_power <- power_ppi_mean(
      delta = delta_10, N = N_10, n = n_star, alpha = alpha_10,
      sigma_y2 = sigma_y2_10, sigma_f2 = sigma_f2, cov_y_f = cov_yf
    )

    sim <- simulate_ppi_mean(
      R = R_10, n = n_star, N = N_10, delta = delta_10,
      var_f = sigma_f2, var_res = sigma_y2_10 - 2 * cov_yf + sigma_f2,
      alpha = alpha_10, seed = 10000 + counter
    )

    results_10[[counter]] <- data.frame(
      design = "mean_continuous", rho = rho, target_power = tp,
      n_star = n_star,
      analytical_achieved = analytical_power,
      empirical_achieved = sim$empirical_power,
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# Mean (Binary)
# =============================================================================
p0_10_bin    <- 0.3
delta_10_bin <- 0.05
N_10_bin     <- 500
sens_spec_10 <- list(c(0.70, 0.70), c(0.85, 0.85), c(0.95, 0.95))

for (ss in sens_spec_10) {
  sens_i <- ss[1]; spec_i <- ss[2]

  for (tp in target_powers) {
    counter <- counter + 1
    progress_bar(counter, total_10, "Setting 10",
                 sprintf("mean-bin sens=%.2f, target=%.2f: ", sens_i, tp),
                 start_time = setting10_start)

    n_star <- tryCatch(
      power_ppi_mean(
        delta = delta_10_bin, N = N_10_bin, n = NULL, power = tp, alpha = alpha_10,
        metrics = list(
          sensitivity = sens_i, specificity = spec_i,
          p_y = p0_10_bin, m_obs = 500
        ),
        metric_type = "classification",
        lambda_mode = "oracle"
      ),
      error = function(e) NA_integer_
    )

    if (is.na(n_star)) {
      results_10[[counter]] <- data.frame(
        design = "mean_binary", rho = NA, target_power = tp,
        n_star = NA, analytical_achieved = NA, empirical_achieved = NA,
        stringsAsFactors = FALSE
      )
      next
    }

    analytical_power <- power_ppi_mean(
      delta = delta_10_bin, N = N_10_bin, n = n_star, alpha = alpha_10,
      metrics = list(
        sensitivity = sens_i, specificity = spec_i,
        p_y = p0_10_bin, m_obs = n_star
      ),
      metric_type = "classification"
    )

    # MC check via direct simulation
    moments_bin <- binary_moments_from_sens_spec(p = p0_10_bin + delta_10_bin,
                                                  sens = sens_i, spec = spec_i)
    lambda_oracle_bin <- moments_bin$cov_y_f /
      ((1 + n_star / N_10_bin) * moments_bin$sigma_f2)

    rej_mc <- 0
    for (r in seq_len(R_10)) {
      set.seed(10200 + counter * R_10 + r)
      Y_L <- rbinom(n_star, 1, p0_10_bin + delta_10_bin)
      f_L <- ifelse(Y_L == 1, rbinom(n_star, 1, sens_i),
                     rbinom(n_star, 1, 1 - spec_i))
      Y_U_h <- rbinom(N_10_bin, 1, p0_10_bin + delta_10_bin)
      f_U <- ifelse(Y_U_h == 1, rbinom(N_10_bin, 1, sens_i),
                     rbinom(N_10_bin, 1, 1 - spec_i))
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = p0_10_bin,
                                alpha = alpha_10,
                                lambda_oracle = lambda_oracle_bin)
      rej_mc <- rej_mc + test$reject
    }

    results_10[[counter]] <- data.frame(
      design = "mean_binary", rho = NA, target_power = tp,
      n_star = n_star,
      analytical_achieved = analytical_power,
      empirical_achieved = rej_mc / R_10,
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# Two-sample t-test (Continuous) — solve via per-group n
# =============================================================================
delta_tt_10  <- 0.3
N_tt_10      <- 500
sigma_Y2_tt  <- 1.0
sigma_f2_tt  <- 1.0

for (rho in rho_values_10) {
  cov_tt <- rho * sqrt(sigma_Y2_tt * sigma_f2_tt)

  for (tp in target_powers) {
    counter <- counter + 1
    progress_bar(counter, total_10, "Setting 10",
                 sprintf("ttest-cont rho=%.1f, target=%.2f: ", rho, tp),
                 start_time = setting10_start)

    # Approximate n per group: 2 * sigma_Y2 * (1 - rho^2) / S^2
    z_alpha <- qnorm(1 - alpha_10 / 2)
    z_beta  <- qnorm(tp)
    S2_tt   <- (delta_tt_10 / (z_alpha + z_beta))^2
    n_star  <- ceiling(2 * sigma_Y2_tt * (1 - rho^2) / S2_tt)
    n_star  <- max(n_star, 10L)

    # Analytical check
    analytical_power <- tryCatch(
      power_ppi_ttest(
        delta = delta_tt_10,
        n_A = n_star, n_B = n_star,
        N_A = N_tt_10, N_B = N_tt_10,
        alpha = alpha_10,
        sigma_y2_A = sigma_Y2_tt, sigma_f2_A = sigma_f2_tt, cov_yf_A = cov_tt,
        sigma_y2_B = sigma_Y2_tt, sigma_f2_B = sigma_f2_tt, cov_yf_B = cov_tt
      ),
      error = function(e) NA_real_
    )

    # MC check
    Sigma_tt <- matrix(c(sigma_Y2_tt, cov_tt, cov_tt, sigma_f2_tt), 2, 2)
    L_tt     <- chol(Sigma_tt)
    lambda_oracle_tt <- cov_tt / ((1 + n_star / N_tt_10) * sigma_f2_tt)

    rej_mc <- 0
    for (r in seq_len(R_10)) {
      set.seed(10300 + counter * R_10 + r)
      Z_A <- matrix(rnorm(n_star * 2), n_star, 2) %*% L_tt
      Y_A <- Z_A[, 1] + delta_tt_10 / 2
      f_A_L <- Z_A[, 2] + delta_tt_10 / 2
      Z_AU <- matrix(rnorm(N_tt_10 * 2), N_tt_10, 2) %*% L_tt
      f_A_U <- Z_AU[, 2] + delta_tt_10 / 2

      Z_B <- matrix(rnorm(n_star * 2), n_star, 2) %*% L_tt
      Y_B <- Z_B[, 1] - delta_tt_10 / 2
      f_B_L <- Z_B[, 2] - delta_tt_10 / 2
      Z_BU <- matrix(rnorm(N_tt_10 * 2), N_tt_10, 2) %*% L_tt
      f_B_U <- Z_BU[, 2] - delta_tt_10 / 2

      test <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                            delta0 = 0, alpha = alpha_10,
                            lambda_A_oracle = lambda_oracle_tt,
                            lambda_B_oracle = lambda_oracle_tt)
      rej_mc <- rej_mc + test$reject
    }

    results_10[[counter]] <- data.frame(
      design = "ttest_continuous", rho = rho, target_power = tp,
      n_star = n_star,
      analytical_achieved = analytical_power,
      empirical_achieved = rej_mc / R_10,
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# Paired (Continuous)
# =============================================================================
rho_D_10      <- 0.7
N_pair_10     <- 500
delta_pair_10 <- 0.3

for (tp in target_powers) {
  counter <- counter + 1
  progress_bar(counter, total_10, "Setting 10",
               sprintf("paired-cont target=%.2f: ", tp),
               start_time = setting10_start)

  n_star <- tryCatch(
    power_ppi_paired(
      delta = delta_pair_10, N = N_pair_10, n = NULL, power = tp, alpha = alpha_10,
      sigma_D2 = 1.0, rho_D = rho_D_10
    ),
    error = function(e) NA_integer_
  )

  if (is.na(n_star)) {
    results_10[[counter]] <- data.frame(
      design = "paired_continuous", rho = rho_D_10, target_power = tp,
      n_star = NA, analytical_achieved = NA, empirical_achieved = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  analytical_power <- power_ppi_paired(
    delta = delta_pair_10, N = N_pair_10, n = n_star, alpha = alpha_10,
    sigma_D2 = 1.0, rho_D = rho_D_10
  )

  # MC check: simulate paired differences
  sigma_D2_p  <- 1.0
  sigma_fD2_p <- 1.0
  cov_DfD     <- rho_D_10 * sqrt(sigma_D2_p * sigma_fD2_p)
  Sigma_p     <- matrix(c(sigma_D2_p, cov_DfD, cov_DfD, sigma_fD2_p), 2, 2)
  L_p         <- chol(Sigma_p)
  lambda_p    <- cov_DfD / ((1 + n_star / N_pair_10) * sigma_fD2_p)

  rej_mc <- 0
  for (r in seq_len(R_10)) {
    set.seed(10400 + counter * R_10 + r)
    Z_L <- matrix(rnorm(n_star * 2), n_star, 2) %*% L_p
    D_L   <- Z_L[, 1] + delta_pair_10
    fD_L  <- Z_L[, 2] + delta_pair_10
    Z_U <- matrix(rnorm(N_pair_10 * 2), N_pair_10, 2) %*% L_p
    fD_U  <- Z_U[, 2] + delta_pair_10

    test <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0,
                              alpha = alpha_10, lambda_oracle = lambda_p)
    rej_mc <- rej_mc + test$reject
  }

  results_10[[counter]] <- data.frame(
    design = "paired_continuous", rho = rho_D_10, target_power = tp,
    n_star = n_star,
    analytical_achieved = analytical_power,
    empirical_achieved = rej_mc / R_10,
    stringsAsFactors = FALSE
  )
}

results_n_inversion <- do.call(rbind, results_10)

cat(sprintf("Setting 10 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting10_start, units = "secs")))
cat("  Sample size inversion results:\n")
print(results_n_inversion)
