# =============================================================================
# Simulation utilities for power validation scripts
# =============================================================================

# Progress tracking helper with elapsed time and ETA
progress_bar <- function(current, total, setting_name, prefix = "",
                         start_time = NULL) {
  pct <- round(100 * current / total)
  bar_width <- 30
  filled <- round(bar_width * current / total)
  bar <- paste0("[", paste(rep("=", filled), collapse = ""),
                paste(rep(" ", bar_width - filled), collapse = ""), "]")

  time_str <- ""
  if (!is.null(start_time) && current > 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    eta     <- elapsed / current * (total - current)
    time_str <- sprintf(" [%s elapsed, ~%s left]",
                        format_seconds(elapsed), format_seconds(eta))
  }

  cat(sprintf("\r%s%s %s %d%% (%d/%d)%s    ",
              prefix, setting_name, bar, pct, current, total, time_str))
  if (current == total) cat("\n")
  flush.console()
}

# Format seconds into human-readable string
format_seconds <- function(secs) {
  if (is.na(secs) || secs < 0) return("?")
  if (secs < 60)   return(sprintf("%.0fs", secs))
  if (secs < 3600) return(sprintf("%.0fm%02.0fs", secs %/% 60, secs %% 60))
  return(sprintf("%.0fh%02.0fm", secs %/% 3600, (secs %% 3600) %/% 60))
}

# =============================================================================
# HELPER WRAPPERS (thin wrappers around package functions for simulation)
# =============================================================================

#' Wrapper for ppi_mean_test that accepts lambda_oracle parameter name
run_ppi_mean_test <- function(Y_L, f_L, f_U, theta0 = 0, alpha = 0.05, lambda_oracle = NULL) {
  result <- pppower::ppi_mean_test(Y_L, f_L, f_U, theta0 = theta0, alpha = alpha,
                                  lambda = lambda_oracle)
  # Return with legacy field names for compatibility
  list(
    theta_hat = result$estimate,
    se = result$se,
    z_stat = result$z_stat,
    p_value = result$p_value,
    reject = result$reject,
    lambda = result$lambda
  )
}

#' Wrapper for ppi_ttest that accepts lambda_A_oracle/lambda_B_oracle parameter names
run_ppi_ttest <- function(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                          delta0 = 0, alpha = 0.05,
                          lambda_A_oracle = NULL, lambda_B_oracle = NULL) {
  result <- pppower::ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                              delta0 = delta0, alpha = alpha,
                              lambda_A = lambda_A_oracle, lambda_B = lambda_B_oracle)
  # Return with legacy field names for compatibility
  list(
    delta_hat = result$estimate,
    se = result$se,
    z_stat = result$z_stat,
    p_value = result$p_value,
    reject = result$reject,
    lambda_A = result$lambda_A,
    lambda_B = result$lambda_B
  )
}

#' Classical one-sample mean test (wraps t.test)
classical_mean_test <- function(Y, theta0 = 0, alpha = 0.05) {
  test <- t.test(Y, mu = theta0, conf.level = 1 - alpha)
  list(
    theta_hat = as.numeric(test$estimate),
    se = test$stderr,
    reject = test$p.value < alpha
  )
}

#' Classical two-sample t-test (wraps t.test)
classical_ttest <- function(Y_A, Y_B, delta0 = 0, alpha = 0.05) {
  test <- t.test(Y_A, Y_B, mu = delta0, conf.level = 1 - alpha, var.equal = FALSE)
  list(
    delta_hat = as.numeric(diff(rev(test$estimate))),
    se = test$stderr,
    reject = test$p.value < alpha
  )
}

# =============================================================================
# THEORETICAL POWER FUNCTIONS
# Wrapper functions that call pppower and pwr packages
# =============================================================================

#' PPI++ power for mean estimation (wraps pppower::power_ppi_mean)
#' @param rho Correlation between Y and f (used to derive cov_y_f)
theo_power_ppi_mean <- function(delta, n, N, sigma_Y2, sigma_f2, cov_y_f, alpha = 0.05) {
  pppower::power_ppi_mean(
    delta    = delta,
    N        = N,
    n        = n,
    alpha    = alpha,
    sigma_y2 = sigma_Y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f
  )
}

#' PPI++ power for paired t-test (wraps pppower::power_ppi_paired)
theo_power_ppi_paired <- function(delta, n, N, sigma_D2, rho_D, alpha = 0.05) {
  pppower::power_ppi_paired(
    delta    = delta,
    N        = N,
    n        = n,
    alpha    = alpha,
    sigma_D2 = sigma_D2,
    rho_D    = rho_D
  )
}

#' PPI++ power for paired binary test (wraps pppower::power_ppi_paired_binary)
theo_power_ppi_paired_binary <- function(delta, n, N, p_A, p_B,
                                         rho_within, sens, spec, alpha = 0.05) {
  pppower::power_ppi_paired_binary(
    delta = delta,
    N = N,
    n = n,
    alpha = alpha,
    p_A = p_A,
    p_B = p_B,
    rho_within = rho_within,
    sens = sens,
    spec = spec
  )
}

#' Vanilla PPI power for mean estimation (lambda = 1, no tuning)
theo_power_vanilla_ppi_mean <- function(delta, n, N, sigma_Y2, sigma_f2, cov_y_f, alpha = 0.05) {
  pppower::power_ppi_mean(
    delta    = delta,
    N        = N,
    n        = n,
    alpha    = alpha,
    sigma_y2 = sigma_Y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f,
    lambda   = 1,
    lambda_type = "user"
  )
}

#' Classical power for one-sample mean (wraps pwr::pwr.t.test)
theo_power_classical_onesample <- function(delta, n, sigma_Y2, alpha = 0.05) {
  d <- abs(delta) / sqrt(sigma_Y2)
  pwr::pwr.t.test(n = n, d = d, sig.level = alpha, type = "one.sample",
                 alternative = "two.sided")$power
}

#' Classical power for two-sample t-test (wraps pwr::pwr.t.test)
theo_power_classical_twosample <- function(delta, n, sigma_Y2, alpha = 0.05) {
  d <- abs(delta) / sqrt(sigma_Y2)
  pwr::pwr.t.test(n = n, d = d, sig.level = alpha, type = "two.sample",
                 alternative = "two.sided")$power
}

#' PPI++ power for two-sample test (wraps pppower::power_ppi_ttest)
theo_power_ppi_twosample <- function(delta, n_A, n_B, N_A, N_B,
                                     sigma_Y2_A, sigma_f2_A, cov_A,
                                     sigma_Y2_B, sigma_f2_B, cov_B,
                                     alpha = 0.05) {
  pppower::power_ppi_ttest(
    delta = delta,
    n_A = n_A, n_B = n_B,
    N_A = N_A, N_B = N_B,
    alpha = alpha,
    sigma_y2_A = sigma_Y2_A, sigma_f2_A = sigma_f2_A, cov_yf_A = cov_A,
    sigma_y2_B = sigma_Y2_B, sigma_f2_B = sigma_f2_B, cov_yf_B = cov_B
  )
}

#' Monte Carlo power for two-sample binary proportion tests (PPI++)
#'
#' Uses sensitivity/specificity to derive population moments, computes
#' theoretical power via `power_ppi_ttest_binary()`, and estimates
#' empirical power via repeated simulations with oracle lambdas.
simulate_twosample_binary_power <- function(p0, delta, sens, spec,
                                            n_values, N_values,
                                            R = 1000,
                                            alpha = 0.05,
                                            seed = NULL,
                                            progress = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  p_A <- p0 + delta / 2
  p_B <- p0 - delta / 2

  moments_A <- pppower::binary_moments_from_sens_spec(p = p_A, sens = sens, spec = spec)
  moments_B <- pppower::binary_moments_from_sens_spec(p = p_B, sens = sens, spec = spec)

  grid <- expand.grid(N = N_values, n = n_values, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  total <- nrow(grid)

  out <- vector("list", total)
  for (i in seq_len(total)) {
    N <- grid$N[i]; n <- grid$n[i]
    if (progress) {
      progress_bar(i, total, "Setting 4",
                   sprintf("sens=%.0f%%, N=%d, n=%d: ", sens * 100, N, n))
    }

    lambda_A <- moments_A$cov_y_f / ((1 + n / N) * moments_A$sigma_f2)
    lambda_B <- moments_B$cov_y_f / ((1 + n / N) * moments_B$sigma_f2)

    # theoretical power
    power_theo_ppi <- pppower::power_ppi_ttest_binary(
      p_A = p_A, p_B = p_B,
      n_A = n, n_B = n,
      N_A = N, N_B = N,
      sens_A = sens, spec_A = spec,
      sens_B = sens, spec_B = spec,
      delta = delta,
      alpha = alpha
    )

    # classical (use average variance)
    sigma_Y2_avg <- (moments_A$sigma_y2 + moments_B$sigma_y2) / 2
    power_theo_classical <- theo_power_classical_twosample(delta, n, sigma_Y2_avg, alpha)

    # Monte Carlo
    rej_ppi <- 0
    rej_classical <- 0

    for (r in seq_len(R)) {
      Y_A <- rbinom(n, 1, p_A)
      f_A_L <- ifelse(Y_A == 1, rbinom(n, 1, sens), rbinom(n, 1, 1 - spec))
      Y_A_U <- rbinom(N, 1, p_A)
      f_A_U <- ifelse(Y_A_U == 1, rbinom(N, 1, sens), rbinom(N, 1, 1 - spec))

      Y_B <- rbinom(n, 1, p_B)
      f_B_L <- ifelse(Y_B == 1, rbinom(n, 1, sens), rbinom(n, 1, 1 - spec))
      Y_B_U <- rbinom(N, 1, p_B)
      f_B_U <- ifelse(Y_B_U == 1, rbinom(N, 1, sens), rbinom(N, 1, 1 - spec))

      test_ppi <- pppower::ppi_ttest(
        Y_A, f_A_L, f_A_U,
        Y_B, f_B_L, f_B_U,
        delta0 = 0,
        alpha = alpha,
        lambda_A = lambda_A,
        lambda_B = lambda_B
      )
      rej_ppi <- rej_ppi + as.integer(test_ppi$reject)

      test_classical <- classical_ttest(Y_A, Y_B, delta0 = 0, alpha = alpha)
      rej_classical <- rej_classical + as.integer(test_classical$reject)
    }

    out[[i]] <- data.frame(
      setting = "Proportion Test (Binary)",
      sens = sens,
      spec = spec,
      n = n,
      N = N,
      delta = delta,
      rho_A = moments_A$cov_y_f / sqrt(moments_A$sigma_y2 * moments_A$sigma_f2),
      rho_B = moments_B$cov_y_f / sqrt(moments_B$sigma_y2 * moments_B$sigma_f2),
      power_theo_ppi = power_theo_ppi,
      power_emp_ppi = rej_ppi / R,
      power_theo_classical = power_theo_classical,
      power_emp_classical = rej_classical / R
    )
  }

  do.call(rbind, out)
}
