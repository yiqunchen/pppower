#!/usr/bin/env Rscript
# =============================================================================
# Setting 10: Sample Size Inversion Validation
# =============================================================================
# Verifies that n_required_*() functions return sample sizes that achieve
# (or exceed) the target power, both analytically and via Monte Carlo.

cat("Setting 10: Sample Size Inversion Validation\n")
cat("-----------------------------------------------\n")
setting10_start <- Sys.time()

library(pppower)

target_powers <- c(0.60, 0.70, 0.80, 0.90)
alpha_10      <- 0.05
R_10          <- 2000

# ---- Mean (Continuous) ----
rho_values_10 <- c(0.5, 0.7, 0.9)
N_10          <- 5000
delta_10      <- 0.2
sigma_y2_10   <- 1.0

results_10 <- list()
counter <- 0
total_10 <- length(target_powers) * length(rho_values_10) + length(target_powers) * 2

for (rho in rho_values_10) {
  sigma_f2 <- rho^2 * sigma_y2_10
  cov_yf   <- rho * sqrt(sigma_y2_10 * sigma_f2)

  for (tp in target_powers) {
    counter <- counter + 1
    progress_bar(counter, total_10, "Setting 10",
                 sprintf("mean-cont rho=%.1f, target=%.2f: ", rho, tp))

    n_star <- tryCatch(
      n_required_ppi_pp(
        delta = delta_10, N = N_10, power = tp, alpha = alpha_10,
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

    # Analytical check
    analytical_power <- power_ppi_pp_mean(
      delta = delta_10, N = N_10, n = n_star, alpha = alpha_10,
      sigma_y2 = sigma_y2_10, sigma_f2 = sigma_f2, cov_y_f = cov_yf
    )

    # MC check
    sim <- simulate_power_ppiplus_mean(
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

# ---- EIF Binary ----
p_10     <- 0.5
sens_10  <- 0.85
spec_10  <- 0.85
N_eif_10 <- 5000
delta_eif_10 <- 0.05

for (tp in target_powers) {
  counter <- counter + 1
  progress_bar(counter, total_10, "Setting 10",
               sprintf("eif-binary target=%.2f: ", tp))

  m_star <- tryCatch(
    n_required_eif_binary(
      delta = delta_eif_10, N = N_eif_10, power = tp, alpha = alpha_10,
      p = p_10, sens = sens_10, spec = spec_10
    ),
    error = function(e) NA_integer_
  )

  if (is.na(m_star)) {
    results_10[[counter]] <- data.frame(
      design = "eif_binary", rho = NA, target_power = tp,
      n_star = NA, analytical_achieved = NA, empirical_achieved = NA,
      stringsAsFactors = FALSE
    )
    next
  }

  analytical_power <- power_eif_binary(
    delta = delta_eif_10, N = N_eif_10, m_cal = m_star, alpha = alpha_10,
    p = p_10, sens = sens_10, spec = spec_10
  )

  sim <- simulate_power_eif_binary(
    R = R_10, delta = delta_eif_10, N = N_eif_10, m_cal = m_star,
    alpha = alpha_10, p = p_10, sens = sens_10, spec = spec_10,
    seed = 10100 + counter
  )

  results_10[[counter]] <- data.frame(
    design = "eif_binary", rho = NA, target_power = tp,
    n_star = m_star,
    analytical_achieved = as.numeric(analytical_power),
    empirical_achieved = sim$empirical_power,
    stringsAsFactors = FALSE
  )
}

# ---- Paired (Continuous) ----
rho_D_10  <- 0.7
N_pair_10 <- 5000
delta_pair_10 <- 0.3

for (tp in target_powers) {
  counter <- counter + 1
  progress_bar(counter, total_10, "Setting 10",
               sprintf("paired-cont target=%.2f: ", tp))

  n_star <- tryCatch(
    n_required_ppi_pp_paired(
      delta = delta_pair_10, N = N_pair_10, power = tp, alpha = alpha_10,
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

  analytical_power <- power_ppi_pp_paired(
    delta = delta_pair_10, N = N_pair_10, n = n_star, alpha = alpha_10,
    sigma_D2 = 1.0, rho_D = rho_D_10
  )

  results_10[[counter]] <- data.frame(
    design = "paired_continuous", rho = rho_D_10, target_power = tp,
    n_star = n_star,
    analytical_achieved = analytical_power,
    empirical_achieved = NA_real_,
    stringsAsFactors = FALSE
  )
}

results_n_inversion <- do.call(rbind, results_10)

cat(sprintf("Setting 10 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting10_start, units = "secs")))
cat("  Sample size inversion results:\n")
print(results_n_inversion)
