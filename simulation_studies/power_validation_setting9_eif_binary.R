#!/usr/bin/env Rscript
# =============================================================================
# Setting 9: EIF Binary Validation
# =============================================================================
# Validates the analytical power formula for the EIF estimator with binary
# surrogates (e.g., LLM-as-a-judge) against Monte Carlo simulations.

cat("Setting 9: EIF Binary Validation\n")
cat("----------------------------------\n")
setting9_start <- Sys.time()

library(pppower)

# Parameters
p_values     <- c(0.3, 0.5)
sens_spec    <- list(c(0.70, 0.70), c(0.85, 0.85), c(0.95, 0.95))
N_values     <- c(200, 500)
m_cal_values <- c(20, 40, 60, 80, 100)
delta_eif    <- 0.05
R_eif        <- 1000
alpha_eif    <- 0.05

grid_9 <- expand.grid(
  p     = p_values,
  ss_ix = seq_along(sens_spec),
  N     = N_values,
  m_cal = m_cal_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_9 <- nrow(grid_9)

results_eif_binary <- vector("list", total_9)

for (i in seq_len(total_9)) {
  row <- grid_9[i, ]
  ss  <- sens_spec[[row$ss_ix]]
  sens_i <- ss[1]; spec_i <- ss[2]

  progress_bar(i, total_9, "Setting 9",
               sprintf("p=%.1f, sens=%.2f, N=%d, m=%d: ",
                       row$p, sens_i, row$N, row$m_cal),
               start_time = setting9_start)

  # Analytical power (PPI EIF)
  power_theo <- tryCatch(
    power_eif_binary(
      delta = delta_eif, N = row$N, m_cal = row$m_cal,
      alpha = alpha_eif, p = row$p, sens = sens_i, spec = spec_i
    ),
    error = function(e) NA_real_
  )

  # Empirical power via Monte Carlo
  sim_res <- tryCatch(
    simulate_eif_binary(
      R = R_eif, delta = delta_eif, N = row$N, m_cal = row$m_cal,
      alpha = alpha_eif, p = row$p, sens = sens_i, spec = spec_i,
      seed = 9000 + i
    ),
    error = function(e) list(empirical_power = NA_real_)
  )
  power_emp <- sim_res$empirical_power

  # Classical power (binomial proportion test, no auxiliary data)
  se_classical <- sqrt(row$p * (1 - row$p) / row$m_cal)
  z_alpha <- qnorm(1 - alpha_eif / 2)
  if (se_classical > 0) {
    mu_cl <- abs(delta_eif) / se_classical
    power_classical <- 1 - pnorm(z_alpha - mu_cl) + pnorm(-z_alpha - mu_cl)
  } else {
    power_classical <- NA_real_
  }

  results_eif_binary[[i]] <- data.frame(
    setting = "EIF Binary",
    p = row$p,
    sens = sens_i,
    spec = spec_i,
    N = row$N,
    m_cal = row$m_cal,
    delta = delta_eif,
    power_theo = as.numeric(power_theo),
    power_emp = power_emp,
    power_classical = power_classical,
    stringsAsFactors = FALSE
  )
}

results_eif_binary <- do.call(rbind, results_eif_binary)

cat(sprintf("Setting 9 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting9_start, units = "secs")))
cat(sprintf("  Max |theo - emp| for EIF: %.4f\n",
            max(abs(results_eif_binary$power_theo - results_eif_binary$power_emp),
                na.rm = TRUE)))
