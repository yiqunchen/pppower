#!/usr/bin/env Rscript
# =============================================================================
# Setting 11: Rule-of-Thumb Validation
# =============================================================================
# Validates the approximation n_PPI++ / n_classical ≈ 1 - rho^2 for mean
# estimation across a sweep of correlation values and unlabeled pool sizes.

cat("Setting 11: Rule-of-Thumb Validation\n")
cat("--------------------------------------\n")
setting11_start <- Sys.time()

library(pppower)

rho_grid     <- seq(0.10, 0.99, by = 0.05)
N_rot_values <- c(200, 500, 1000, 5000)
delta_rot    <- 0.3
power_rot    <- 0.80
alpha_rot    <- 0.05
sigma_y2_rot <- 1.0

z_alpha <- qnorm(1 - alpha_rot / 2)
z_beta  <- qnorm(power_rot)
S2 <- (delta_rot / (z_alpha + z_beta))^2
n_classical <- ceiling(sigma_y2_rot / S2)

grid_11 <- expand.grid(
  rho = rho_grid,
  N   = N_rot_values,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_11 <- nrow(grid_11)

results_rule_of_thumb <- vector("list", total_11)

for (i in seq_len(total_11)) {
  rho <- grid_11$rho[i]
  N   <- grid_11$N[i]

  progress_bar(i, total_11, "Setting 11",
               sprintf("rho=%.2f, N=%d: ", rho, N),
               start_time = setting11_start)

  sigma_f2 <- rho^2 * sigma_y2_rot
  cov_yf   <- rho * sqrt(sigma_y2_rot * sigma_f2)

  n_ppi <- tryCatch(
    power_ppi_mean(
      delta = delta_rot, N = N, n = NULL, power = power_rot, alpha = alpha_rot,
      sigma_y2 = sigma_y2_rot, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
      lambda_mode = "oracle"
    ),
    error = function(e) NA_real_
  )

  ratio       <- n_ppi / n_classical
  theoretical <- 1 - rho^2

  results_rule_of_thumb[[i]] <- data.frame(
    rho = rho,
    rho_sq = rho^2,
    N = N,
    n_ppi = n_ppi,
    n_classical = n_classical,
    ratio = ratio,
    theoretical_ratio = theoretical,
    abs_diff = abs(ratio - theoretical),
    stringsAsFactors = FALSE
  )
}

results_rule_of_thumb <- do.call(rbind, results_rule_of_thumb)

cat(sprintf("Setting 11 completed. (%.1f seconds)\n",
            difftime(Sys.time(), setting11_start, units = "secs")))
cat(sprintf("  Max |ratio - (1-rho^2)|: %.4f\n",
            max(results_rule_of_thumb$abs_diff, na.rm = TRUE)))
