#!/usr/bin/env Rscript
# =============================================================================
# Head-to-head comparison table: Classical vs Vanilla PPI vs PPI++ (oracle/plugin)
# =============================================================================

library(pppower)
source(file.path("simulation_studies", "sim_utils.R"))

set.seed(20260303)

# Representative planning scenario (continuous one-sample mean)
alpha <- 0.05
target_power <- 0.80
R <- 2000

N <- 5000
delta <- 0.2
sigma_y2 <- 1
sigma_f2 <- 1
rho <- 0.7
cov_yf <- rho * sqrt(sigma_y2 * sigma_f2)

z_alpha <- qnorm(1 - alpha / 2)
z_beta <- qnorm(target_power)
S2 <- (delta / (z_alpha + z_beta))^2
n_classical <- ceiling(sigma_y2 / S2)

n_vanilla <- ceiling(power_ppi_mean(
  delta = delta, N = N, n = NULL, power = target_power, alpha = alpha,
  sigma_y2 = sigma_y2, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
  lambda_mode = "vanilla"
))

n_oracle <- ceiling(power_ppi_mean(
  delta = delta, N = N, n = NULL, power = target_power, alpha = alpha,
  sigma_y2 = sigma_y2, sigma_f2 = sigma_f2, cov_y_f = cov_yf,
  lambda_mode = "oracle"
))

# Practical plugin planning: use oracle closed-form n and evaluate achieved performance
n_plugin <- n_oracle

simulate_rejection_rate <- function(n, effect, method) {
  Sigma <- matrix(c(sigma_y2, cov_yf, cov_yf, sigma_f2), 2, 2)
  L <- chol(Sigma)
  lambda_oracle <- cov_yf / ((1 + n / N) * sigma_f2)

  rej <- 0
  for (r in seq_len(R)) {
    Z_L <- matrix(rnorm(n * 2), n, 2) %*% L
    Y_L <- Z_L[, 1] + effect
    f_L <- Z_L[, 2] + effect

    Z_U <- matrix(rnorm(N * 2), N, 2) %*% L
    f_U <- Z_U[, 2] + effect

    if (method == "Classical") {
      test <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha)
    } else if (method == "Vanilla PPI") {
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha,
                                lambda_oracle = 1)
    } else if (method == "PPI++ (Oracle)") {
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha,
                                lambda_oracle = lambda_oracle)
    } else if (method == "PPI++ (Plugin)") {
      cov_hat <- cov(Y_L, f_L)
      varf_hat <- var(f_L)
      lambda_hat <- if (is.finite(varf_hat) && varf_hat > 0) {
        cov_hat / ((1 + n / N) * varf_hat)
      } else {
        0
      }
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha,
                                lambda_oracle = lambda_hat)
    } else {
      stop("Unknown method")
    }

    rej <- rej + as.integer(test$reject)
  }

  rej / R
}

methods <- data.frame(
  method = c("Classical", "Vanilla PPI", "PPI++ (Oracle)", "PPI++ (Plugin)"),
  n_required = c(n_classical, n_vanilla, n_oracle, n_plugin),
  stringsAsFactors = FALSE
)

methods$achieved_power <- NA_real_
methods$type1_error <- NA_real_

for (i in seq_len(nrow(methods))) {
  m <- methods$method[i]
  n_i <- methods$n_required[i]
  methods$achieved_power[i] <- simulate_rejection_rate(n_i, effect = delta, method = m)
  methods$type1_error[i] <- simulate_rejection_rate(n_i, effect = 0, method = m)
}

methods$reduction_vs_classical <- round(100 * (1 - methods$n_required / n_classical), 1)
methods$achieved_power <- round(methods$achieved_power, 3)
methods$type1_error <- round(methods$type1_error, 3)

out_dir <- file.path("simulation_studies", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

csv_path <- file.path(out_dir, "head_to_head_one_sample.csv")
write.csv(methods, csv_path, row.names = FALSE)

cat("Saved:\n")
cat("  ", csv_path, "\n")
print(methods)
