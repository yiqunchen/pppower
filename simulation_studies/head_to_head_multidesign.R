#!/usr/bin/env Rscript
# =============================================================================
# Head-to-head comparison across designs
# Classical vs Vanilla PPI vs PPI++ (oracle/plugin)
# =============================================================================

library(pppower)
source(file.path("simulation_studies", "sim_utils.R"))

set.seed(20260303)

alpha <- 0.05
target_power <- 0.80
R <- 5000

power_two_sided <- function(delta, se, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
}

find_min_n <- function(power_fn, target, lo = 5L, hi = 5000L) {
  while (power_fn(hi) < target && hi < 200000L) hi <- hi * 2L
  while (lo < hi) {
    mid <- (lo + hi) %/% 2L
    if (power_fn(mid) >= target) hi <- mid else lo <- mid + 1L
  }
  lo
}

# ===================== One-sample mean =====================
N1 <- 5000
delta1 <- 0.2
sigma_y2_1 <- 1
sigma_f2_1 <- 1
rho1 <- 0.7
cov1 <- rho1 * sqrt(sigma_y2_1 * sigma_f2_1)

z_alpha <- qnorm(1 - alpha / 2)
z_beta <- qnorm(target_power)
S2 <- (delta1 / (z_alpha + z_beta))^2

n1_classical <- ceiling(sigma_y2_1 / S2)
n1_vanilla <- ceiling(power_ppi_mean(delta = delta1, N = N1, n = NULL, power = target_power,
                                     alpha = alpha, sigma_y2 = sigma_y2_1, sigma_f2 = sigma_f2_1,
                                     cov_y_f = cov1, lambda_mode = "vanilla"))
n1_oracle <- ceiling(power_ppi_mean(delta = delta1, N = N1, n = NULL, power = target_power,
                                    alpha = alpha, sigma_y2 = sigma_y2_1, sigma_f2 = sigma_f2_1,
                                    cov_y_f = cov1, lambda_mode = "oracle"))
n1_plugin <- n1_oracle

sim_rate_one <- function(n, effect, method) {
  Sigma <- matrix(c(sigma_y2_1, cov1, cov1, sigma_f2_1), 2, 2)
  L <- chol(Sigma)
  lambda_oracle <- cov1 / ((1 + n / N1) * sigma_f2_1)
  rej <- 0
  for (r in seq_len(R)) {
    Z_L <- matrix(rnorm(n * 2), n, 2) %*% L
    Y_L <- Z_L[, 1] + effect
    f_L <- Z_L[, 2] + effect
    Z_U <- matrix(rnorm(N1 * 2), N1, 2) %*% L
    f_U <- Z_U[, 2] + effect
    if (method == "Classical") {
      test <- classical_mean_test(Y_L, theta0 = 0, alpha = alpha)
    } else if (method == "Vanilla PPI") {
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha, lambda_oracle = 1)
    } else if (method == "PPI++ (Oracle)") {
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha, lambda_oracle = lambda_oracle)
    } else {
      cov_hat <- cov(Y_L, f_L); varf_hat <- var(f_L)
      lambda_hat <- if (is.finite(varf_hat) && varf_hat > 0) cov_hat / ((1 + n / N1) * varf_hat) else 0
      test <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = alpha, lambda_oracle = lambda_hat)
    }
    rej <- rej + as.integer(test$reject)
  }
  rej / R
}

# ===================== Two-sample =====================
N2 <- 5000
delta2 <- 0.3
sigma_y2_2 <- 1
sigma_f2_2 <- 1
rho2 <- 0.7
cov2 <- rho2 * sqrt(sigma_y2_2 * sigma_f2_2)

# per-group n
n2_classical <- ceiling(2 * sigma_y2_2 / ((delta2 / (z_alpha + z_beta))^2))

power_ttest_with_lambda <- function(delta, n, N, sigma_y2, sigma_f2, cov_yf, lambda, alpha = 0.05) {
  var_g <- sigma_y2 / n + lambda^2 * sigma_f2 * (1 / N + 1 / n) - 2 * lambda * cov_yf / n
  se <- sqrt(2 * max(var_g, 0))
  power_two_sided(delta, se, alpha)
}

lambda2_oracle_fn <- function(n) cov2 / ((1 + n / N2) * sigma_f2_2)

n2_vanilla <- find_min_n(function(n) power_ttest_with_lambda(delta2, n, N2, sigma_y2_2, sigma_f2_2, cov2, 1, alpha), target_power)
n2_oracle <- find_min_n(function(n) {
  l <- lambda2_oracle_fn(n)
  power_ttest_with_lambda(delta2, n, N2, sigma_y2_2, sigma_f2_2, cov2, l, alpha)
}, target_power)
n2_plugin <- n2_oracle

sim_rate_ttest <- function(n, effect, method) {
  Sigma <- matrix(c(sigma_y2_2, cov2, cov2, sigma_f2_2), 2, 2)
  L <- chol(Sigma)
  lambda_oracle <- lambda2_oracle_fn(n)
  rej <- 0
  for (r in seq_len(R)) {
    Z_A <- matrix(rnorm(n * 2), n, 2) %*% L
    Z_B <- matrix(rnorm(n * 2), n, 2) %*% L
    Y_A <- Z_A[, 1] + effect / 2
    Y_B <- Z_B[, 1] - effect / 2
    f_A_L <- Z_A[, 2] + effect / 2
    f_B_L <- Z_B[, 2] - effect / 2
    Z_AU <- matrix(rnorm(N2 * 2), N2, 2) %*% L
    Z_BU <- matrix(rnorm(N2 * 2), N2, 2) %*% L
    f_A_U <- Z_AU[, 2] + effect / 2
    f_B_U <- Z_BU[, 2] - effect / 2

    if (method == "Classical") {
      test <- classical_ttest(Y_A, Y_B, delta0 = 0, alpha = alpha)
    } else if (method == "Vanilla PPI") {
      test <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                            delta0 = 0, alpha = alpha,
                            lambda_A_oracle = 1, lambda_B_oracle = 1)
    } else if (method == "PPI++ (Oracle)") {
      test <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                            delta0 = 0, alpha = alpha,
                            lambda_A_oracle = lambda_oracle,
                            lambda_B_oracle = lambda_oracle)
    } else {
      covA <- cov(Y_A, f_A_L); varA <- var(f_A_L)
      covB <- cov(Y_B, f_B_L); varB <- var(f_B_L)
      lambdaA <- if (is.finite(varA) && varA > 0) covA / ((1 + n / N2) * varA) else 0
      lambdaB <- if (is.finite(varB) && varB > 0) covB / ((1 + n / N2) * varB) else 0
      test <- run_ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                            delta0 = 0, alpha = alpha,
                            lambda_A_oracle = lambdaA,
                            lambda_B_oracle = lambdaB)
    }
    rej <- rej + as.integer(test$reject)
  }
  rej / R
}

# ===================== Paired =====================
N3 <- 5000
delta3 <- 0.3
sigma_D2_3 <- 1
sigma_fD2_3 <- 1
rho3 <- 0.7
cov3 <- rho3 * sqrt(sigma_D2_3 * sigma_fD2_3)

n3_classical <- ceiling(sigma_D2_3 / ((delta3 / (z_alpha + z_beta))^2))
n3_vanilla <- ceiling(power_ppi_mean(delta = delta3, N = N3, n = NULL, power = target_power,
                                     alpha = alpha, sigma_y2 = sigma_D2_3, sigma_f2 = sigma_fD2_3,
                                     cov_y_f = cov3, lambda_mode = "vanilla"))
n3_oracle <- ceiling(power_ppi_paired(delta = delta3, N = N3, n = NULL, power = target_power,
                                      alpha = alpha, sigma_D2 = sigma_D2_3, rho_D = rho3,
                                      sigma_fD2 = sigma_fD2_3))
n3_plugin <- n3_oracle

sim_rate_paired <- function(n, effect, method) {
  Sigma <- matrix(c(sigma_D2_3, cov3, cov3, sigma_fD2_3), 2, 2)
  L <- chol(Sigma)
  lambda_oracle <- cov3 / ((1 + n / N3) * sigma_fD2_3)
  rej <- 0
  for (r in seq_len(R)) {
    Z_L <- matrix(rnorm(n * 2), n, 2) %*% L
    D_L <- Z_L[, 1] + effect
    fD_L <- Z_L[, 2] + effect
    Z_U <- matrix(rnorm(N3 * 2), N3, 2) %*% L
    fD_U <- Z_U[, 2] + effect

    if (method == "Classical") {
      test <- classical_mean_test(D_L, theta0 = 0, alpha = alpha)
    } else if (method == "Vanilla PPI") {
      test <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0, alpha = alpha, lambda_oracle = 1)
    } else if (method == "PPI++ (Oracle)") {
      test <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0, alpha = alpha, lambda_oracle = lambda_oracle)
    } else {
      cov_hat <- cov(D_L, fD_L); varf_hat <- var(fD_L)
      lambda_hat <- if (is.finite(varf_hat) && varf_hat > 0) cov_hat / ((1 + n / N3) * varf_hat) else 0
      test <- run_ppi_mean_test(D_L, fD_L, fD_U, theta0 = 0, alpha = alpha, lambda_oracle = lambda_hat)
    }
    rej <- rej + as.integer(test$reject)
  }
  rej / R
}

build_rows <- function(design, n_vals, sim_fn) {
  methods <- c("Classical", "Vanilla PPI", "PPI++ (Oracle)", "PPI++ (Plugin)")
  out <- data.frame(
    design = design,
    method = methods,
    n_required = as.integer(n_vals),
    achieved_power = NA_real_,
    type1_error = NA_real_,
    stringsAsFactors = FALSE
  )
  n_class <- n_vals[1]
  for (i in seq_len(nrow(out))) {
    out$achieved_power[i] <- sim_fn(out$n_required[i], effect = if (design == "Two-sample") delta2 else if (design == "Paired") delta3 else delta1, method = out$method[i])
    out$type1_error[i] <- sim_fn(out$n_required[i], effect = 0, method = out$method[i])
  }
  out$reduction_vs_classical <- 100 * (1 - out$n_required / n_class)
  out$achieved_power <- round(out$achieved_power, 3)
  out$type1_error <- round(out$type1_error, 3)
  out$reduction_vs_classical <- round(out$reduction_vs_classical, 1)
  out
}

rows_one <- build_rows("One-sample", c(n1_classical, n1_vanilla, n1_oracle, n1_plugin), sim_rate_one)
rows_two <- build_rows("Two-sample", c(n2_classical, n2_vanilla, n2_oracle, n2_plugin), sim_rate_ttest)
rows_pair <- build_rows("Paired", c(n3_classical, n3_vanilla, n3_oracle, n3_plugin), sim_rate_paired)

results <- rbind(rows_one, rows_two, rows_pair)

out_dir <- file.path("simulation_studies", "outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_csv <- file.path(out_dir, "head_to_head_multidesign.csv")
write.csv(results, out_csv, row.names = FALSE)

cat("Saved:\n  ", out_csv, "\n", sep = "")
print(results)
