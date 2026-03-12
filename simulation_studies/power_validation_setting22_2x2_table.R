# =============================================================================
# SETTING 22: 2x2 CONTINGENCY TABLE (ODDS RATIO AND RISK RATIO)
# =============================================================================
# Validates power_ppi_2x2() for both OR (logistic regression) and RR
# (delta-method on group PPI++ means) against Monte Carlo simulations.
#
# DGP:  X ~ Bernoulli(prev_exp),  Y|X ~ Bernoulli(p_exp or p_ctrl)
#        f|Y ~ Bernoulli(sens if Y=1, 1-spec if Y=0)
# =============================================================================
cat("Setting 22: 2x2 Contingency Table (OR and RR)\n")
cat("------------------------------------------------\n")
setting22_start <- Sys.time()

library(pppower)

# Parameters
p_ctrl_22    <- 0.20
p_exp_values <- c(0.30, 0.35, 0.40)   # OR ~ 1.71, 2.15, 2.67
accuracy_values_22 <- c(0.7, 0.8, 0.9)  # sens = spec = accuracy
n_values_22  <- c(20, 40, 60, 80, 100)
N_values_22  <- c(200, 500)
prev_exp_22  <- 0.5
alpha_22     <- 0.05
R_22         <- 3000

# ---- OR validation ----
grid_22_or <- expand.grid(
  p_exp = p_exp_values,
  accuracy = accuracy_values_22,
  N = N_values_22,
  n = n_values_22,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_22_or <- nrow(grid_22_or)

results_2x2_or <- vector("list", total_22_or)

for (i in seq_len(total_22_or)) {
  row <- grid_22_or[i, ]
  sens_i <- row$accuracy
  spec_i <- row$accuracy

  progress_bar(i, total_22_or, "Setting 22 (OR)",
               sprintf("p_exp=%.2f, acc=%.0f%%, N=%d, n=%d: ",
                       row$p_exp, row$accuracy * 100, row$N, row$n),
               start_time = setting22_start)

  # Theoretical power (OR)
  power_theo_or <- tryCatch(
    power_ppi_2x2(
      p_exp = row$p_exp, p_ctrl = p_ctrl_22,
      N = row$N, n = row$n, alpha = alpha_22,
      prev_exp = prev_exp_22,
      sens = sens_i, spec = spec_i,
      effect_measure = "OR",
      lambda_mode = "oracle",
      warn_smallN = FALSE
    ),
    error = function(e) NA_real_
  )

  # Classical logistic power (no surrogate): Fisher-based
  beta0 <- log(p_ctrl_22 / (1 - p_ctrl_22))
  beta1 <- log(row$p_exp / (1 - row$p_exp)) - beta0
  OR_true <- exp(beta1)

  pi <- prev_exp_22
  w0 <- p_ctrl_22 * (1 - p_ctrl_22)
  w1 <- row$p_exp * (1 - row$p_exp)
  H_pop <- matrix(c((1 - pi) * w0 + pi * w1, pi * w1,
                     pi * w1, pi * w1), 2, 2)
  H_inv <- solve(H_pop)
  se_cl <- sqrt(H_inv[2, 2] / row$n)
  z_alpha2 <- qnorm(1 - alpha_22 / 2)
  power_classical <- 1 - pnorm(z_alpha2 - abs(beta1) / se_cl) +
    pnorm(-z_alpha2 - abs(beta1) / se_cl)

  # Surrogate prevalence per group
  p_f_ctrl <- sens_i * p_ctrl_22 + (1 - spec_i) * (1 - p_ctrl_22)
  p_f_exp  <- sens_i * row$p_exp  + (1 - spec_i) * (1 - row$p_exp)

  # Oracle lambda for contrast c = (0,1) from logistic blocks
  blk <- pppower:::.ppi_2x2_logistic_blocks(row$p_exp, p_ctrl_22, prev_exp_22, sens_i, spec_i)
  Hinv_pop <- solve(blk$H)
  cvec <- c(0, 1)
  A_pop <- as.numeric(t(cvec) %*% Hinv_pop %*% blk$Sigma_ff %*% Hinv_pop %*% cvec)
  C_pop <- as.numeric(t(cvec) %*% Hinv_pop %*% blk$Sigma_Yf %*% Hinv_pop %*% cvec)
  lambda_oracle <- max(0, min(1, C_pop / ((1 + row$n / row$N) * A_pop)))

  # Monte Carlo (OR via PPI++ logistic regression)
  rejections_or <- 0
  valid_or <- 0

  for (r in seq_len(R_22)) {
    # Generate data
    n_exp <- round(row$n * prev_exp_22)
    n_ctrl <- row$n - n_exp
    N_exp <- round(row$N * prev_exp_22)
    N_ctrl <- row$N - N_exp

    X_L <- c(rep(0, n_ctrl), rep(1, n_exp))
    Y_L <- c(rbinom(n_ctrl, 1, p_ctrl_22), rbinom(n_exp, 1, row$p_exp))
    f_L <- ifelse(Y_L == 1, rbinom(row$n, 1, sens_i), rbinom(row$n, 1, 1 - spec_i))

    X_U <- c(rep(0, N_ctrl), rep(1, N_exp))
    Y_U_lat <- c(rbinom(N_ctrl, 1, p_ctrl_22), rbinom(N_exp, 1, row$p_exp))
    f_U <- ifelse(Y_U_lat == 1, rbinom(row$N, 1, sens_i), rbinom(row$N, 1, 1 - spec_i))

    # Design matrices: [1, X]
    DM_L <- cbind(1, X_L)
    DM_U <- cbind(1, X_U)

    # PPI++ logistic regression via Newton-Raphson
    ppi_ok <- tryCatch({
      theta <- suppressWarnings(glm.fit(DM_L, Y_L, family = binomial())$coefficients)
      if (any(is.na(theta))) stop("init failed")

      for (iter in 1:30) {
        mu_L <- plogis(DM_L %*% theta)
        mu_U <- plogis(DM_U %*% theta)
        w_L <- as.vector(mu_L * (1 - mu_L))
        w_U <- as.vector(mu_U * (1 - mu_U))

        g <- crossprod(DM_L, Y_L - mu_L) / row$n -
          lambda_oracle * (crossprod(DM_L, f_L - mu_L) / row$n -
                             crossprod(DM_U, f_U - mu_U) / row$N)

        H_mix <- (1 - lambda_oracle) * crossprod(DM_L * sqrt(w_L)) / row$n +
          lambda_oracle * crossprod(DM_U * sqrt(w_U)) / row$N

        d_theta <- solve(H_mix, g)
        theta <- theta + d_theta
        if (max(abs(d_theta)) < 1e-8) break
      }

      mu_L <- plogis(DM_L %*% theta)
      mu_U <- plogis(DM_U %*% theta)
      w_L <- as.vector(mu_L * (1 - mu_L))
      w_U <- as.vector(mu_U * (1 - mu_U))

      H_mix <- (1 - lambda_oracle) * crossprod(DM_L * sqrt(w_L)) / row$n +
        lambda_oracle * crossprod(DM_U * sqrt(w_U)) / row$N
      H_mix_inv <- solve(H_mix)

      phi_L <- DM_L * (as.vector(Y_L - mu_L) -
                          lambda_oracle * as.vector(f_L - mu_L))
      psi_U <- lambda_oracle * DM_U * as.vector(f_U - mu_U)

      V_ppi <- H_mix_inv %*%
        (crossprod(phi_L) / row$n^2 + crossprod(psi_U) / row$N^2) %*%
        H_mix_inv

      se_ppi <- sqrt(as.numeric(t(cvec) %*% V_ppi %*% cvec))
      z_ppi <- as.numeric(t(cvec) %*% theta) / se_ppi
      list(reject = abs(z_ppi) > z_alpha2)
    }, error = function(e) NULL)

    if (!is.null(ppi_ok)) {
      valid_or <- valid_or + 1
      rejections_or <- rejections_or + ppi_ok$reject
    }
  }

  results_2x2_or[[i]] <- data.frame(
    setting = "2x2 Table (OR)",
    p_exp = row$p_exp,
    p_ctrl = p_ctrl_22,
    OR = OR_true,
    accuracy = row$accuracy,
    n = row$n,
    N = row$N,
    delta = beta1,
    lambda_oracle = lambda_oracle,
    power_theo_ppi = as.numeric(power_theo_or),
    power_emp_ppi = rejections_or / max(valid_or, 1),
    power_theo_classical = power_classical,
    stringsAsFactors = FALSE
  )
}

results_2x2_or <- do.call(rbind, results_2x2_or)

cat(sprintf("\n  OR validation done. Max |theo - emp|: %.4f\n",
            max(abs(results_2x2_or$power_theo_ppi - results_2x2_or$power_emp_ppi), na.rm = TRUE)))

# ---- RR validation ----
cat("  Starting RR validation...\n")
setting22_rr_start <- Sys.time()

grid_22_rr <- expand.grid(
  p_exp = p_exp_values,
  accuracy = accuracy_values_22,
  N = N_values_22,
  n = n_values_22,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
total_22_rr <- nrow(grid_22_rr)

results_2x2_rr <- vector("list", total_22_rr)

for (i in seq_len(total_22_rr)) {
  row <- grid_22_rr[i, ]
  sens_i <- row$accuracy
  spec_i <- row$accuracy

  progress_bar(i, total_22_rr, "Setting 22 (RR)",
               sprintf("p_exp=%.2f, acc=%.0f%%, N=%d, n=%d: ",
                       row$p_exp, row$accuracy * 100, row$N, row$n),
               start_time = setting22_rr_start)

  RR_true <- row$p_exp / p_ctrl_22
  delta_rr <- log(RR_true)

  # Theoretical power (RR)
  power_theo_rr <- tryCatch(
    power_ppi_2x2(
      p_exp = row$p_exp, p_ctrl = p_ctrl_22,
      N = row$N, n = row$n, alpha = alpha_22,
      prev_exp = prev_exp_22,
      sens = sens_i, spec = spec_i,
      effect_measure = "RR",
      lambda_mode = "oracle",
      warn_smallN = FALSE
    ),
    error = function(e) NA_real_
  )

  # Classical RR power (delta method, no surrogate)
  se_cl_rr <- sqrt((1 - row$p_exp) / (row$p_exp * row$n * prev_exp_22) +
                      (1 - p_ctrl_22) / (p_ctrl_22 * row$n * (1 - prev_exp_22)))
  z_alpha2 <- qnorm(1 - alpha_22 / 2)
  power_classical_rr <- 1 - pnorm(z_alpha2 - abs(delta_rr) / se_cl_rr) +
    pnorm(-z_alpha2 - abs(delta_rr) / se_cl_rr)

  # Monte Carlo (RR via per-group PPI++ means + delta method)
  mom_exp  <- binary_moments_from_sens_spec(p = row$p_exp, sens = sens_i, spec = spec_i)
  mom_ctrl <- binary_moments_from_sens_spec(p = p_ctrl_22, sens = sens_i, spec = spec_i)

  rejections_rr <- 0

  for (r in seq_len(R_22)) {
    n_exp <- round(row$n * prev_exp_22)
    n_ctrl <- row$n - n_exp
    N_exp <- round(row$N * prev_exp_22)
    N_ctrl <- row$N - N_exp

    # Exposed group
    Y_exp_L <- rbinom(n_exp, 1, row$p_exp)
    f_exp_L <- ifelse(Y_exp_L == 1, rbinom(n_exp, 1, sens_i), rbinom(n_exp, 1, 1 - spec_i))
    Y_exp_U <- rbinom(N_exp, 1, row$p_exp)
    f_exp_U <- ifelse(Y_exp_U == 1, rbinom(N_exp, 1, sens_i), rbinom(N_exp, 1, 1 - spec_i))

    # Control group
    Y_ctrl_L <- rbinom(n_ctrl, 1, p_ctrl_22)
    f_ctrl_L <- ifelse(Y_ctrl_L == 1, rbinom(n_ctrl, 1, sens_i), rbinom(n_ctrl, 1, 1 - spec_i))
    Y_ctrl_U <- rbinom(N_ctrl, 1, p_ctrl_22)
    f_ctrl_U <- ifelse(Y_ctrl_U == 1, rbinom(N_ctrl, 1, sens_i), rbinom(N_ctrl, 1, 1 - spec_i))

    # PPI++ mean estimates per group with oracle lambda
    lam_exp <- mom_exp$cov_y_f / ((1 + n_exp / N_exp) * mom_exp$sigma_f2)
    lam_ctrl <- mom_ctrl$cov_y_f / ((1 + n_ctrl / N_ctrl) * mom_ctrl$sigma_f2)

    p_hat_exp <- mean(Y_exp_L) - lam_exp * (mean(f_exp_L) - mean(f_exp_U))
    p_hat_ctrl <- mean(Y_ctrl_L) - lam_ctrl * (mean(f_ctrl_L) - mean(f_ctrl_U))

    # Clamp to avoid log of non-positive
    p_hat_exp  <- max(p_hat_exp, 1e-6)
    p_hat_ctrl <- max(p_hat_ctrl, 1e-6)

    log_rr_hat <- log(p_hat_exp) - log(p_hat_ctrl)

    # Variance estimate via delta method
    # Var(p_hat_x) from PPI++ sandwich
    res_Y_exp <- Y_exp_L - p_hat_exp
    res_f_exp_L <- f_exp_L - p_hat_exp
    res_f_exp_U <- f_exp_U - p_hat_exp

    phi_exp <- res_Y_exp - lam_exp * res_f_exp_L
    psi_exp <- lam_exp * res_f_exp_U

    var_p_exp <- var(phi_exp) / n_exp + var(psi_exp) / N_exp

    res_Y_ctrl <- Y_ctrl_L - p_hat_ctrl
    res_f_ctrl_L <- f_ctrl_L - p_hat_ctrl
    res_f_ctrl_U <- f_ctrl_U - p_hat_ctrl

    phi_ctrl <- res_Y_ctrl - lam_ctrl * res_f_ctrl_L
    psi_ctrl <- lam_ctrl * res_f_ctrl_U

    var_p_ctrl <- var(phi_ctrl) / n_ctrl + var(psi_ctrl) / N_ctrl

    var_log_rr <- var_p_exp / p_hat_exp^2 + var_p_ctrl / p_hat_ctrl^2
    se_rr <- sqrt(max(var_log_rr, 0))

    if (se_rr > 0) {
      z_rr <- log_rr_hat / se_rr
      rejections_rr <- rejections_rr + (abs(z_rr) > z_alpha2)
    }
  }

  results_2x2_rr[[i]] <- data.frame(
    setting = "2x2 Table (RR)",
    p_exp = row$p_exp,
    p_ctrl = p_ctrl_22,
    RR = RR_true,
    accuracy = row$accuracy,
    n = row$n,
    N = row$N,
    delta = delta_rr,
    power_theo_ppi = as.numeric(power_theo_rr),
    power_emp_ppi = rejections_rr / R_22,
    power_theo_classical = power_classical_rr,
    stringsAsFactors = FALSE
  )
}

results_2x2_rr <- do.call(rbind, results_2x2_rr)

cat(sprintf("\n  RR validation done. Max |theo - emp|: %.4f\n",
            max(abs(results_2x2_rr$power_theo_ppi - results_2x2_rr$power_emp_ppi), na.rm = TRUE)))
cat(sprintf("Setting 22 completed. (%.1f seconds)\n\n",
            difftime(Sys.time(), setting22_start, units = "secs")))
