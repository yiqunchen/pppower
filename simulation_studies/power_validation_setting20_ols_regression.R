# =============================================================================
# SETTING 20: OLS REGRESSION CONTRAST
# =============================================================================
# Validates power_ppi_regression() for an OLS beta-coefficient contrast.
#
# DGP:  X ~ N(0, I_2),  Y = X beta + eps,  f = X beta + nu
#        eps ~ N(0,1),  nu = rho*eps + sqrt(1-rho^2)*eta,  eta ~ N(0,1)
# Contrast: a = (1, -1), testing H0: beta_1 - beta_2 = 0
# Population blocks:  H = I_2,  Sigma_YY = I_2,  Sigma_ff = I_2,  Sigma_Yf = rho*I_2
# =============================================================================
cat("Setting 20: OLS Regression Contrast (a = (1,-1))\n")
cat("--------------------------------------------------\n")
setting20_start <- Sys.time()

# Parameters
p <- 2
a <- c(1, -1)
n_values <- c(20, 40, 60, 80, 100)
N_values <- c(200, 500)
rho_values <- c(0.5, 0.7, 0.9)
delta <- 0.3
alpha <- 0.05
R <- 1000

beta_true <- c(delta, 0)   # a'beta = delta

# Population blocks (closed-form for X ~ N(0, I_p), eps ~ N(0,1))
H_pop    <- diag(p)
Sigma_YY <- diag(p)

total_configs_20 <- length(rho_values) * length(n_values) * length(N_values)
current_config_20 <- 0

results_ols_regression <- data.frame()

for (rho in rho_values) {
  Sigma_ff <- diag(p)          # Var(nu) = 1
  Sigma_Yf <- rho * diag(p)   # Cov(eps, nu) = rho

  for (N in N_values) {
    for (n in n_values) {
      current_config_20 <- current_config_20 + 1
      progress_bar(current_config_20, total_configs_20, "Setting 20",
                   sprintf("rho=%.1f, N=%d, n=%d: ", rho, N, n),
                   start_time = setting20_start)

      # Oracle lambda for this contrast
      A_pop <- as.numeric(t(a) %*% Sigma_ff %*% a)   # a'a = 2  (H=I cancels)
      C_pop <- as.numeric(t(a) %*% Sigma_Yf %*% a)   # rho*a'a = 2*rho
      lambda_oracle <- C_pop / ((1 + n / N) * A_pop)
      lambda_oracle <- max(0, min(1, lambda_oracle))

      # Theoretical PPI++ power
      power_theo_ppi <- pppower::power_ppi_regression(
        delta = delta, N = N, n = n, alpha = alpha,
        lambda_mode = "oracle",
        c = a,
        H_L = H_pop, H_U = H_pop,
        Sigma_YY   = Sigma_YY,
        Sigma_ff_l = Sigma_ff, Sigma_ff_u = Sigma_ff,
        Sigma_Yf   = Sigma_Yf,
        warn_smallN = FALSE
      )

      # Theoretical classical power: Var(a'beta_ols) = ||a||^2 / n = 2/n
      se_cl_pop <- sqrt(sum(a^2) / n)
      z_alpha2 <- qnorm(1 - alpha / 2)
      power_theo_classical <- 1 - pnorm(z_alpha2 - abs(delta) / se_cl_pop) +
                                  pnorm(-z_alpha2 - abs(delta) / se_cl_pop)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0

      for (r in 1:R) {
        X_L <- matrix(rnorm(n * p), n, p)
        X_U <- matrix(rnorm(N * p), N, p)

        # Correlated errors
        eps_L <- rnorm(n)
        nu_L  <- rho * eps_L + sqrt(1 - rho^2) * rnorm(n)
        eps_U <- rnorm(N)
        nu_U  <- rho * eps_U + sqrt(1 - rho^2) * rnorm(N)

        Y_L <- X_L %*% beta_true + eps_L
        f_L <- X_L %*% beta_true + nu_L
        f_U <- X_U %*% beta_true + nu_U

        # ---- PPI++ OLS estimator ----
        H_hat_L <- crossprod(X_L) / n
        H_hat_U <- crossprod(X_U) / N
        H_hat_mix <- (1 - lambda_oracle) * H_hat_L + lambda_oracle * H_hat_U
        H_hat_inv <- solve(H_hat_mix)

        beta_ppi <- H_hat_inv %*% (
          crossprod(X_L, Y_L) / n -
          lambda_oracle * (crossprod(X_L, f_L) / n - crossprod(X_U, f_U) / N)
        )

        # Sandwich variance
        res_Y  <- as.vector(Y_L - X_L %*% beta_ppi)
        res_fL <- as.vector(f_L - X_L %*% beta_ppi)
        res_fU <- as.vector(f_U - X_U %*% beta_ppi)

        phi_L <- X_L * (res_Y - lambda_oracle * res_fL)
        psi_U <- lambda_oracle * X_U * res_fU

        V_ppi <- H_hat_inv %*%
          (crossprod(phi_L) / n^2 + crossprod(psi_U) / N^2) %*%
          H_hat_inv

        se_ppi <- sqrt(as.numeric(t(a) %*% V_ppi %*% a))
        z_ppi  <- as.numeric(t(a) %*% beta_ppi) / se_ppi
        rejections_ppi <- rejections_ppi + (abs(z_ppi) > z_alpha2)

        # ---- Classical OLS ----
        beta_ols <- solve(crossprod(X_L)) %*% crossprod(X_L, Y_L)
        res_ols  <- as.vector(Y_L - X_L %*% beta_ols)
        s2_hat   <- sum(res_ols^2) / (n - p)
        V_ols    <- s2_hat * solve(crossprod(X_L))

        se_ols <- sqrt(as.numeric(t(a) %*% V_ols %*% a))
        z_ols  <- as.numeric(t(a) %*% beta_ols) / se_ols
        rejections_classical <- rejections_classical + (abs(z_ols) > z_alpha2)
      }

      results_ols_regression <- rbind(results_ols_regression, data.frame(
        setting = "OLS Regression",
        rho = rho, n = n, N = N, delta = delta,
        lambda_oracle = lambda_oracle,
        power_theo_ppi       = power_theo_ppi,
        power_emp_ppi        = rejections_ppi / R,
        power_theo_classical = power_theo_classical,
        power_emp_classical  = rejections_classical / R
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n",
            difftime(Sys.time(), setting20_start, units = "secs")))
