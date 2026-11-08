# =============================================================================
# SETTING 21: GLM LOGISTIC REGRESSION CONTRAST
# =============================================================================
# Validates power_ppi_regression() for a logistic-regression beta contrast.
#
# DGP:  X ~ N(0, I_2),  P(Y=1|X) = expit(X'beta)
#        beta_true = (delta, 0),  contrast a = (1, -1)
#        Surrogate: noisy binary classifier with sens = spec = accuracy
# Population blocks estimated from a large pilot via compute_ppi_blocks()
# for the analytical curve; empirical power uses a 2-fold cross-fitted
# plug-in lambda estimated within each replicate.
# =============================================================================
cat("Setting 21: GLM Logistic Regression Contrast (a = (1,-1))\n")
cat("-----------------------------------------------------------\n")
setting21_start <- Sys.time()

# Parameters
p <- 2
a <- c(1, -1)
n_values <- c(20, 40, 60, 80, 100)
N_values <- c(200, 500)
accuracy_values <- c(0.7, 0.8, 0.9)   # sens = spec = accuracy
delta_log <- 0.5
alpha <- 0.05
R <- 2000  # R1 revision: was 1000
M_pilot <- 100000   # pilot size for population blocks

beta_true <- c(delta_log, 0)

clip_prob <- function(x, eps = 1e-6) {
  pmin(pmax(as.numeric(x), eps), 1 - eps)
}

safe_solve <- function(M, ridge = 1e-8) {
  solve(M + diag(ridge, nrow(M)))
}

estimate_plugin_lambda_glm_cf <- function(X_L, Y_L, f_L, X_U, f_U, a, K = 2) {
  n <- nrow(X_L)
  N <- nrow(X_U)
  fold_id <- sample(rep(seq_len(K), length.out = n))

  mu_L_cf <- numeric(n)
  mu_U_cf_mat <- matrix(NA_real_, N, K)

  for (k in seq_len(K)) {
    te <- which(fold_id == k)
    tr <- setdiff(seq_len(n), te)

    fit_k <- suppressWarnings(glm.fit(
      x = X_L[tr, , drop = FALSE],
      y = Y_L[tr],
      family = binomial()
    ))
    beta_k <- fit_k$coefficients
    if (any(!is.finite(beta_k))) {
      return(0)
    }

    mu_L_cf[te] <- clip_prob(plogis(X_L[te, , drop = FALSE] %*% beta_k))
    mu_U_cf_mat[, k] <- clip_prob(plogis(X_U %*% beta_k))
  }

  mu_U_cf <- rowMeans(mu_U_cf_mat)
  w_U <- pmax(mu_U_cf * (1 - mu_U_cf), 1e-6)
  H_U_hat <- crossprod(X_U * sqrt(w_U)) / N

  res_Y_l <- Y_L - mu_L_cf
  res_f_l <- f_L - mu_L_cf
  res_f_u <- f_U - mu_U_cf

  Sigma_ff_u_hat <- crossprod(X_U * res_f_u) / N
  Sigma_Yf_hat <- crossprod(X_L * res_Y_l, X_L * res_f_l) / n

  Hinv_hat <- tryCatch(safe_solve(H_U_hat), error = function(e) NULL)
  if (is.null(Hinv_hat)) {
    return(0)
  }

  A_hat <- as.numeric(t(a) %*% Hinv_hat %*% Sigma_ff_u_hat %*% Hinv_hat %*% a)
  C_hat <- as.numeric(t(a) %*% Hinv_hat %*% Sigma_Yf_hat %*% Hinv_hat %*% a)

  if (!is.finite(A_hat) || A_hat <= 1e-10 || !is.finite(C_hat)) {
    return(0)
  }

  lambda_hat <- C_hat / ((1 + n / N) * A_hat)
  max(0, min(1, lambda_hat))
}

# ---- Pre-compute population blocks from large pilot per accuracy level ----
cat("  Pre-computing population blocks from pilot data...\n")
rng_state <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) .Random.seed else NULL
blocks_list <- list()

for (acc in accuracy_values) {
  set.seed(42)   # deterministic pilot

  X_pilot  <- matrix(rnorm(M_pilot * p), M_pilot, p)
  mu_pilot <- plogis(X_pilot %*% beta_true)
  Y_pilot  <- rbinom(M_pilot, 1, mu_pilot)
  f_pilot  <- ifelse(Y_pilot == 1,
                      rbinom(M_pilot, 1, acc),
                      rbinom(M_pilot, 1, 1 - acc))

  half <- M_pilot %/% 2
  idx_L <- 1:half
  idx_U <- (half + 1):M_pilot

  blk <- pppower::compute_ppi_blocks(
    model_type = "glm",
    X_l = X_pilot[idx_L, ], Y_l = Y_pilot[idx_L], f_l = f_pilot[idx_L],
    X_u = X_pilot[idx_U, ],                        f_u = f_pilot[idx_U],
    beta   = beta_true,
    family = "binomial"
  )
  blocks_list[[as.character(acc)]] <- blk
  cat(sprintf("    accuracy = %.0f%%: blocks computed\n", acc * 100))
}

if (!is.null(rng_state)) {
  .Random.seed <- rng_state
}

# ---- Main simulation loop ----
total_configs_21 <- length(accuracy_values) * length(n_values) * length(N_values)
current_config_21 <- 0

results_glm_logistic <- data.frame()

for (acc in accuracy_values) {
  blk <- blocks_list[[as.character(acc)]]
  H_inv_pop <- solve(blk$H_U)

  for (N in N_values) {
    for (n in n_values) {
      current_config_21 <- current_config_21 + 1
      progress_bar(current_config_21, total_configs_21, "Setting 21",
                   sprintf("acc=%.0f%%, N=%d, n=%d: ", acc * 100, N, n),
                   start_time = setting21_start)

      # Oracle lambda from pilot blocks
      A_pop <- as.numeric(t(a) %*% H_inv_pop %*% blk$Sigma_ff_u %*% H_inv_pop %*% a)
      C_pop <- as.numeric(t(a) %*% H_inv_pop %*% blk$Sigma_Yf   %*% H_inv_pop %*% a)
      lambda_oracle <- max(0, min(1, C_pop / ((1 + n / N) * A_pop)))

      # Theoretical PPI++ power
      power_theo_ppi <- pppower::power_ppi_regression(
        delta = delta_log, N = N, n = n, alpha = alpha,
        lambda_mode = "oracle",
        c = a,
        H_L = blk$H_L, H_U = blk$H_U,
        Sigma_YY   = blk$Sigma_YY,
        Sigma_ff_l = blk$Sigma_ff_l, Sigma_ff_u = blk$Sigma_ff_u,
        Sigma_Yf   = blk$Sigma_Yf,
        warn_smallN = FALSE
      )

      # Theoretical classical power (Fisher-based)
      V_cl_pop <- H_inv_pop / n
      se_cl_pop <- sqrt(as.numeric(t(a) %*% V_cl_pop %*% a))
      z_alpha2 <- qnorm(1 - alpha / 2)
      power_theo_classical <- 1 - pnorm(z_alpha2 - abs(delta_log) / se_cl_pop) +
                                  pnorm(-z_alpha2 - abs(delta_log) / se_cl_pop)

      # Monte Carlo simulation
      rejections_ppi <- 0
      rejections_classical <- 0
      valid_ppi <- 0
      valid_cl  <- 0
      lambda_plugin_sum <- 0

      for (r in 1:R) {
        X_L <- matrix(rnorm(n * p), n, p)
        X_U <- matrix(rnorm(N * p), N, p)

        mu_L <- plogis(X_L %*% beta_true)
        mu_U <- plogis(X_U %*% beta_true)
        Y_L  <- rbinom(n, 1, mu_L)
        Y_U_latent <- rbinom(N, 1, mu_U)

        f_L <- ifelse(Y_L == 1,        rbinom(n, 1, acc), rbinom(n, 1, 1 - acc))
        f_U <- ifelse(Y_U_latent == 1, rbinom(N, 1, acc), rbinom(N, 1, 1 - acc))

        lambda_plugin <- estimate_plugin_lambda_glm_cf(
          X_L = X_L, Y_L = Y_L, f_L = f_L,
          X_U = X_U, f_U = f_U, a = a, K = 2
        )

        # ---- PPI++ GLM logistic (Newton-Raphson) ----
        ppi_ok <- tryCatch({
          theta <- suppressWarnings(glm.fit(X_L, Y_L, family = binomial())$coefficients)
          if (any(is.na(theta))) stop("init failed")

          for (iter in 1:30) {
            mu_hat_L <- plogis(X_L %*% theta)
            mu_hat_U <- plogis(X_U %*% theta)
            w_L <- as.vector(mu_hat_L * (1 - mu_hat_L))
            w_U <- as.vector(mu_hat_U * (1 - mu_hat_U))

            g <- crossprod(X_L, Y_L - mu_hat_L) / n -
                 lambda_plugin * (crossprod(X_L, f_L - mu_hat_L) / n -
                                  crossprod(X_U, f_U - mu_hat_U) / N)

            H_mix <- (1 - lambda_plugin) * crossprod(X_L * sqrt(w_L)) / n +
                      lambda_plugin      * crossprod(X_U * sqrt(w_U)) / N

            d_theta <- solve(H_mix, g)
            theta   <- theta + d_theta
            if (max(abs(d_theta)) < 1e-8) break
          }

          mu_hat_L <- plogis(X_L %*% theta)
          mu_hat_U <- plogis(X_U %*% theta)
          w_L <- as.vector(mu_hat_L * (1 - mu_hat_L))
          w_U <- as.vector(mu_hat_U * (1 - mu_hat_U))

          H_mix <- (1 - lambda_plugin) * crossprod(X_L * sqrt(w_L)) / n +
                    lambda_plugin      * crossprod(X_U * sqrt(w_U)) / N
          H_mix_inv <- solve(H_mix)

          phi_L <- X_L * (as.vector(Y_L - mu_hat_L) -
                           lambda_plugin * as.vector(f_L - mu_hat_L))
          psi_U <- lambda_plugin * X_U * as.vector(f_U - mu_hat_U)

          V_ppi <- H_mix_inv %*%
            (crossprod(phi_L) / n^2 + crossprod(psi_U) / N^2) %*%
            H_mix_inv

          se_ppi <- sqrt(as.numeric(t(a) %*% V_ppi %*% a))
          z_ppi  <- as.numeric(t(a) %*% theta) / se_ppi
          list(reject = abs(z_ppi) > z_alpha2)
        }, error = function(e) NULL)

        if (!is.null(ppi_ok)) {
          valid_ppi <- valid_ppi + 1
          rejections_ppi <- rejections_ppi + ppi_ok$reject
          lambda_plugin_sum <- lambda_plugin_sum + lambda_plugin
        }

        # ---- Classical GLM (Fisher-based SE) ----
        cl_ok <- tryCatch({
          fit_cl  <- suppressWarnings(glm.fit(X_L, Y_L, family = binomial()))
          beta_cl <- fit_cl$coefficients
          if (any(is.na(beta_cl))) stop("glm failed")
          mu_cl   <- plogis(X_L %*% beta_cl)
          w_cl    <- as.vector(mu_cl * (1 - mu_cl))
          V_cl    <- solve(crossprod(X_L * sqrt(w_cl)))
          se_cl <- sqrt(as.numeric(t(a) %*% V_cl %*% a))
          z_cl  <- as.numeric(t(a) %*% beta_cl) / se_cl
          list(reject = abs(z_cl) > z_alpha2)
        }, error = function(e) NULL)

        if (!is.null(cl_ok)) {
          valid_cl <- valid_cl + 1
          rejections_classical <- rejections_classical + cl_ok$reject
        }
      }

      results_glm_logistic <- rbind(results_glm_logistic, data.frame(
        setting  = "GLM Logistic",
        accuracy = acc, n = n, N = N, delta = delta_log,
        lambda_oracle        = lambda_oracle,
        lambda_plugin_mean   = lambda_plugin_sum / max(valid_ppi, 1),
        power_theo_ppi       = power_theo_ppi,
        power_emp_ppi        = rejections_ppi / max(valid_ppi, 1),
        power_theo_classical = power_theo_classical,
        power_emp_classical  = rejections_classical / max(valid_cl, 1),
        mc_reps = max(valid_ppi, 1)
      ))
    }
  }
}

cat(sprintf("Done. (%.1f seconds)\n\n",
            difftime(Sys.time(), setting21_start, units = "secs")))
