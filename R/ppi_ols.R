ppi_ols <- function(X_l, Y_l, f_l, X_u, f_u) {
  # labeled (n): regress f_l - Y_l on X_l  -> delta_hat
  fit_l <- ols_fit(X_l, f_l - Y_l)
  delta_hat <- fit_l$coef

  # unlabeled (N): regress f_u on X_u      -> theta_tilde
  fit_u <- ols_fit(X_u, f_u)
  theta_tilde <- fit_u$coef

  # PPI–OLS estimator
  theta_hat <- theta_tilde - delta_hat

  # sandwich variance pieces
  n <- nrow(X_l); N <- nrow(X_u)

  # labeled rectifier: V = iSigma * M * iSigma
  iSigma_l <- bread_inv(X_l)
  M_l      <- meat_matrix(X_l, fit_l$resid)   # uses resid_l of (f_l - Y_l) ~ X_l
  V_l      <- iSigma_l %*% M_l %*% iSigma_l   # per-observation covariance

  # unlabeled imputation: V_tilde = iSigma * M * iSigma
  iSigma_u <- bread_inv(X_u)
  M_u      <- meat_matrix(X_u, fit_u$resid)   # uses resid_u of f_u ~ X_u
  V_u      <- iSigma_u %*% M_u %*% iSigma_u   # per-observation covariance

  # full cov of theta_hat
  V_theta  <- V_u / N + V_l / n
  se       <- sqrt(diag(V_theta))

  list(theta_hat = as.numeric(theta_hat),
       V_theta   = V_theta,
       se        = se,
       pieces    = list(V_u = V_u, V_l = V_l, n = n, N = N))
}


# Analytical power for OLS
# Under H1: c' beta = c' beta0 + delta; test H0: c' beta = c' beta0 (two-sided)
power_ppi_ols <- function(delta, V_u, V_l, N, n, c, alpha = 0.05) {
  # variance of c' theta_hat: v = c'(V_u/N + V_l/n)c
  v <- as.numeric(drop(c %*% (V_u/N + V_l/n) %*% c))
  se <- sqrt(v)
  z_alpha <- qnorm(1 - alpha/2)
  mu <- abs(delta) / se
  # exact two-sided normal power
  1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
}