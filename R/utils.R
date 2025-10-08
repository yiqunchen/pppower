# Compute PP estimate + SE from a single draw (superpopulation)
pp_once <- function(df, N, n, use_sample_var = TRUE, seed = NULL) {
  # df must have columns: y, fhat_cf
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(c("y", "fhat_cf") %in% names(df)))

  idx <- sample.int(nrow(df), size = N + n, replace = FALSE)
  idx_u <- idx[1:N]
  idx_l <- idx[(N+1):(N+n)]

  y_u  <- df$y[idx_u]         # not used, unlabeled in PP estimator
  f_u  <- df$fhat_cf[idx_u]
  y_l  <- df$y[idx_l]
  f_l  <- df$fhat_cf[idx_l]

  # PP estimator
  A_N <- mean(f_u)
  B_n <- mean(y_l - f_l)
  theta_hat <- A_N + B_n

  # SE: superpopulation variance estimate
  if (use_sample_var) {
    # use sample variances from this draw
    var_f   <- stats::var(f_u)         # unbiased (ddof=1)
    var_res <- stats::var(y_l - f_l)
  } else {
    # use "population" variances from the whole df (plug-in)
    var_f   <- stats::var(df$fhat_cf)
    var_res <- stats::var(df$y - df$fhat_cf)
  }
  se <- sqrt(var_f / N + var_res / n)

  list(theta_hat = theta_hat, se = se,
       var_f = var_f, var_res = var_res,
       A_N = A_N, B_n = B_n)
}

# stable least-squares via lm.fit
ols_fit <- function(X, y) {
  # X: n x p (already includes intercept column if desired)
  # y: n vector
  fit <- lm.fit(X, y)
  # residuals and QR object come for free
  list(coef = fit$coefficients,
       resid = fit$residuals,
       qr    = fit$qr)
}

# OLS sandwich estimator
meat_matrix <- function(X, resid) {
  # returns M = E[(eps^2) x x^T] estimated by sample mean
  # M_hat = (1/n) sum_i (resid_i^2) x_i x_i^T
  n <- nrow(X)
  # vectorized: crossprod with residual weights
  WX <- X * as.numeric(resid)
  crossprod(WX) / n
}

bread_inv <- function(X) {
  # (X'X / n)^{-1}
  n <- nrow(X)
  solve(crossprod(X) / n, tol = 1e-12)
}

# k-fold splitting
kfold_split <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  idx <- sample.int(n)
  folds <- split(idx, rep(1:K, length.out = n))
  lapply(folds, sort)
}