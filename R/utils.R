# Compute PP estimate + SE from a single draw (superpopulation view)
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