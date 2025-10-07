
#' Power for prediction-powered linear regression (OLS) estimator
#' @description
#' Computes the analytical (normal-theory) power for a two-sided Wald test on a
#' linear contrast \eqn{c^\top \beta} in the \emph{prediction-powered linear regression (PPI-OLS)} framework.
#' @param delta Numeric scalar. True effect size for the linear contrast
#' \eqn{c^\top (\beta - \beta_0)}.
#' @param V_u Numeric \eqn{p \times p} matrix.
#' Sandwich variance for the unlabeled (imputed) OLS term.
#' @param V_l Numeric \eqn{p \times p} matrix.
#' Sandwich variance for the labeled (rectifier) OLS term.
#' @param N Integer. Unlabeled sample size.
#' @param n Integer. Labeled sample size.
#' @param c Numeric vector of length \eqn{p}.
#' Linear contrast to test; e.g. to test coefficient \eqn{j}, set \code{c = e_j}.
#' @param alpha Numeric scalar in (0, 1). Two-sided significance level.
#'
#' @return
#' Numeric scalar: the analytical (normal-theory) power \eqn{1-\beta}
#' of the two-sided test.
#' 
#' @examples
#' # Example: 2D regression with modest effect
#' set.seed(1)
#' p <- 2
#' V_u <- matrix(c(1.0, 0.2, 0.2, 0.8), ncol = p)
#' V_l <- matrix(c(0.9, 0.1, 0.1, 0.7), ncol = p)
#' N <- 2000
#' n <- 200
#' c <- c(1, 0)     # test first coefficient
#' delta <- 0.15
#'
#' power_ppi_ols(delta, V_u, V_l, N, n, c, alpha = 0.05)
#'
#' # Test second coefficient (different direction)
#' c2 <- c(0, 1)
#' power_ppi_ols(delta = 0.10, V_u, V_l, N, n, c2)
#' 
#' @export

power_ppi_ols <- function(delta, V_u, V_l, N, n, c, alpha = 0.05) {
  # variance of c' theta_hat: v = c'(V_u/N + V_l/n)c
  v <- as.numeric(drop(c %*% (V_u/N + V_l/n) %*% c))
  se <- sqrt(v)
  z_alpha <- qnorm(1 - alpha/2)
  mu <- abs(delta) / se
  # exact two-sided normal power
  1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
}

# Monte Carlo power for PPI-OLS (superpopulation view)
simulate_power_ppi_ols <- function(df,
                                   formula,
                                   N,
                                   n,
                                   c,
                                   theta0,
                                   y_col = "y",
                                   f_col = "fhat_cf",
                                   alpha = 0.05,
                                   R = 2000,
                                   seed = 1) {
  if (!is.null(seed)) set.seed(seed)

  df <- as.data.frame(df)
  if (!all(c(y_col, f_col) %in% names(df))) {
    stop("Data frame must contain columns '", y_col, "' and '", f_col, "'.", call. = FALSE)
  }
  if (N <= 0L || n <= 0L) stop("N and n must be positive integers.", call. = FALSE)
  if (N + n > nrow(df)) {
    stop("Requested N + n exceeds number of available observations in df.", call. = FALSE)
  }
  if (length(theta0) != 1L || !is.finite(theta0)) {
    stop("theta0 must be a finite numeric scalar.", call. = FALSE)
  }

  X_full <- stats::model.matrix(formula, df)
  p <- ncol(X_full)
  c_vec <- as.numeric(c)
  if (length(c_vec) != p) {
    stop("Contrast vector c must have length equal to number of columns in model matrix.", call. = FALSE)
  }

  y_full <- df[[y_col]]
  f_full <- df[[f_col]]

  beta_true <- ols_fit(X_full, y_full)$coef
  contrast_true <- as.numeric(crossprod(c_vec, beta_true))
  delta <- contrast_true - theta0

  # Plug-in V_u and V_l using full data
  iSigma_full <- bread_inv(X_full)
  fit_u_full <- ols_fit(X_full, f_full)
  fit_l_full <- ols_fit(X_full, f_full - y_full)
  V_u_pop <- iSigma_full %*% meat_matrix(X_full, fit_u_full$resid) %*% iSigma_full
  V_l_pop <- iSigma_full %*% meat_matrix(X_full, fit_l_full$resid) %*% iSigma_full

  power_exact <- power_ppi_ols(delta, V_u_pop, V_l_pop, N, n, c_vec, alpha)

  rej <- logical(R)
  contrasts_hat <- numeric(R)
  ses <- numeric(R)
  theta_hats <- matrix(NA_real_, nrow = R, ncol = p)
  colnames(theta_hats) <- colnames(X_full)

  z_alpha <- stats::qnorm(1 - alpha/2)

  for (r in seq_len(R)) {
    idx <- sample.int(nrow(df), size = N + n, replace = FALSE)
    idx_u <- idx[seq_len(N)]
    idx_l <- idx[(N + 1):(N + n)]

    X_u <- stats::model.matrix(formula, df[idx_u, , drop = FALSE])
    X_l <- stats::model.matrix(formula, df[idx_l, , drop = FALSE])
    Y_l <- df[[y_col]][idx_l]
    f_l <- df[[f_col]][idx_l]
    f_u <- df[[f_col]][idx_u]

    fit <- ppi_ols(X_l, Y_l, f_l, X_u, f_u)
    theta_hat <- fit$theta_hat
    theta_hats[r, ] <- theta_hat

    contrast_hat <- as.numeric(crossprod(c_vec, theta_hat))
    var_c <- as.numeric(t(c_vec) %*% fit$V_theta %*% c_vec)
    if (var_c < 0 && abs(var_c) < 1e-12) var_c <- 0
    se_c <- sqrt(var_c)

    z <- (contrast_hat - theta0) / se_c
    rej[r] <- abs(z) > z_alpha
    contrasts_hat[r] <- contrast_hat
    ses[r] <- se_c
  }

  list(
    beta = beta_true,
    contrast_true = contrast_true,
    theta0 = theta0,
    delta = delta,
    empirical_power = mean(rej),
    analytical_power = power_exact,
    avg_SE = mean(ses),
    details = data.frame(contrast_hat = contrasts_hat, se = ses, reject = rej),
    theta_hat = theta_hats
  )
}
