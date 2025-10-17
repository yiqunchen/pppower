#' Power for the PPI++ OLS estimator
#'
#' @param delta Effect size for the contrast \eqn{c^\top(\beta - \beta_0)}.
#' @param contrast Numeric vector (length = number of regressors) identifying the contrast.
#' @param H_L Empirical Gram matrices for labeled designs (per observation).
#' @param H_U Empirical Gram matrices for unlabeled designs (per observation).
#' @param Sigma_YY Variance of the labeled score \eqn{X(Y - X^\top\beta^\star)}.
#' @param Sigma_ff_l Variances of \eqn{X(f(X) - X^\top\beta^\star)} under labeled draws.
#' @param Sigma_ff_u Variances of \eqn{X(f(X) - X^\top\beta^\star)} under unlabeled draws.
#' @param Sigma_Yf Covariance between the labeled score and prediction residual score.
#' @param N Unlabeled sample sizes.
#' @param n Labeled sample sizes.
#' @param lambda Blend weight. Supply the same value used in estimation.
#' @param alpha Two-sided test level.
#' @param tol Numerical tolerance passed to linear solves.
#' 
#' @return Scalar power in \eqn{[0,1]}.
#' @export
power_ppi_pp_ols <- function(delta,
                             contrast,
                             H_L,
                             H_U,
                             Sigma_YY,
                             Sigma_ff_l,
                             Sigma_ff_u,
                             Sigma_Yf,
                             N,
                             n,
                             lambda,
                             alpha = 0.05,
                             tol = 1e-12) {
  stopifnot(length(lambda) == 1L, is.finite(lambda))
  H_mix <- (1 - lambda) * H_L + lambda * H_U
  bread_inv <- solve(H_mix, tol = tol)

  middle <- Sigma_YY / n +
    lambda^2 * (Sigma_ff_u / N + Sigma_ff_l / n) -
    2 * lambda * Sigma_Yf / n

  var_c <- as.numeric(t(contrast) %*% bread_inv %*% middle %*% bread_inv %*% contrast)
  if (var_c < 0 && abs(var_c) < 1e-12) var_c <- 0
  se <- sqrt(var_c)

  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}

#' Monte Carlo power for the PPI++ OLS estimator
#'
#' @param df Data frame with columns `y` and `fhat_cf`.
#' @param formula RHS-only formula describing the regressors.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param c Contrast vector (length equals number of regressors).
#' @param theta0 Null value for the contrast.
#' @param y_col,f_col Column names for the outcome and predictions.
#' @param alpha Two-sided significance level.
#' @param R Number of Monte Carlo replicates.
#' @param seed RNG seed (or `NULL` to leave unchanged).
#' @param lambda Optional user weight (used when `lambda_type = "user"`).
#' @param lambda_type `"plugin"` (default) recomputes \eqn{\lambda^\star} each draw,
#'   `"user"` uses the supplied `lambda`.
#' @param clip If `TRUE`, clips lambda to `[0,1]` inside the solver.
#'
#' @return List with empirical/analytical power, lambda summary, and per-draw details.
#' @export
simulate_power_ppi_pp_ols <- function(df,
                                        formula,
                                        N,
                                        n,
                                        c,
                                        theta0,
                                        y_col = "y",
                                        f_col = "fhat_cf",
                                        alpha = 0.05,
                                        R = 2000,
                                        lambda = NULL,
                                        lambda_type = c("plugin", "user"),
                                        clip = TRUE,
                                        seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  lambda_type <- match.arg(lambda_type)

  df <- as.data.frame(df)
  if (!all(c(y_col, f_col) %in% names(df))) {
    stop("Data frame must contain columns '", y_col, "' and '", f_col, "'.", call. = FALSE)
  }
  if (N <= 0L || n <= 0L) stop("N and n must be positive integers.", call. = FALSE)
  if (N + n > nrow(df)) stop("Requested N + n exceeds number of available observations in df.", call. = FALSE)

  X_full <- stats::model.matrix(formula, df)
  p <- ncol(X_full)
  contrast <- as.numeric(c)
  if (length(contrast) != p) {
    stop("Contrast vector c must have length equal to number of columns in model matrix.", call. = FALSE)
  }

  y_full <- df[[y_col]]
  f_full <- df[[f_col]]

  beta_true <- ols_fit(X_full, y_full)$coef
  contrast_true <- as.numeric(crossprod(contrast, beta_true))
  if (length(theta0) != 1L || !is.finite(theta0)) {
    stop("theta0 must be a finite numeric scalar.", call. = FALSE)
  }
  delta <- contrast_true - theta0

  # Population proxies for variance pieces
  res_y <- y_full - drop(X_full %*% beta_true)
  res_f <- f_full - drop(X_full %*% beta_true)
  H_pop <- crossprod(X_full) / nrow(X_full)
  Sigma_YY_pop <- meat_matrix(X_full, res_y)
  Sigma_ff_pop <- meat_matrix(X_full, res_f)
  Sigma_Yf_pop <- crossprod(X_full * as.numeric(res_y),
                            X_full * as.numeric(res_f)) / nrow(X_full)

  lambda_pop <- switch(
    lambda_type,
    plugin = {
      denom <- (1 + n / N) *
        as.numeric(t(contrast) %*% solve(H_pop) %*% Sigma_ff_pop %*% solve(H_pop) %*% contrast)
      numer <- as.numeric(t(contrast) %*% solve(H_pop) %*% Sigma_Yf_pop %*% solve(H_pop) %*% contrast)
      l_star <- if (denom <= 0) 0 else numer / denom
      if (clip) pmax(pmin(l_star, 1), 0) else l_star
    },
    user = {
      if (is.null(lambda)) stop("Supply lambda when lambda_type = 'user'.", call. = FALSE)
      lambda
    }
  )

  power_exact <- power_ppi_pp_ols(
    delta = delta,
    contrast = contrast,
    H_L = H_pop,
    H_U = H_pop,
    Sigma_YY = Sigma_YY_pop,
    Sigma_ff_l = Sigma_ff_pop,
    Sigma_ff_u = Sigma_ff_pop,
    Sigma_Yf = Sigma_Yf_pop,
    N = N,
    n = n,
    lambda = lambda_pop,
    alpha = alpha
  )

  z_alpha <- stats::qnorm(1 - alpha / 2)
  rej <- logical(R)
  contrasts_hat <- numeric(R)
  ses <- numeric(R)
  lambdas <- numeric(R)

  for (r in seq_len(R)) {
    idx <- sample.int(nrow(df), size = N + n, replace = FALSE)
    idx_u <- idx[seq_len(N)]
    idx_l <- idx[(N + 1):(N + n)]

    X_u <- stats::model.matrix(formula, df[idx_u, , drop = FALSE])
    X_l <- stats::model.matrix(formula, df[idx_l, , drop = FALSE])
    Y_l <- df[[y_col]][idx_l]
    f_l <- df[[f_col]][idx_l]
    f_u <- df[[f_col]][idx_u]

    fit <- ppi_plus_ols(
      X_l = X_l,
      Y_l = Y_l,
      f_l = f_l,
      X_u = X_u,
      f_u = f_u,
      lambda = lambda,
      lambda_type = lambda_type,
      contrast = contrast,
      clip = clip
    )

    contrast_hat <- as.numeric(crossprod(contrast, fit$theta_hat))
    var_c <- as.numeric(t(contrast) %*% fit$V_theta %*% contrast)
    if (var_c < 0 && abs(var_c) < 1e-12) var_c <- 0
    se_c <- sqrt(var_c)

    z <- (contrast_hat - theta0) / se_c
    rej[r] <- abs(z) > z_alpha
    contrasts_hat[r] <- contrast_hat
    ses[r] <- se_c
    lambdas[r] <- fit$lambda
  }

  list(
    beta = beta_true,
    contrast_true = contrast_true,
    theta0 = theta0,
    delta = delta,
    empirical_power = mean(rej),
    analytical_power = power_exact,
    avg_SE = mean(ses),
    lambda_population = lambda_pop,
    avg_lambda_draw = mean(lambdas),
    details = data.frame(
      contrast_hat = contrasts_hat,
      se = ses,
      lambda = lambdas,
      reject = rej
    )
  )
}
