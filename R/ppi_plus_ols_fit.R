
#' Internal solver for PPI++ OLS
#'
#' @param X_l Labeled design matrix (n × p).
#' @param Y_l Labeled outcomes (length n).
#' @param f_l Predictions on labeled covariates (length n).
#' @param X_u Unlabeled design matrix (N × p).
#' @param f_u Predictions on unlabeled covariates (length N).
#' @param lambda Optional numeric weight; ignored unless `lambda_type = "user"`.
#' @param lambda_type One of `"plugin"` (default) or `"user"`.
#' @param contrast Optional contrast vector used to derive the plugin lambda.
#' @param clip Logical; if `TRUE`, caps `lambda` inside \eqn{[0,1]} after estimation.
#' @param tol Numerical tolerance passed to linear solves.
#'
#' @return List with `theta_hat`, `V_theta`, `se`, `lambda`, and moment pieces.
#' @keywords internal
ppi_plus_ols <- function(X_l,
                         Y_l,
                         f_l,
                         X_u,
                         f_u,
                         lambda = NULL,
                         lambda_type = c("plugin", "user"),
                         contrast = NULL,
                         clip = TRUE,
                         tol = 1e-12) {
  lambda_type <- match.arg(lambda_type)
  n <- nrow(X_l); N <- nrow(X_u)
  if (ncol(X_l) != ncol(X_u)) stop("X_l and X_u must have matching column counts.", call. = FALSE)

  base_fit <- ppi_ols(X_l, Y_l, f_l, X_u, f_u)
  theta_init <- base_fit$theta_hat
  p <- length(theta_init)

  if (lambda_type == "user") {
    if (is.null(lambda) || length(lambda) != 1L || !is.finite(lambda)) {
      stop("Provide a finite scalar lambda when lambda_type = 'user'.", call. = FALSE)
    }
  } else {
    if (is.null(contrast)) stop("Supply a contrast vector to plug in lambda.", call. = FALSE)
    a <- as.numeric(contrast)
    if (length(a) != p) stop("contrast must have length equal to number of columns in X.", call. = FALSE)

    H_total <- crossprod(rbind(X_l, X_u)) / (n + N)
    H_total_inv <- solve(H_total, tol = tol)

    res_Y  <- Y_l - drop(X_l %*% theta_init)
    res_f_l <- f_l - drop(X_l %*% theta_init)
    res_f_u <- f_u - drop(X_u %*% theta_init)

    Sigma_YY <- meat_matrix(X_l, res_Y)
    Sigma_ff_l <- meat_matrix(X_l, res_f_l)
    Sigma_ff_u <- meat_matrix(X_u, res_f_u)
    Sigma_ff <- (Sigma_ff_l + Sigma_ff_u) / 2

    WX_Y <- X_l * as.numeric(res_Y)
    WX_f <- X_l * as.numeric(res_f_l)
    Sigma_Yf <- crossprod(WX_Y, WX_f) / n

    denom <- (1 + n / N) * as.numeric(t(a) %*% H_total_inv %*% Sigma_ff %*% H_total_inv %*% a)
    if (denom <= 0) {
      lambda <- 0
    } else {
      numer <- as.numeric(t(a) %*% H_total_inv %*% Sigma_Yf %*% H_total_inv %*% a)
      lambda <- numer / denom
    }
  }

  if (isTRUE(clip)) lambda <- max(min(lambda, 1), 0)

  H_L <- crossprod(X_l) / n
  H_U <- crossprod(X_u) / N
  G_L <- crossprod(X_l, Y_l) / n
  G_U <- crossprod(X_u, f_u) / N

  H_mix <- (1 - lambda) * H_L + lambda * H_U
  rhs <- (1 - lambda) * G_L + lambda * G_U
  theta_hat <- drop(solve(H_mix, rhs, tol = tol))

  res_Y  <- Y_l - drop(X_l %*% theta_hat)
  res_f_l <- f_l - drop(X_l %*% theta_hat)
  res_f_u <- f_u - drop(X_u %*% theta_hat)

  Sigma_YY <- meat_matrix(X_l, res_Y)
  Sigma_ff_l <- meat_matrix(X_l, res_f_l)
  Sigma_ff_u <- meat_matrix(X_u, res_f_u)

  WX_Y <- X_l * as.numeric(res_Y)
  WX_f <- X_l * as.numeric(res_f_l)
  Sigma_Yf <- crossprod(WX_Y, WX_f) / n

  middle <- Sigma_YY / n +
    lambda^2 * (Sigma_ff_u / N + Sigma_ff_l / n) -
    2 * lambda * Sigma_Yf / n

  bread_inv <- solve(H_mix, tol = tol)
  V_theta <- bread_inv %*% middle %*% bread_inv
  V_theta <- (V_theta + t(V_theta)) / 2
  se <- sqrt(pmax(diag(V_theta), 0))

  list(
    theta_hat = theta_hat,
    V_theta = V_theta,
    se = se,
    lambda = lambda,
    lambda_type = lambda_type,
    contrast = if (!is.null(contrast)) as.numeric(contrast) else NULL,
    pieces = list(
      H_L = H_L,
      H_U = H_U,
      Sigma_YY = Sigma_YY,
      Sigma_ff_l = Sigma_ff_l,
      Sigma_ff_u = Sigma_ff_u,
      Sigma_Yf = Sigma_Yf,
      n = n,
      N = N
    )
  )
}

#' Fit PPI++ OLS with sandwich variability
#'
#' @param y Outcome column name in `data_l`.
#' @param f_col Column name with predictions in both data sets.
#' @param formula RHS-only formula for covariates.
#' @param data_l Labeled data frame.
#' @param data_u Unlabeled data frame.
#' @param lambda Optional blend weight when `lambda_type = "user"`.
#' @param lambda_type `"plugin"` (default) chooses lambda via the plug-in lambda\eqn{^\star}
#'   for a supplied `contrast`; `"user"` uses the provided `lambda`.
#' @param contrast Numeric vector (length equals number of regressors) specifying
#'   the contrast used for lambda estimation when `lambda_type = "plugin"`.
#' @param clip Logical; if `TRUE`, clips lambda to \eqn{[0,1]}.
#' @param tol Numerical tolerance for linear solves.
#'
#' @return An object of class `c("ppi_plus_ols_fit", "ppi_ols_fit")` containing
#'   coefficients, covariance, standard errors, lambda, and moment pieces.
#' @export
ppi_plus_ols_fit <- function(y,
                             f_col,
                             formula,
                             data_l,
                             data_u,
                             lambda = NULL,
                             lambda_type = c("plugin", "user"),
                             contrast = NULL,
                             clip = TRUE,
                             tol = 1e-12) {
  y <- rlang::as_name(rlang::ensym(y))
  lambda_type <- match.arg(lambda_type)

  X_l <- stats::model.matrix(formula, data_l)
  X_u <- stats::model.matrix(formula, data_u)
  Y_l <- data_l[[y]]
  f_l <- data_l[[f_col]]
  f_u <- data_u[[f_col]]

  fit <- ppi_plus_ols(
    X_l = X_l,
    Y_l = Y_l,
    f_l = f_l,
    X_u = X_u,
    f_u = f_u,
    lambda = lambda,
    lambda_type = lambda_type,
    contrast = contrast,
    clip = clip,
    tol = tol
  )

  structure(
    list(
      coef = fit$theta_hat,
      vcov = fit$V_theta,
      se = fit$se,
      pieces = fit$pieces,
      n = nrow(X_l),
      N = nrow(X_u),
      lambda = fit$lambda,
      lambda_type = fit$lambda_type,
      contrast = fit$contrast,
      call = match.call(),
      terms = stats::terms(formula)
    ),
    class = c("ppi_plus_ols_fit", "ppi_ols_fit") # Existing S3 methods still works
  )
}
