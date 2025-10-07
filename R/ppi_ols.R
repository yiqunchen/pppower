#' PPI–OLS: Internal estimator and sandwich variance
#'
#' Fits the Prediction-Powered Inference (PPI) OLS estimator
#' \eqn{\hat\theta = \tilde\theta - \hat\delta}, where
#' \eqn{\hat\delta} comes from regressing \eqn{f_l - Y_l} on \eqn{X_l}
#' (labeled data) and \eqn{\tilde\theta} comes from regressing \eqn{f_u} on
#' \eqn{X_u} (unlabeled data). Returns the estimate and a sandwich covariance
#' \eqn{\mathrm{Var}(\hat\theta) = V_u/N + V_l/n}.
#'
#' @param X_l Matrix of labeled covariates (n × p).
#' @param Y_l Vector of labeled outcomes (length n).
#' @param f_l Vector of predictions on labeled covariates (length n).
#' @param X_u Matrix of unlabeled covariates (N × p).
#' @param f_u Vector of predictions on unlabeled covariates (length N).
#'
#' @return A list with:
#' \item{theta_hat}{Numeric vector of coefficients \eqn{\hat\theta}.}
#' \item{V_theta}{Estimated covariance matrix of \eqn{\hat\theta}.}
#' \item{se}{Standard errors (diagonal of \code{V_theta}).}
#' \item{pieces}{List with \code{V_u}, \code{V_l}, \code{n}, \code{N}.}
#'
#' @details
#' Internals rely on helpers: \code{ols_fit()}, \code{bread_inv()},
#' and \code{meat_matrix()}.
#'
#' @keywords internal
#' @noRd

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
