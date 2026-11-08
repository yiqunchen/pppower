#' PPI++ One-Sample Mean Test
#'
#' @description
#' Performs a one-sample hypothesis test for the population mean using the
#' PPI++ (Prediction-Powered Inference Plus) estimator with EIF-optimal lambda.
#'
#' @param Y_L Numeric vector of labeled outcomes.
#' @param f_L Numeric vector of predictions for labeled samples.
#' @param f_U Numeric vector of predictions for unlabeled samples.
#' @param theta0 Null hypothesis value (default 0).
#' @param alpha Significance level (default 0.05).
#' @param lambda Optional user-specified lambda. If NULL, uses EIF-optimal lambda.
#' @param se_type `"standard"` (default) uses the usual plug-in standard error
#'   and a normal reference. `"corrected"` adds the second-order term
#'   \eqn{c\,\hat\sigma_\varepsilon^2\hat\kappa/\{n(n-3)\}} for the estimated
#'   \eqn{\lambda}, with \eqn{\hat\kappa} bounded below by 1, and uses a
#'   \eqn{t_{n-2}} reference. It is meant for small `n`; it is ignored when
#'   `lambda` is supplied, because a fixed \eqn{\lambda} needs no correction.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{The PPI++ point estimate of the mean.}
#'   \item{se}{Standard error of the estimate.}
#'   \item{z_stat}{Test statistic (compared with a \eqn{t_{n-2}} distribution
#'     when `se_type = "corrected"`).}
#'   \item{df}{Degrees of freedom of the reference distribution (`Inf` for normal).}
#'   \item{p_value}{Two-sided p-value.}
#'   \item{reject}{Logical indicating whether null is rejected at level alpha.}
#'   \item{lambda}{The lambda value used.}
#'   \item{ci}{Confidence interval at level 1-alpha.}
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 100; N <- 1000
#' Y_L <- rnorm(n, mean = 0.5)
#' f_L <- Y_L + rnorm(n, sd = 0.5)
#' f_U <- rnorm(N, mean = 0.5) + rnorm(N, sd = 0.5)
#' ppi_mean_test(Y_L, f_L, f_U, theta0 = 0)
#'
#' @export
ppi_mean_test <- function(Y_L, f_L, f_U, theta0 = 0, alpha = 0.05, lambda = NULL,
                          se_type = c("standard", "corrected")) {
  se_type <- match.arg(se_type)
  estimated_lambda <- is.null(lambda)
  n <- length(Y_L)
  N <- length(f_U)

  # Sample means
  Y_bar <- mean(Y_L)
  f_L_bar <- mean(f_L)
  f_U_bar <- mean(f_U)

  # Estimate variances and covariance

  sigma_Y2_hat <- stats::var(Y_L)
  sigma_f2_hat <- stats::var(f_L)
  cov_Yf_hat <- stats::cov(Y_L, f_L)

  # Use provided lambda or compute EIF-optimal
  if (is.null(lambda)) {
    lambda <- plugin_lambda(cov_Yf_hat, sigma_f2_hat, n, N)
  }

  # PPI++ estimator
  theta_hat <- Y_bar + lambda * (f_U_bar - f_L_bar)

  # Variance estimate
  var_hat <- sigma_Y2_hat / n +
    lambda^2 * sigma_f2_hat * (1 / N + 1 / n) -
    2 * lambda * cov_Yf_hat / n
  var_hat <- max(var_hat, 1e-10)
  df <- Inf
  if (se_type == "corrected" && estimated_lambda) {
    if (n <= 3) stop("se_type = \"corrected\" needs n > 3.", call. = FALSE)
    s_eps2 <- max(sigma_Y2_hat - cov_Yf_hat^2 / sigma_f2_hat, 0) * (n - 1) / (n - 2)
    var_hat <- var_hat + N / (n + N) * s_eps2 * kappa_hat(Y_L, f_L) / (n * (n - 3))
    df <- n - 2
  }
  se <- sqrt(var_hat)

  # Test statistic and p-value
  z_stat <- (theta_hat - theta0) / se
  p_value <- 2 * stats::pt(-abs(z_stat), df)
  reject <- p_value < alpha


  # Confidence interval
  z_alpha <- stats::qt(1 - alpha / 2, df)
  ci <- c(theta_hat - z_alpha * se, theta_hat + z_alpha * se)

  list(
    estimate = theta_hat,
    se = se,
    z_stat = z_stat,
    df = df,
    p_value = p_value,
    reject = reject,
    lambda = lambda,
    ci = ci
  )
}

#' PPI++ Two-Sample t-Test
#'
#' @description
#' Performs a two-sample hypothesis test for the difference in means between
#' two independent groups using PPI++ estimators.
#'
#' @param Y_A Numeric vector of labeled outcomes for group A.
#' @param f_A_L Numeric vector of predictions for labeled samples in group A.
#' @param f_A_U Numeric vector of predictions for unlabeled samples in group A.
#' @param Y_B Numeric vector of labeled outcomes for group B.
#' @param f_B_L Numeric vector of predictions for labeled samples in group B.
#' @param f_B_U Numeric vector of predictions for unlabeled samples in group B.
#' @param delta0 Null hypothesis value for the difference (default 0).
#' @param alpha Significance level (default 0.05).
#' @param lambda_A Optional lambda for group A. If NULL, uses EIF-optimal.
#' @param lambda_B Optional lambda for group B. If NULL, uses EIF-optimal.
#'
#' @return A list with components:
#' \describe{
#'   \item{estimate}{The PPI++ estimate of the difference (mean_A - mean_B).}
#'   \item{se}{Standard error of the difference estimate.}
#'   \item{z_stat}{Z test statistic.}
#'   \item{p_value}{Two-sided p-value.}
#'   \item{reject}{Logical indicating whether null is rejected at level alpha.}
#'   \item{lambda_A}{Lambda used for group A.}
#'   \item{lambda_B}{Lambda used for group B.}
#'   \item{ci}{Confidence interval at level 1-alpha.}
#' }
#'
#' @examples
#' set.seed(123)
#' n <- 50; N <- 500
#' # Group A (treatment)
#' Y_A <- rnorm(n, mean = 0.3)
#' f_A_L <- Y_A + rnorm(n, sd = 0.5)
#' f_A_U <- rnorm(N, mean = 0.3) + rnorm(N, sd = 0.5)
#' # Group B (control)
#' Y_B <- rnorm(n, mean = 0)
#' f_B_L <- Y_B + rnorm(n, sd = 0.5)
#' f_B_U <- rnorm(N, mean = 0) + rnorm(N, sd = 0.5)
#'
#' ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U)
#'
#' @export
ppi_ttest <- function(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U,
                      delta0 = 0, alpha = 0.05,
                      lambda_A = NULL, lambda_B = NULL) {
  # Group A
  n_A <- length(Y_A)
  N_A <- length(f_A_U)
  Y_A_bar <- mean(Y_A)
  f_A_L_bar <- mean(f_A_L)
  f_A_U_bar <- mean(f_A_U)
  sigma_YA2 <- stats::var(Y_A)
  sigma_fA2 <- stats::var(f_A_L)
  cov_A <- stats::cov(Y_A, f_A_L)

  if (is.null(lambda_A)) {
    lambda_A <- plugin_lambda(cov_A, sigma_fA2, n_A, N_A)
  }
  theta_A <- Y_A_bar + lambda_A * (f_A_U_bar - f_A_L_bar)

  # Group B
  n_B <- length(Y_B)
  N_B <- length(f_B_U)
  Y_B_bar <- mean(Y_B)
  f_B_L_bar <- mean(f_B_L)
  f_B_U_bar <- mean(f_B_U)
  sigma_YB2 <- stats::var(Y_B)
  sigma_fB2 <- stats::var(f_B_L)
  cov_B <- stats::cov(Y_B, f_B_L)

  if (is.null(lambda_B)) {
    lambda_B <- plugin_lambda(cov_B, sigma_fB2, n_B, N_B)
  }
  theta_B <- Y_B_bar + lambda_B * (f_B_U_bar - f_B_L_bar)

  # Difference estimate
  delta_hat <- theta_A - theta_B

  # Variance of difference
  var_A <- sigma_YA2 / n_A + lambda_A^2 * sigma_fA2 * (1 / N_A + 1 / n_A) -
    2 * lambda_A * cov_A / n_A
  var_B <- sigma_YB2 / n_B + lambda_B^2 * sigma_fB2 * (1 / N_B + 1 / n_B) -
    2 * lambda_B * cov_B / n_B
  var_hat <- max(var_A, 0) + max(var_B, 0)
  var_hat <- max(var_hat, 1e-10)
  se <- sqrt(var_hat)

  # Test statistic and p-value
  z_stat <- (delta_hat - delta0) / se
  p_value <- 2 * stats::pnorm(-abs(z_stat))
  reject <- p_value < alpha

  # Confidence interval
  z_alpha <- stats::qnorm(1 - alpha / 2)
  ci <- c(delta_hat - z_alpha * se, delta_hat + z_alpha * se)

  list(
    estimate = delta_hat,
    se = se,
    z_stat = z_stat,
    p_value = p_value,
    reject = reject,
    lambda_A = lambda_A,
    lambda_B = lambda_B,
    ci = ci
  )
}

#' PPI++ Paired t-Test
#'
#' @description
#' Performs a paired hypothesis test for the mean difference using PPI++.
#' This is equivalent to a one-sample test on the differences D = Y_A - Y_B.
#'
#' @param D_L Numeric vector of labeled differences (Y_A - Y_B).
#' @param fD_L Numeric vector of predicted differences for labeled samples.
#' @param fD_U Numeric vector of predicted differences for unlabeled samples.
#' @param delta0 Null hypothesis value for the mean difference (default 0).
#' @param alpha Significance level (default 0.05).
#' @param lambda Optional lambda. If NULL, uses EIF-optimal.
#'
#' @return A list with components (same as \code{ppi_mean_test}).
#'
#' @examples
#' set.seed(123)
#' n <- 50; N <- 500
#' D_L <- rnorm(n, mean = 0.3)  # True difference
#' fD_L <- D_L + rnorm(n, sd = 0.5)
#' fD_U <- rnorm(N, mean = 0.3) + rnorm(N, sd = 0.5)
#' ppi_paired_test(D_L, fD_L, fD_U)
#'
#' @export
ppi_paired_test <- function(D_L, fD_L, fD_U, delta0 = 0, alpha = 0.05, lambda = NULL) {
  # Paired test is just a one-sample test on differences
  ppi_mean_test(Y_L = D_L, f_L = fD_L, f_U = fD_U,
                theta0 = delta0, alpha = alpha, lambda = lambda)
}

# Sample estimate of the moment factor kappa for the estimated tuning parameter,
# bounded to [1, 50] (the formula assumes a linear E[Y | f]).
# Plug-in lambda = c * slope; the slope is undefined when all labeled
# predictions are equal (e.g., a binary f with one class), so use lambda = 0.
plugin_lambda <- function(cov_yf, sigma_f2, n, N) {
  if (!is.finite(sigma_f2) || sigma_f2 <= 0) return(0)
  cov_yf / ((1 + n / N) * sigma_f2)
}

kappa_hat <- function(Y, f) {
  fc <- f - mean(f)
  sf2 <- mean(fc^2)
  if (sf2 <= 0) return(1)
  b <- sum((Y - mean(Y)) * fc) / sum(fc^2)
  e <- Y - mean(Y) - b * fc
  se2 <- mean(e^2)
  if (se2 <= 0) return(1)
  k <- 2 - mean(fc^2 * e^2) / (sf2 * se2) + 2 * mean(fc * (e^2 - se2)) * mean(fc^3) / (sf2^2 * se2)
  if (!is.finite(k)) return(1)
  min(max(k, 1), 50)
}
