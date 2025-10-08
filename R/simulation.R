#' Monte Carlo power for the PP mean estimator
#'
#' @description
#' Simulates two-sided test power under the normal approximation, returning both
#' the empirical Monte Carlo estimate and the analytical power from `power_ppi_mean()`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param var_f Variance of the predictor function \eqn{f(X)}.
#' @param var_res Variance of the residuals \eqn{Y - f(X)}.
#' @param N Unlabeled sample size used for imputation.
#' @param n Labeled sample size used for rectification.
#' @param alpha Two-sided significance level (default 0.05).
#' @param R Number of Monte Carlo draws (default 100000).
#'
#' @return Named numeric vector of length two with entries `Exact` (analytical power)
#'   and `Empirical` (Monte Carlo estimate).
#' 
#' @export 
simulate_power <- function(delta, var_f, var_res, N, n, alpha = 0.05, R = 100000) {
  se <- sqrt(var_f / N + var_res / n)
  mu <- abs(delta) / se
  Z <- rnorm(R, mean = mu, sd = 1)
  c_empirical = mean(abs(Z) > qnorm(1 - alpha / 2))
  c_exact     = power_ppi_mean(delta, var_f, var_res, N, n, alpha)
  c(Exact = c_exact, Empirical = c_empirical)
}

#' Simulate cross-fitted predictive data
#'
#' @description
#' Generates a synthetic regression/classification dataset with covariates, outcomes,
#' and out-of-fold predictions from a K-fold cross-fitted GLM.
#'
#' @param n Number of observations to simulate.
#' @param p Number of covariates (columns in the design matrix).
#' @param family GLM family (defaults to `stats::binomial()`; `stats::gaussian()` gives a Gaussian response).
#' @param K Number of folds used for cross-fitting.
#' @param seed RNG seed for reproducibility.
#'
#' @return A data frame containing the outcome `y`, covariate columns `x1, ..., xp`, and
#' `fhat_cf`, the cross-fitted predictions. Attributes `family` and `K` record the generating family and fold count.
#' 
#' @export
simulate_crossfit_data <- function(n = 2000, p = 5,
                                   family = stats::binomial(),
                                   K = 5, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)

  beta <- seq(0.3, length.out = p, by = 0.5)

  if (identical(family$family, "binomial")) {
    # Logistic DGP
    eta <- drop(X %*% beta)
    prob <- plogis(eta)
    y <- rbinom(n, size = 1, prob = prob)
  } else {
    # Gaussian DGP
    mu <- drop(X %*% beta)
    y <- mu + rnorm(n, sd = 2.0)
  }

  # Cross-fitted predictions (out-of-fold)
  fhat_cf <- crossfit_glm(X, y, K = K, family = family, seed = seed)

  # Return tidy data.frame
  df <- data.frame(y = y, X, fhat_cf = fhat_cf)
  attr(df, "family") <- family$family
  attr(df, "K") <- K
  df
}
