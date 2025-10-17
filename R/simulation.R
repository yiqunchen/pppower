#' Monte Carlo power for the PP mean estimator
#'
#' @description
#' Simulates two-sided test power under the normal approximation, returning both
#' the empirical Monte Carlo estimate and the analytical power from `power_ppi_mean()` and `power_ppplus_mean()`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size used for imputation.
#' @param n Labeled sample size used for rectification.
#' @param alpha Two-sided significance level (default 0.05).
#' @param R Number of Monte Carlo draws (default 100000).
#' @param var_f Variance of the predictor function \eqn{f(X)}.
#' @param var_res Variance of the residuals \eqn{Y - f(X)}.
#' @param sigma_y2 Optional outcome variance; overrides anything implied by `metrics`.
#' @param sigma_f2 Optional prediction variance; overrides anything implied by `metrics`.
#' @param cov_y_f Optional covariance \eqn{\Cov(Y, f(X))}. When supplied (directly
#'   or via `metrics$cov_y_f`) the PPI++ power is returned in addition to the PP
#'   quantities.
#' @param metrics Optional list of predictive-performance summaries (e.g.,
#'   `mse`, `var_y`, `cov_y_f`, etc.) used to back out missing variance pieces.
#' @param metric_type Character string identifying the metric bundle supplied in
#'   `metrics` (e.g., `"continuous"`, `"hard"`, `"prob"`, `"precision_recall"`).
#' @param m_labeled Sample size associated with the metrics (defaults to `n`);
#'   needed for finite-sample corrections when pulling variances from metrics.
#' @param correction Logical; apply the unbiased finite-sample correction to
#'   the variance estimates (default `TRUE`).
#' @return Named numeric vector of length two with entries `Exact` (analytical power)
#'   and `Empirical` (Monte Carlo estimate).
#' 
#' @export 
simulate_power <- function(delta,
                           N,
                           n,
                           alpha = 0.05,
                           R = 100000,
                           var_f = NULL,
                           var_res = NULL,
                           sigma_y2 = NULL,
                           sigma_f2 = NULL,
                           cov_y_f = NULL,
                           metrics = NULL,
                           metric_type = NULL,
                           m_labeled = n,
                           correction = TRUE) {

  # vanilla PP variance pieces
  comps <- resolve_ppi_variances(
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )

  se_pp <- sqrt(comps$var_f / N + comps$var_res / n)

  # Monte Carlo for empirical power (same Gaussian draw)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu_pp <- abs(delta) / se_pp
  draws <- stats::rnorm(R, mean = mu_pp, sd = 1)
  empirical_pp <- mean(abs(draws) > z_alpha)

  # analytical powers
  exact_pp <- 1 - stats::pnorm(z_alpha - mu_pp) + stats::pnorm(-z_alpha - mu_pp)
  out <- c(Exact_PP = exact_pp, Empirical_PP = empirical_pp)

  cov_y_f_input <- cov_y_f %||% metrics$cov_y_f 
  if (!is.null(cov_y_f_input)) {  # Only add ppi++ if user supply cov_y_f
    ppplus <- resolve_ppi_pp_moments(
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f_input,
      var_f = comps$var_f,
      var_res = comps$var_res,
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = m_labeled,
      correction = correction
    )
    se_ppplus <- sqrt(ppi_pp_variance(
      n = n,
      N = N,
      sigma_y2 = ppplus$sigma_y2,
      sigma_f2 = ppplus$sigma_f2,
      cov_y_f = ppplus$cov_y_f,
      lambda_type = "oracle"
    ))
    mu_ppplus <- abs(delta) / se_ppplus
    draws_ppplus <- stats::rnorm(R, mean = mu_ppplus, sd = 1)
    empirical_ppplus <- mean(abs(draws_ppplus) > z_alpha)
    exact_ppplus <- 1 - stats::pnorm(z_alpha - mu_ppplus) +
                    stats::pnorm(-z_alpha - mu_ppplus)
    out <- c(out, Exact_PPplus = exact_ppplus, Empirical_PPplus = empirical_ppplus)
  }
  out
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
simulate_crossfit_data <- function(n = 20000, p = 5,
                                   family = stats::binomial(),
                                   K = 5, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)

  beta <- seq(0.3, length.out = p, by = 0.5)

  if (identical(family$family, "binomial")) {
    # Logistic DGP
    eta <- drop(X %*% beta)
    prob <- stats::plogis(eta)
    y <- stats::rbinom(n, size = 1, prob = prob)
  } else {
    # Gaussian DGP
    mu <- drop(X %*% beta)
    y <- mu + stats::rnorm(n, sd = 2.0)
  }

  # Cross-fitted predictions (out-of-fold)
  fhat_cf <- crossfit_glm(X, y, K = K, family = family, seed = seed)

  # Return tidy data.frame
  df <- data.frame(y = y, X, fhat_cf = fhat_cf)
  attr(df, "family") <- family$family
  attr(df, "K") <- K
  df
}
