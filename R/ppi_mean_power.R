#' Power for prediction-powered mean estimator
#'
#' @description
#' Computes two-sided test power for the PP mean estimator under a normal
#' approximation. Variance components can be supplied directly or recovered from
#' predictive metrics via `resolve_ppi_variances()`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param alpha Two-sided significance level.
#' @param var_f Optional variance of \eqn{f(X)}.
#' @param var_res Optional residual variance of residual \eqn{Y - f(X)}.
#' @param metrics Optional list of predictive-performance summaries.
#' @param metric_type Character string describing `metrics` (e.g. `"continuous"`, `"hard"`, `"prob"`, `"precision_recall"`).
#' @param m_labeled Labeled sample size associated with `metrics` (defaults to `n`).
#' @param correction Logical; apply finite-sample variance corrections when
#'   recovering moments from metrics.
#'
#' @return Numeric scalar power in \[0, 1\].
#' @export
power_ppi_mean <- function(delta,
                           N,
                           n,
                           alpha = 0.05,
                           var_f = NULL,
                           var_res = NULL,
                           metrics = NULL,
                           metric_type = NULL,
                           m_labeled = n,
                           correction = TRUE) {
  comps <- resolve_ppi_variances(
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )
  se <- sqrt(comps$var_f / N + comps$var_res / n)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}

#' @keywords Internal
#' Monte Carlo power using the df (superpopulation view)
simulate_power_ppi_mean <- function(df, N, n, alpha = 0.05,
                              R = 2000, theta0 = NULL,
                              use_sample_var = TRUE, seed = 1) {
  set.seed(seed)
  stopifnot(all(c("y", "fhat_cf") %in% names(df)))

  # "True" theta implied by df (superpopulation proxy)
  theta <- mean(df$y)
  if (is.null(theta0)) stop("Please provide theta0 (null).")
  delta <- theta - theta0

  rej <- logical(R)
  thetas <- numeric(R); ses <- numeric(R)

  # Plug-in "population" variances for analytical power (from full df)
  var_f_pop   <- stats::var(df$fhat_cf)
  var_res_pop <- stats::var(df$y - df$fhat_cf)
  power_exact <- power_ppi_mean(
    delta = delta,
    N = N,
    n = n,
    alpha = alpha,
    var_f = var_f_pop,
    var_res = var_res_pop
  )

  for (r in seq_len(R)) {
    out       <- pp_once(df, N, n, use_sample_var = use_sample_var, seed = sample.int(1e9,1))
    z         <- (out$theta_hat - theta0) / out$se
    rej[r]    <- (abs(z) > stats::qnorm(1 - alpha/2))
    thetas[r] <- out$theta_hat
    ses[r]    <- out$se
  }

  list(
    theta       = theta,
    theta0      = theta0,
    delta       = delta,
    empirical_power = mean(rej),
    analytical_power = power_exact,
    avg_SE      = mean(ses),
    details     = data.frame(theta_hat = thetas, se = ses, reject = rej)
  )
}
