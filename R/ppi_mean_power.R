#' Power for prediction-powered mean estimator
#'
#' @description 
#' Computes two-sided test power for the PP mean estimator under
#' normal approximation.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param var_f Variance of \eqn{f(X)}.
#' @param var_res Variance of residuals \eqn{Y - f(X)}.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param alpha Two-sided significance level.
#'
#' @return Numeric scalar power in \[0, 1\].
#'
#' @details
#' Standard error is \eqn{\sqrt{\mathrm{Var}(f)/N + \mathrm{Var}(Y-f)/n}}.
#'
#' @examples
#' power_ppi_mean(delta = 0.2, var_f = 0.4, var_res = 1.0, N = 1000, n = 200)
#'
#' @export
power_ppi_mean <- function(delta, var_f, var_res, N, n, alpha = 0.05) {
  se <- sqrt(var_f / N + var_res / n)
  z_alpha <- qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  # Exact two-sided normal-theory power
  power <- 1 - pnorm(z_alpha - mu) + pnorm(-z_alpha - mu)
  power
}

# Monte Carlo power using the df (superpopulation view)
# theta0: null; delta = theta - theta0 (if you prefer to pass delta directly)
simulate_power_PP <- function(df, N, n, alpha = 0.05,
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
  power_exact <- power_ppi_mean(delta, var_f_pop, var_res_pop, N, n, alpha)

  for (r in 1:R) {
    out <- pp_once(df, N, n, use_sample_var = use_sample_var, seed = sample.int(1e9,1))
    z   <- (out$theta_hat - theta0) / out$se
    rej[r] <- (abs(z) > stats::qnorm(1 - alpha/2))
    thetas[r] <- out$theta_hat; ses[r] <- out$se
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
