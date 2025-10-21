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
#' @param metric_type Character string describing `metrics` (e.g. `"continuous"`, `"prob"`, `"classification"`).
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

#' Monte Carlo power for the PPI mean estimator
#'
#' @param R Integer, number of Monte Carlo draws.
#' @param n Labeled sample size.
#' @param N Unlabeled sample size.
#' @param alpha Two-sided test level.
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param var_f Variance of the predictor \eqn{f(X)}.
#' @param var_res Variance of the residual \eqn{Y - f(X)}.
#' @param moments Optional list with entries `delta`, `var_f`, and `var_res`.
#' @param seed RNG seed (set to `NULL` to skip seeding).
#' @param family A `stats::family()` object (`gaussian()` or `binomial()`), default is `gaussian()`.
#' @param theta0 Null hypothesis value \eqn{\theta_0}.
#'
#' @return A list with empirical and analytical power, average standard error, and simulation details.
#' @export
simulate_power_ppi_mean <- function(R,
                                    n,
                                    N,
                                    alpha = 0.05,
                                    delta = NULL,
                                    var_f = NULL,
                                    var_res = NULL,
                                    moments = NULL,
                                    seed = 1,
                                    family = stats::gaussian(),
                                    theta0 = 0) {
  if (!is.null(seed)) set.seed(seed)

  # ---- Resolve population moments ----
  if (!is.null(moments)) {
    delta   <- moments$delta   %||% delta
    var_f   <- moments$var_f   %||% var_f
    var_res <- moments$var_res %||% var_res
  }

  if (is.null(delta) || is.null(var_f)) {
    stop("simulate_power_ppi_mean requires `delta` and `var_f`.", call. = FALSE)
  }

  theta_true <- theta0 + delta

  if (identical(family$family, "binomial")) {
    if (theta_true <= 0 || theta_true >= 1) {
      stop("For binomial mean simulation, theta_true must lie in (0,1).", call. = FALSE)
    }
    var_res <- var_res %||% (theta_true * (1 - theta_true) - var_f)
  } else if (identical(family$family, "gaussian")) {
    if (is.null(var_res)) stop("Gaussian mean simulation requires `var_res`.", call. = FALSE)
  } else {
    stop("Unsupported family: must be gaussian() or binomial().", call. = FALSE)
  }

  # Analytical power
  se_pop <- sqrt(var_f / N + var_res / n)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  if (se_pop == 0) {
    analytical_power <- if (delta == 0) alpha else 1
  } else {
    mu_pop <- abs(delta) / se_pop
    analytical_power <- 1 - stats::pnorm(z_alpha - mu_pop) + stats::pnorm(-z_alpha - mu_pop)
  }

  # Monte Carlo simulation
  seeds <- sample.int(.Machine$integer.max, R)
  reps <- lapply(
    seq_len(R),
    function(i) simulate_ppi_mean_rep(
      n_l = n,
      n_u = N,
      family = family,
      method = "ppi",
      alpha = alpha,
      theta0 = theta0,
      delta = delta,
      var_f = var_f,
      var_res = var_res,
      seed = seeds[i]
    )
  )

  empirical_power <- mean(vapply(reps, `[[`, logical(1), "reject"))
  avg_SE <- mean(vapply(reps, function(x) x$se, numeric(1)))
  mean_theta <- mean(vapply(reps, function(x) x$theta_hat, numeric(1)))

  list(
    empirical_power = empirical_power,
    analytical_power = analytical_power,
    avg_SE = avg_SE,
    delta_used = delta,
    avg_theta_hat = mean_theta,
    details = reps
  )
}