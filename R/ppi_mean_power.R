# Internal: vanilla PPI power (lambda=1, no covariance)
# Previously exported as power_ppi_mean(); now superseded by the unified
# power_ppi_mean() in ppi_plus_mean_power.R.
#
# @keywords internal
power_ppi_vanilla_mean <- function(delta,
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

# Internal: vanilla PPI Monte Carlo simulation
# Previously exported as simulate_power_ppi_mean().
#
# @keywords internal
simulate_ppi_vanilla_mean <- function(R,
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

  # Resolve population moments
  if (!is.null(moments)) {
    delta   <- moments$delta   %||% delta
    var_f   <- moments$var_f   %||% var_f
    var_res <- moments$var_res %||% var_res
  }

  if (is.null(delta) || is.null(var_f)) {
    stop("simulate_ppi_vanilla_mean requires `delta` and `var_f`.", call. = FALSE)
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

# ---- Deprecated shims (old names) ----

#' @rdname deprecated
#' @export
simulate_power_ppi_mean <- function(...) {
  .Deprecated("simulate_ppi_vanilla_mean", package = "pppower",
              msg = "simulate_power_ppi_mean() is deprecated. Use simulate_ppi_vanilla_mean() for vanilla PPI or simulate_ppi_mean() for PPI++.")
  simulate_ppi_vanilla_mean(...)
}
