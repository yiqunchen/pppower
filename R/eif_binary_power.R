#' Resolve EIF inputs for binary surrogates
#'
#' Accepts either the trio `(mu0, mu1, p_hat)` or a measurement-error
#' parameterization via `(p, sens, spec)` and returns the quantities needed
#' by EIF power/sample-size routines.
#'
#' @param mu0,mu1 Conditional means \eqn{E[Y \mid \hat Y = 0]} and
#'   \eqn{E[Y \mid \hat Y = 1]}.
#' @param p_hat Prevalence of predicted positives \eqn{P(\hat Y = 1)}.
#' @param p True prevalence \eqn{P(Y = 1)}.
#' @param sens,spec Sensitivity and specificity of the surrogate classifier.
#'
#' @return List with `mu0`, `mu1`, `p_hat`, and `theta` (prevalence).
#' @keywords internal
resolve_eif_binary_inputs <- function(mu0 = NULL,
                                      mu1 = NULL,
                                      p_hat = NULL,
                                      p = NULL,
                                      sens = NULL,
                                      spec = NULL) {
  has_direct <- !is.null(mu0) && !is.null(mu1) && !is.null(p_hat)
  has_meas_error <- !is.null(p) && !is.null(sens) && !is.null(spec)

  if (!has_direct && !has_meas_error) {
    stop("Supply either (mu0, mu1, p_hat) or (p, sens, spec).", call. = FALSE)
  }

  if (has_direct) {
    vals <- c(mu0, mu1, p_hat)
    if (any(!is.numeric(vals)) || any(vals < 0) || any(vals > 1)) {
      stop("mu0, mu1, and p_hat must lie in [0, 1].", call. = FALSE)
    }
    theta <- (1 - p_hat) * mu0 + p_hat * mu1
    return(list(mu0 = mu0, mu1 = mu1, p_hat = p_hat, theta = theta))
  }

  # Measurement-error path: derive from prevalence and sens/spec
  moments <- binary_moments_from_sens_spec(p = p, sens = sens, spec = spec)
  if (moments$p_hat <= 0 || moments$p_hat >= 1) {
    stop("Derived p_hat from sens/spec is degenerate; adjust inputs.", call. = FALSE)
  }

  mu1 <- (sens * p) / moments$p_hat
  mu0 <- ((1 - sens) * p) / (1 - moments$p_hat)

  list(
    mu0 = mu0,
    mu1 = mu1,
    p_hat = moments$p_hat,
    theta = p
  )
}

#' Power or required calibration size for EIF estimator with binary surrogates
#'
#' @description
#' Unified interface for the efficient influence function (EIF) estimator used
#' with binary surrogate predictions (e.g., LLM-as-a-judge hard votes).
#' Follows the \pkg{pwr} convention: supply **either** `m_cal` or `power`
#' (but not both). The function computes whichever is left `NULL`.
#'
#' You can provide either direct conditional means `(mu0, mu1, p_hat)` or
#' sensitivity/specificity plus prevalence `(p, sens, spec)`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled/test sample size with surrogate predictions.
#' @param m_cal Calibration (labeled) sample size. Set to `NULL` to solve for
#'   the required `m_cal`.
#' @param power Target power. Set to `NULL` (default) to compute power from
#'   `m_cal`.
#' @param alpha Two-sided significance level (default 0.05).
#' @inheritParams resolve_eif_binary_inputs
#'
#' @return When computing power: scalar in \[0, 1\] with attributes `theta`,
#'   `p_hat`, `var_unlabeled`, and `var_cal`.
#'   When computing sample size: integer required calibration size.
#' @export
power_eif_binary <- function(delta,
                             N,
                             m_cal = NULL,
                             power = NULL,
                             alpha = 0.05,
                             mu0 = NULL,
                             mu1 = NULL,
                             p_hat = NULL,
                             p = NULL,
                             sens = NULL,
                             spec = NULL) {
  # Validate: exactly one of m_cal, power must be NULL
  if (is.null(m_cal) && is.null(power)) {
    stop("Exactly one of 'm_cal' and 'power' must be NULL.", call. = FALSE)
  }
  if (!is.null(m_cal) && !is.null(power)) {
    stop("Exactly one of 'm_cal' and 'power' must be NULL.", call. = FALSE)
  }

  if (!is.numeric(delta) || length(delta) != 1L) {
    stop("delta must be a numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || N <= 0) stop("N must be positive.", call. = FALSE)

  inp <- resolve_eif_binary_inputs(mu0 = mu0, mu1 = mu1, p_hat = p_hat,
                                   p = p, sens = sens, spec = spec)

  if (!is.null(m_cal)) {
    # ---- Compute power given m_cal ----
    if (!is.numeric(m_cal) || m_cal <= 0) stop("m_cal must be positive.", call. = FALSE)

    var_unlabeled <- (inp$mu1 - inp$mu0)^2 * inp$p_hat * (1 - inp$p_hat) / N
    var_cal <- ((1 - inp$p_hat) * inp$mu0 * (1 - inp$mu0) +
                  inp$p_hat * inp$mu1 * (1 - inp$mu1)) / m_cal
    se <- sqrt(var_unlabeled + var_cal)

    z_alpha <- stats::qnorm(1 - alpha / 2)
    if (se == 0) {
      pwr <- if (delta == 0) alpha else 1
    } else {
      mu <- abs(delta) / se
      pwr <- 1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
    }

    attr(pwr, "theta") <- inp$theta
    attr(pwr, "p_hat") <- inp$p_hat
    attr(pwr, "var_unlabeled") <- var_unlabeled
    attr(pwr, "var_cal") <- var_cal
    return(pwr)

  } else {
    # ---- Solve for m_cal given target power ----
    if (power <= 0 || power >= 1) stop("power must lie in (0, 1).", call. = FALSE)
    z_alpha <- stats::qnorm(1 - alpha / 2)
    z_beta <- stats::qnorm(power)
    S2 <- (delta / (z_alpha + z_beta))^2

    var_unlabeled <- (inp$mu1 - inp$mu0)^2 * inp$p_hat * (1 - inp$p_hat) / N

    if (var_unlabeled >= S2) {
      stop("Infeasible: unlabeled variance alone exceeds target variance.", call. = FALSE)
    }

    numer <- (1 - inp$p_hat) * inp$mu0 * (1 - inp$mu0) +
      inp$p_hat * inp$mu1 * (1 - inp$mu1)
    denom <- S2 - var_unlabeled
    m_needed <- ceiling(numer / denom)
    return(as.integer(m_needed))
  }
}

#' Monte Carlo power for the EIF binary estimator
#'
#' @param R Number of Monte Carlo replicates.
#' @param m_cal Calibration sample size.
#' @param N Unlabeled/test sample size with surrogate predictions.
#' @param alpha Two-sided significance level.
#' @param seed RNG seed for reproducibility.
#' @inheritParams power_eif_binary
#' @return List with empirical and analytical power plus simulation details.
#' @export
simulate_eif_binary <- function(R,
                                delta,
                                N,
                                m_cal,
                                alpha = 0.05,
                                mu0 = NULL,
                                mu1 = NULL,
                                p_hat = NULL,
                                p = NULL,
                                sens = NULL,
                                spec = NULL,
                                seed = 1) {
  if (!is.null(seed)) set.seed(seed)
  inp <- resolve_eif_binary_inputs(mu0 = mu0, mu1 = mu1, p_hat = p_hat,
                                   p = p, sens = sens, spec = spec)

  theta_true <- inp$theta
  theta0 <- theta_true - delta
  analytical <- power_eif_binary(
    delta = delta,
    N = N,
    m_cal = m_cal,
    alpha = alpha,
    mu0 = inp$mu0,
    mu1 = inp$mu1,
    p_hat = inp$p_hat
  )

  seeds <- sample.int(.Machine$integer.max, R)
  rejects <- vapply(
    seq_len(R),
    function(i) {
      if (!is.null(seed)) set.seed(seeds[i])
      # Calibration draws
      yhat_cal <- stats::rbinom(m_cal, 1, inp$p_hat)
      y_cal <- ifelse(
        yhat_cal == 1,
        stats::rbinom(m_cal, 1, inp$mu1),
        stats::rbinom(m_cal, 1, inp$mu0)
      )
      m1 <- sum(yhat_cal == 1)
      m0 <- m_cal - m1
      if (m0 == 0 || m1 == 0) return(FALSE)

      mu1_hat <- sum(y_cal[yhat_cal == 1]) / m1
      mu0_hat <- sum(y_cal[yhat_cal == 0]) / m0

      # Unlabeled predictions
      yhat_test <- stats::rbinom(N, 1, inp$p_hat)
      p_hat_test <- mean(yhat_test)

      theta_hat <- (1 - p_hat_test) * mu0_hat + p_hat_test * mu1_hat

      var_hat <- (mu1_hat - mu0_hat)^2 * p_hat_test * (1 - p_hat_test) / N +
        (1 - p_hat_test)^2 * mu0_hat * (1 - mu0_hat) / m0 +
        p_hat_test^2 * mu1_hat * (1 - mu1_hat) / m1

      se_hat <- sqrt(max(var_hat, 0))
      if (se_hat == 0) return(FALSE)
      z_stat <- (theta_hat - theta0) / se_hat
      abs(z_stat) > stats::qnorm(1 - alpha / 2)
    },
    logical(1)
  )

  empirical <- mean(rejects)
  mc_se <- sqrt(analytical * (1 - analytical) / R)

  list(
    empirical_power = empirical,
    analytical_power = analytical,
    mc_se = mc_se,
    theta_true = theta_true,
    theta0 = theta0,
    delta_used = delta,
    inputs = inp,
    R = R
  )
}


# ---- Deprecated shims ----

#' @rdname deprecated
#' @export
n_required_eif_binary <- function(delta, N, power = 0.80, alpha = 0.05,
                                  mu0 = NULL, mu1 = NULL, p_hat = NULL,
                                  p = NULL, sens = NULL, spec = NULL) {
  .Deprecated("power_eif_binary", package = "pppower",
              msg = "n_required_eif_binary() is deprecated. Use power_eif_binary(..., m_cal = NULL, power = 0.80) instead.")
  power_eif_binary(delta = delta, N = N, m_cal = NULL, power = power,
                   alpha = alpha, mu0 = mu0, mu1 = mu1, p_hat = p_hat,
                   p = p, sens = sens, spec = spec)
}

#' @rdname deprecated
#' @export
simulate_power_eif_binary <- function(...) {
  .Deprecated("simulate_eif_binary", package = "pppower")
  simulate_eif_binary(...)
}
