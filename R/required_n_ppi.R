#' Required labeled sample size for desired PPI power
#'
#' @description
#' Computes the minimum labeled sample size \eqn{n} that attains a target
#' two-sided power for the vanilla prediction-powered mean estimator, a PP-OLS
#' contrast, or a custom variance pair. For the mean case you can supply
#' \eqn{\Var(f)} and \eqn{\Var(Y-f)} directly or recover them from predictive
#' metrics via `resolve_ppi_variances()`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param alpha Two-sided significance level.
#' @param power Target power (\eqn{1-\beta}).
#' @param type `"mean"`, `"ols"`, or `"custom"`.
#' @param var_f Optional variance components of \eqn{f(X)} 
#' (used when `type = "mean"` and `metrics` is `NULL`).
#' @param var_res Optional variance components of \eqn{Y-f(X)} 
#' (used when `type = "mean"` and `metrics` is `NULL`).
#' @param metrics Optional list of predictive-performance summaries (passed to
#'   `resolve_ppi_variances()`).
#' @param metric_type Character string describing `metrics` (e.g. `"continuous"`,
#'   `"hard"`, `"prob"`, `"precision_recall"`).
#' @param m_labeled Labeled sample size associated with `metrics`; defaults to
#'   `metrics$m_obs` when present.
#' @param correction Logical; apply finite-sample variance corrections when
#'   deriving moments from metrics (default `TRUE`).
#' @param V_u,V_l Sandwich covariance matrices for PP-OLS (required when
#'   `type = "ols"`).
#' @param c Numeric contrast vector for PP-OLS (same length as the regressors).
#' @param var_unlabeled,var_labeled Optional variance components for
#'   `type = "custom"` corresponding to \eqn{\Var_u} and \eqn{\Var_l}.
#' @param warn_smallN Logical; warn if `N < smallN_threshold`.
#' @param smallN_threshold Threshold at which the small-\eqn{N} warning fires.
#' @param mode `"error"` (default) to error when the required \eqn{n} exceeds `N`,
#'   or `"cap"` to return `n = N` with the achieved power attached as an attribute.
#'
#' @return Integer labeled sample size meeting the target power; when `mode = "cap"`
#'   and the cap is hit, the result carries an `"achieved_power"` attribute.
#'
#' @examples
#' # Mean estimator
#' n_required_pp(
#'   delta = 0.1, N = 5000, var_f = 0.2, var_res = 0.5,
#'   alpha = 0.05, power = 0.8, type = "mean"
#' )
#'
#' # PP-OLS contrast
#' V_u <- diag(c(0.8, 0.6)); V_l <- diag(c(0.5, 0.7))
#' n_required_pp(
#'   delta = 0.15, N = 4000, V_u = V_u, V_l = V_l, c = c(1, 0),
#'   alpha = 0.05, power = 0.85, type = "ols"
#' )
#' @export
n_required_pp <- function(delta,
                          N,
                          alpha = 0.05,
                          power = 0.80,
                          type = c("mean", "ols", "custom"),
                          var_f = NULL,
                          var_res = NULL,
                          metrics = NULL,
                          metric_type = NULL,
                          m_labeled = NULL,
                          correction = TRUE,
                          V_u = NULL,
                          V_l = NULL,
                          c = NULL,
                          var_unlabeled = NULL,
                          var_labeled = NULL,
                          warn_smallN = TRUE,
                          smallN_threshold = 500,
                          mode = c("error", "cap")) {
  type <- match.arg(type)
  mode <- match.arg(mode)

  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta == 0) {
    stop("delta must be a non-zero finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || N <= 0) {
    stop("N must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(power) || power <= 0 || power >= 1) {
    stop("power must lie in (0, 1).", call. = FALSE)
  }

  if (warn_smallN && N < smallN_threshold) {
    warning(
      "Unlabeled N is quite small. The finite-N term may dominate the variance.",
      call. = FALSE
    )
  }

  components <- switch(
    type,
    mean = {
      if (is.null(metrics)) {
        if (is.null(var_f)) {
          stop("var_f must be supplied when metrics are not provided.", call. = FALSE)
        }
        if (is.null(var_res)) {
          stop("var_res must be supplied when metrics are not provided.", call. = FALSE)
        }
        if (!is.numeric(var_f) || length(var_f) != 1L || var_f < 0) {
          stop("var_f must be >= 0 when supplied directly.", call. = FALSE)
        }
        if (!is.numeric(var_res) || length(var_res) != 1L || var_res <= 0) {
          stop("var_res must be > 0 when supplied directly.", call. = FALSE)
        }
      }
      comps <- resolve_ppi_variances(
        var_f = var_f,
        var_res = var_res,
        metrics = metrics,
        metric_type = metric_type,
        m_labeled = m_labeled,
        correction = correction
      )
      list(
        var_u = comps$var_f,
        var_l = comps$var_res,
        power_fn = function(n) power_ppi_mean(
          delta = delta,
          var_f = comps$var_f,
          var_res = comps$var_res,
          N = N,
          n = n,
          alpha = alpha
        )
      )
    },
    ols = {
      if (is.null(V_u) || is.null(V_l) || is.null(c)) {
        stop("V_u, V_l, and c must be supplied when type = 'ols'.", call. = FALSE)
      }
      if (!is.matrix(V_u) || !is.matrix(V_l)) {
        stop("V_u and V_l must be matrices.", call. = FALSE)
      }
      c_vec <- as.numeric(c)
      p <- length(c_vec)
      if (p == 0L || anyNA(c_vec)) {
        stop("Contrast c must be a non-empty numeric vector without NA.", call. = FALSE)
      }
      if (!all(dim(V_u) == p) || !all(dim(V_l) == p)) {
        stop("Dimensions of V_u/V_l must match length of c.", call. = FALSE)
      }
      var_u <- as.numeric(drop(t(c_vec) %*% V_u %*% c_vec))
      var_l <- as.numeric(drop(t(c_vec) %*% V_l %*% c_vec))
      if (var_u < 0) stop("c' V_u c must be non-negative.", call. = FALSE)
      if (var_l <= 0) stop("c' V_l c must be positive.", call. = FALSE)
      list(
        var_u = var_u,
        var_l = var_l,
        power_fn = function(n) power_ppi_ols(
          delta = delta,
          V_u = V_u,
          V_l = V_l,
          N = N,
          n = n,
          c = c_vec,
          alpha = alpha
        )
      )
    },
    custom = {
      if (is.null(var_unlabeled) || is.null(var_labeled)) {
        stop("var_unlabeled and var_labeled must be supplied when type = 'custom'.",
             call. = FALSE)
      }
      var_u <- as.numeric(var_unlabeled)
      var_l <- as.numeric(var_labeled)
      if (var_u < 0) stop("var_unlabeled must be non-negative.", call. = FALSE)
      if (var_l <= 0) stop("var_labeled must be positive.", call. = FALSE)
      list(
        var_u = var_u,
        var_l = var_l,
        power_fn = function(n) {
          v <- var_u / N + var_l / n
          se <- sqrt(v)
          z_alpha <- stats::qnorm(1 - alpha / 2)
          mu <- abs(delta) / se
          1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
        }
      )
    }
  )

  var_u <- components$var_u
  var_l <- components$var_l

  c0 <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  denom <- (delta^2 / c0^2) - (var_u / N)

  if (denom <= 0) {
    stop(
      paste0(
        "Infeasible: unlabeled component dominates. Increase N or delta, ",
        "or relax target power/alpha."
      ),
      call. = FALSE
    )
  }

  n_star <- ceiling(var_l / denom)

  if (n_star <= N) {
    return(as.integer(n_star))
  }

  if (mode == "error") {
    stop(
      sprintf(
        "Infeasible: required n = %d exceeds N = %d. Increase N or delta, or relax power/alpha.",
        n_star, as.integer(N)
      ),
      call. = FALSE
    )
  }

  n_capped <- as.integer(N)
  ach_power <- components$power_fn(n_capped)
  warning(
    sprintf(
      "Required n = %d exceeds N = %d. Capping to n = N.\nAchieved power = %.4f (target: %.4f).",
      n_star, n_capped, ach_power, power
    ),
    call. = FALSE
  )

  attr(n_capped, "achieved_power") <- ach_power
  n_capped
}

#' Required labeled sample size for desired PPI++ power
#'
#' @description
#' Finds the minimum labeled sample size \eqn{n} that achieves a target two-sided
#' power for the PP++ mean estimator when the oracle blend \eqn{\lambda^\star} is
#' used. You can supply raw moments (`sigma_y2`, `sigma_f2`, `cov_y_f`) or reuse
#' the metrics interface to infer them.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}; must be non-zero.
#' @param N Unlabeled sample size.
#' @param alpha Two-sided significance level.
#' @param power Target power (\eqn{1-\beta}).
#' @param sigma_y2 Optional outcome variance; overrides anything implied by `metrics`.
#' @param sigma_f2 Optional prediction variance; overrides anything implied by `metrics`.
#' @param cov_y_f Optional covariance \eqn{\Cov(Y, f(X))}. When supplied (directly
#'   or via `metrics$cov_y_f`) the PPI++ power is returned in addition to the PP
#'   quantities.
#' @param var_f Optional variance of \eqn{f(X)}.
#' @param var_res Optional residual variance of residual \eqn{Y - f(X)}.
#' @param metrics Optional list of predictive metrics (passed to 
#'   `resolve_ppi_pp_moments()`).
#' @param metric_type Character string describing `metrics` (e.g. `"continuous"`,
#'   `"hard"`, `"prob"`, `"precision_recall"`).
#' @param m_labeled Labeled sample size associated with the metrics; defaults to
#'   `metrics$m_obs` when present.
#' @param correction Logical; apply finite-sample variance corrections when
#'   deriving moments from metrics (default `TRUE`).
#' @param warn_smallN Logical; warn if `N < smallN_threshold`.
#' @param smallN_threshold Threshold that triggers the small-`N` warning.
#' @param mode `"error"` (default) throws when the target power cannot be met with
#'   `n \le N`; `"cap"` returns `n = N` with the achieved power attached as an attribute.
#'
#' @return Integer labeled sample size; if `mode = "cap"` and the cap is hit, the
#'   returned value carries an `"achieved_power"` attribute.
#' 
#' @examples
#' # PP++ mean estimator using direct moment inputs
#' n_required_ppi_pp(
#'   delta = 0.2,
#'   N = 4000,
#'   alpha = 0.05,
#'   power = 0.8,
#'   sigma_y2 = 1.0,
#'   sigma_f2 = 0.4,
#'   cov_y_f = 0.18,
#'   var_f = 0.4,
#'   var_res = 0.6
#' )
#'
#' # Same calculation, recovering moments from predictive metrics
#' metrics_pp <- list(
#'   type = "continuous",
#'   mse = 0.6,
#'   var_y = 1.0,
#'   cov_y_f = 0.18,
#'   m_obs = 300
#' )
#' n_required_ppi_pp(
#'   delta = 0.2,
#'   N = 4000,
#'   alpha = 0.05,
#'   power = 0.8,
#'   metrics = metrics_pp,
#'   metric_type = "continuous"
#' )
#' @export
n_required_ppi_pp <- function(delta,
                              N,
                              alpha = 0.05,
                              power = 0.80,
                              sigma_y2 = NULL,
                              sigma_f2 = NULL,
                              cov_y_f = NULL,
                              var_f = NULL,
                              var_res = NULL,
                              metrics = NULL,
                              metric_type = NULL,
                              m_labeled = NULL,
                              correction = TRUE,
                              warn_smallN = TRUE,
                              smallN_threshold = 500,
                              mode = c("error", "cap")) {
  mode <- match.arg(mode)

  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta == 0) {
    stop("delta must be a non-zero finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || N <= 0) {
    stop("N must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(power) || power <= 0 || power >= 1) {
    stop("power must lie in (0, 1).", call. = FALSE)
  }

  if (warn_smallN && N < smallN_threshold) {
    warning("Unlabeled N is quite small. The finite-N term may dominate the variance.",
            call. = FALSE)
  }

  moments <- resolve_ppi_pp_moments(
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f,
    var_f    = var_f,
    var_res  = var_res,
    metrics  = metrics,
    metric_type = metric_type,
    m_labeled  = m_labeled,
    correction = correction
  )

  sigma_y2 <- moments$sigma_y2
  sigma_f2 <- moments$sigma_f2
  cov_y_f  <- moments$cov_y_f

  if (!is.finite(sigma_f2) || sigma_f2 <= 0) {
    stop("sigma_f2 must be positive to apply the PP++ oracle variance formula.", call. = FALSE)
  }

  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- stats::qnorm(power)
  c0 <- z_alpha + z_beta
  if (!is.finite(c0) || c0 <= 0) {
    stop("Unable to compute z-scores for the requested alpha/power.", call. = FALSE)
  }

  S2 <- (delta / c0)^2
  term_common <- sigma_y2 - (cov_y_f^2 / sigma_f2)
  discrim <- (sigma_y2 - S2 * N)^2 + 4 * S2 * N * term_common

  if (!is.finite(discrim) || discrim < 0) {
    stop("Infeasible: desired accuracy/power cannot be met with supplied moments and N.", call. = FALSE)
  }

  numerator <- sigma_y2 - S2 * N + sqrt(discrim)
  denom <- 2 * S2
  if (denom <= 0) {
    stop("Infeasible: denominator non-positive. Check delta, alpha, and power inputs.", call. = FALSE)
  }

  n_star <- ceiling(numerator / denom)
  if (!is.finite(n_star) || n_star <= 0) {
    stop("Infeasible: required n is non-positive or undefined. Check inputs.", call. = FALSE)
  }

  if (n_star <= N) {
    return(as.integer(n_star))
  }

  if (mode == "error") {
    stop(
      sprintf(
        "Infeasible: required n = %d exceeds N = %d. Increase N or delta, or relax power/alpha.",
        n_star, as.integer(N)
      ),
      call. = FALSE
    )
  }

  n_capped <- as.integer(N)
  achieved <- power_ppi_pp_mean(
    delta = delta,
    N = N,
    n = n_capped,
    alpha = alpha,
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = moments$var_f,
    var_res = moments$var_res
  )
  warning(
    sprintf(
      "Required n = %d exceeds N = %d. Capping to n = N.\nAchieved power = %.4f (target: %.4f).",
      n_star, n_capped, achieved, power
    ),
    call. = FALSE
  )
  attr(n_capped, "achieved_power") <- achieved
  n_capped
}
