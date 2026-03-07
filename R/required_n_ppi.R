# ============================================================================
# Internal: sample-size solver for PPI++ mean estimation
# Called by the unified power_ppi_mean() when n = NULL
# ============================================================================

#' @importFrom stats glm.fit uniroot
n_required_ppi_pp_mean_internal <- function(
  delta,
  N,
  alpha = 0.05,
  power = 0.80,
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f  = NULL,
  var_f    = NULL,
  var_res  = NULL,
  metrics  = NULL,
  metric_type = NULL,
  m_labeled = NULL,
  correction = TRUE,
  warn_smallN = TRUE,
  smallN_threshold = 500,
  mode = c("error", "cap")
) {
  lambda_mode <- match.arg(lambda_mode)
  mode        <- match.arg(mode)

  ## Basic validation
  if (!is.numeric(delta) || length(delta) != 1L || delta == 0)
    stop("delta must be a non-zero numeric scalar.")
  if (!is.numeric(N) || length(N) != 1L || N <= 0)
    stop("N must be a positive scalar.")

  if (warn_smallN && N < smallN_threshold)
    warning("N is small; finite-N effects may dominate.")

  ## z-values and target variance S^2
  z_alpha <- stats::qnorm(1 - alpha/2)
  z_beta  <- stats::qnorm(power)
  S2      <- (delta / (z_alpha + z_beta))^2

  ## Binary metric helpers
  infer_p_y_if_missing <- function(metrics) {
    if (!is.null(metrics$p_y)) return(metrics$p_y)
    if (!is.null(metrics$tp) && !is.null(metrics$fn) && !is.null(metrics$m_obs)) {
      return((metrics$tp + metrics$fn) / metrics$m_obs)
    }
    NULL
  }

  infer_var_y <- function(metrics) {
    if (!is.null(metrics$var_y)) return(metrics$var_y)
    if (!is.null(metrics$p_y)) return(metrics$p_y * (1 - metrics$p_y))
    NULL
  }

  ## A. Direct-moments mode: user supplies sigma_y2, sigma_f2, cov_y_f
  if (!is.null(sigma_y2) || !is.null(sigma_f2) || !is.null(cov_y_f)) {

    if (is.null(sigma_y2) || is.null(sigma_f2) || is.null(cov_y_f)) {
      stop("If supplying direct moments, provide sigma_y2, sigma_f2, and cov_y_f.",
           call. = FALSE)
    }

    sy2 <- as.numeric(sigma_y2)
    sf2 <- as.numeric(sigma_f2)
    cyf <- as.numeric(cov_y_f)

    if (!is.finite(sf2) || sf2 <= 0) {
      stop("sigma_f2 must be positive.", call. = FALSE)
    }
    if (!is.finite(sy2) || sy2 < 0) {
      stop("sigma_y2 must be non-negative.", call. = FALSE)
    }

    if (cyf^2 > sy2 * sf2 + 1e-8) {
      stop("Infeasible: |Cov(Y,f)| > sqrt(Var(Y)Var(f)).", call. = FALSE)
    }

  } else {

    ## Metrics-based mode: derive moments from metrics/var_f/var_res

    if (is.null(metrics))
      stop("Must supply metrics list when moments are not provided.")

    # infer m_labeled
    if (is.null(m_labeled)) {
      m_labeled <- metrics$m_obs %||% stop("metrics$m_obs is required.")
    }

    # infer metric type
    metric_type <- metric_type %||% metrics$type %||%
      stop("metric_type must be provided or stored in metrics$type.")

    metric_type_clean <- match.arg(tolower(metric_type),
                                   c("continuous", "classification", "prob"))

    # auto-infer p_y if classification
    if (metric_type_clean != "continuous") {
      metrics$p_y <- metrics$p_y %||% infer_p_y_if_missing(metrics)
    }

    vars <- resolve_ppi_variances(
      var_f       = var_f,
      var_res     = var_res,
      metrics     = metrics,
      metric_type = metric_type,
      m_labeled   = m_labeled,
      correction  = correction
    )

    sf2  <- vars$var_f       # Var(f)
    vres <- vars$var_res     # Var(Y - f)

    if (sf2 <= 0) stop("Var(f) must be positive.")

    sy2 <- infer_var_y(metrics)
    if (is.null(sy2))
      stop("metrics must include var_y (continuous) or p_y (binary).")

    cyf <- (sy2 + sf2 - vres) / 2
    if (cyf^2 > sy2 * sf2 + 1e-8) stop("Infeasible: |Cov(Y,f)| > sqrt(Var(Y)Var(f)).")
  }

  ## C. Define variance(n) under lambda_mode
  if (lambda_mode == "vanilla") {

    var_fun <- function(n) {
      sy2 / n + sf2 / N + sf2 / n - 2 * (cyf / n)
    }

  } else if (lambda_mode == "user") {

    if (is.null(lambda_user))
      stop("lambda_user must be provided when lambda_mode = 'user'.", call. = FALSE)
    lambda <- lambda_user

    var_fun <- function(n) {
      lambda^2 * (sf2/N + sf2/n) +
        sy2/n -
        2 * lambda * (cyf/n)
    }

  } else if (lambda_mode == "oracle") {

    ## Exact closed-form oracle variance (Section 3)
    var_fun <- function(n) {
      sy2/n - (cyf^2 / sf2) * (N / (n * (n + N)))
    }
  }

  ## Solve var_fun(n) <= S2
  obj  <- function(n) var_fun(n) - S2
  root <- tryCatch(
    uniroot(obj, interval = c(2, 1e6))$root,
    error = function(e) NA_real_
  )

  if (is.na(root)) {
    if (mode == "error") {
      stop("Infeasible: cannot find n satisfying required power.", call. = FALSE)
    } else {
      n_capped <- as.integer(N)
      varN     <- var_fun(N)
      achieved <- 1 - pnorm(z_alpha - delta / sqrt(varN))
      warning(sprintf(
        "Required n exceeds search range; capping to n = N = %d. Achieved power = %.4f.",
        N, achieved
      ), call. = FALSE)
      attr(n_capped, "achieved_power") <- achieved
      return(n_capped)
    }
  }

  ## Valid root
  n_star <- as.integer(ceiling(root))
  if (n_star <= N) return(n_star)

  ## Root exists but n > N
  if (mode == "error") {
    stop(sprintf("Required n = %d exceeds N = %d.", n_star, N), call. = FALSE)
  }

  n_capped <- as.integer(N)
  varN     <- var_fun(N)
  achieved <- 1 - pnorm(z_alpha - delta / sqrt(varN))
  attr(n_capped, "achieved_power") <- achieved
  return(n_capped)
}

# ============================================================================
# power_ppi_regression(): extracted from the old n_required_ppi_pp(type="regression")
# ============================================================================

#' Power or required sample size for PPI++ regression contrast
#'
#' @description
#' Computes the minimum labeled sample size \eqn{n} for a PPI++ regression
#' contrast, or computes power given `n`. Follows the \pkg{pwr} convention:
#' exactly one of `n` or `power` must be `NULL`.
#'
#' Supports both OLS and canonical GLM contrasts via sandwich/Fisher blocks
#' computed on labeled and unlabeled samples.
#'
#' @param delta Numeric scalar. Effect size \eqn{a^\top \beta - \theta_0}. Must be nonzero when solving for `n`.
#' @param N Positive scalar. Number of unlabeled observations.
#' @param n Labeled sample size. Set to `NULL` to solve for required `n`.
#' @param power Desired power. Set to `NULL` to compute power from `n`.
#' @param alpha Significance level (default `0.05`).
#' @param lambda_mode `"vanilla"`, `"oracle"`, or `"user"`.
#' @param lambda_user Scalar specifying \eqn{\lambda} when `lambda_mode = "user"`.
#' @param c Contrast vector \eqn{a}.
#' @param H_L Labeled Hessian/Fisher block.
#' @param H_U Unlabeled Hessian/Fisher block.
#' @param Sigma_YY Labeled covariance block.
#' @param Sigma_ff_l Labeled prediction covariance block.
#' @param Sigma_ff_u Unlabeled prediction covariance block.
#' @param Sigma_Yf Cross-covariance block.
#' @param warn_smallN Logical. Warn when `N < smallN_threshold`.
#' @param smallN_threshold Numeric. Default `500`.
#' @param mode `"error"` or `"cap"`.
#'
#' @return When computing power: scalar in \[0, 1\].
#'   When computing sample size: integer `n` required.
#' @export
#'
#' @examples
#' p <- 3
#' cvec <- c(1, 0, -1)
#' Hpop <- diag(p)
#' SYY  <- diag(p)
#' Sff  <- 0.5 * diag(p)
#' SYf  <- 0.3 * diag(p)
#'
#' # Compute required sample size
#' power_ppi_regression(
#'   delta = 0.15, N = 2000, power = 0.90,
#'   c = cvec, H_L = Hpop, H_U = Hpop,
#'   Sigma_YY = SYY, Sigma_ff_l = Sff, Sigma_ff_u = Sff,
#'   Sigma_Yf = SYf, lambda_mode = "oracle"
#' )
#'
power_ppi_regression <- function(
  delta,
  N,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  c = NULL,
  H_L = NULL,
  H_U = NULL,
  Sigma_YY   = NULL,
  Sigma_ff_l = NULL,
  Sigma_ff_u = NULL,
  Sigma_Yf   = NULL,
  warn_smallN = TRUE,
  smallN_threshold = 500,
  mode = c("error", "cap")
) {
  lambda_mode <- match.arg(lambda_mode)
  mode        <- match.arg(mode)

  # Validate: exactly one of n, power must be NULL
  if (is.null(n) && is.null(power)) {
    stop("Exactly one of 'n' and 'power' must be NULL.", call. = FALSE)
  }
  if (!is.null(n) && !is.null(power)) {
    stop("Exactly one of 'n' and 'power' must be NULL.", call. = FALSE)
  }

  if (!is.numeric(delta) || length(delta) != 1L || delta == 0)
    stop("delta must be a non-zero numeric scalar.")
  if (!is.numeric(N) || length(N) != 1L || N <= 0)
    stop("N must be a positive scalar.")

  if (warn_smallN && N < smallN_threshold)
    warning("N is small; finite-N effects may dominate.")

  if (is.null(c) || is.null(H_L) || is.null(H_U) ||
      is.null(Sigma_YY) || is.null(Sigma_ff_l) ||
      is.null(Sigma_ff_u) || is.null(Sigma_Yf)) {
    stop(
      "Regression contrast case requires c, H_L, H_U, Sigma_YY, Sigma_ff_l, Sigma_ff_u, Sigma_Yf.\n",
      "You can obtain these via compute_ppi_blocks(model_type, X_l, Y_l, f_l, X_u, f_u, beta, family).",
      call. = FALSE
    )
  }

  H_L <- as.matrix(H_L)
  H_U <- as.matrix(H_U)
  c   <- as.numeric(c)

  ## Population-level Hessian / Fisher and contrast constants
  H_pop    <- H_U
  Hinv_pop <- solve(H_pop)

  A_pop <- as.numeric(t(c) %*% Hinv_pop %*% Sigma_ff_u %*% Hinv_pop %*% c)
  C_pop <- as.numeric(t(c) %*% Hinv_pop %*% Sigma_Yf   %*% Hinv_pop %*% c)

  ## Choose lambda
  if (lambda_mode == "vanilla") {
    lambda <- 1
  } else if (lambda_mode == "user") {
    if (is.null(lambda_user))
      stop("lambda_user must be provided for lambda_mode='user'.", call. = FALSE)
    lambda <- lambda_user
  }

  if (lambda_mode %in% c("vanilla", "user")) {
    lambda <- if (lambda_mode == "vanilla") 1 else lambda_user

    var_fun <- function(n) {
      H_mix <- (1 - lambda) * H_L + lambda * H_U
      Hinv  <- solve(H_mix)
      middle <- Sigma_YY/n +
        lambda^2 * (Sigma_ff_u/N + Sigma_ff_l/n) -
        2 * lambda * (Sigma_Yf/n)
      as.numeric(t(c) %*% Hinv %*% middle %*% Hinv %*% c)
    }

  } else {
    ## Oracle lambda*(a) and associated variance
    var_fun <- function(n) {
      r <- n / N
      lambda_star <- C_pop / ((1 + r) * A_pop)
      lambda_star <- max(0, min(1, lambda_star))

      H_mix <- (1 - lambda_star) * H_L + lambda_star * H_U
      Hinv  <- solve(H_mix)

      middle <- Sigma_YY/n +
        lambda_star^2 * (Sigma_ff_u/N + Sigma_ff_l/n) -
        2 * lambda_star * (Sigma_Yf/n)

      as.numeric(t(c) %*% Hinv %*% middle %*% Hinv %*% c)
    }
  }

  if (!is.null(n)) {
    # ---- Compute power given n ----
    v <- var_fun(n)
    se <- sqrt(v)
    z_alpha <- stats::qnorm(1 - alpha / 2)
    mu <- abs(delta) / se
    return(1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu))
  }

  # ---- Solve for n given target power ----
  z_alpha <- stats::qnorm(1 - alpha/2)
  z_beta  <- stats::qnorm(power)
  S2      <- (delta / (z_alpha + z_beta))^2

  obj  <- function(n) var_fun(n) - S2
  root <- tryCatch(uniroot(obj, interval = c(2, 1e6))$root,
                   error = function(e) NA_real_)

  if (is.na(root)) {
    if (mode == "error") {
      stop("Infeasible: cannot find n achieving required power.", call. = FALSE)
    } else {
      n_capped <- as.integer(N)
      vN       <- var_fun(N)
      achieved <- 1 - pnorm(z_alpha - delta / sqrt(vN))
      warning(sprintf(
        "Required n exceeds search interval; capping to n = N. Achieved power = %.4f (target %.4f).",
        achieved, power
      ), call. = FALSE)
      attr(n_capped, "achieved_power") <- achieved
      return(n_capped)
    }
  }

  n_star <- as.integer(ceiling(root))
  if (n_star <= N) return(n_star)

  if (mode == "error") {
    stop(sprintf("Required n=%d exceeds N=%d.", n_star, N), call. = FALSE)
  }

  n_capped <- as.integer(N)
  vN       <- var_fun(N)
  achieved <- 1 - pnorm(z_alpha - delta / sqrt(vN))
  attr(n_capped, "achieved_power") <- achieved
  return(n_capped)
}

# ============================================================================
# Deprecated shim: n_required_ppi_pp()
# Routes to power_ppi_mean(..., n=NULL) or power_ppi_regression(..., n=NULL)
# ============================================================================

#' @rdname deprecated
#' @export
n_required_ppi_pp <- function(
  delta,
  N,
  alpha = 0.05,
  power = 0.80,
  type = c("mean", "regression"),
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f  = NULL,
  var_f    = NULL,
  var_res  = NULL,
  metrics  = NULL,
  metric_type = NULL,
  m_labeled = NULL,
  correction = TRUE,
  c = NULL,
  H_L = NULL,
  H_U = NULL,
  Sigma_YY   = NULL,
  Sigma_ff_l = NULL,
  Sigma_ff_u = NULL,
  Sigma_Yf   = NULL,
  warn_smallN = TRUE,
  smallN_threshold = 500,
  mode = c("error", "cap")
) {
  type        <- match.arg(type)
  lambda_mode <- match.arg(lambda_mode)
  mode        <- match.arg(mode)

  .Deprecated(
    new = if (type == "mean") "power_ppi_mean" else "power_ppi_regression",
    package = "pppower",
    msg = paste0(
      "n_required_ppi_pp() is deprecated. Use ",
      if (type == "mean") "power_ppi_mean(..., n = NULL)" else "power_ppi_regression(..., n = NULL)",
      " instead."
    )
  )

  if (type == "mean") {
    power_ppi_mean(
      delta = delta,
      N = N,
      n = NULL,
      power = power,
      alpha = alpha,
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f,
      var_f = var_f,
      var_res = var_res,
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = m_labeled,
      correction = correction,
      lambda_mode = lambda_mode,
      lambda_user = lambda_user,
      mode = mode,
      warn_smallN = warn_smallN,
      smallN_threshold = smallN_threshold
    )
  } else {
    power_ppi_regression(
      delta = delta,
      N = N,
      n = NULL,
      power = power,
      alpha = alpha,
      lambda_mode = lambda_mode,
      lambda_user = lambda_user,
      c = c,
      H_L = H_L,
      H_U = H_U,
      Sigma_YY = Sigma_YY,
      Sigma_ff_l = Sigma_ff_l,
      Sigma_ff_u = Sigma_ff_u,
      Sigma_Yf = Sigma_Yf,
      warn_smallN = warn_smallN,
      smallN_threshold = smallN_threshold,
      mode = mode
    )
  }
}
