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

#' Required labeled sample size for PPI++ (mean or OLS)
#'
#' @description
#' Computes the minimum labeled sample size (`n_star`) required to achieve a desired
#' statistical power under the Prediction-Powered Plus (PPI++) framework.
#' Supports both mean estimation and linear regression (PPI++–OLS).
#'
#' @param delta Effect size. For mean estimation, `theta - theta0`; for OLS, `t(c) %*% beta - theta0`.
#' @param N Unlabeled sample size.
#' @param alpha Significance level.
#' @param power Desired power (e.g., 0.8 for 80%).
#' @param type One of `"mean"` (default) or `"ols"`.
#'
#' @param sigma_y2 Variance of Y.
#' @param sigma_f2 Variance of predictions f(X).
#' @param cov_y_f Covariance between Y and f(X).
#' @param var_f Variance of f(X); used interchangeably with `sigma_f2`.
#' @param var_res Residual variance of Y - f(X).
#' @param metrics Optional list of predictive metrics (e.g., MSE, R2) from model validation.
#' @param metric_type Type of predictive metric: `"continuous"` or `"binary"`.
#' @param m_labeled Number of labeled samples used to estimate metrics.
#' @param correction Logical; apply small-sample correction if `TRUE`.
#'
#' @param c Contrast vector for OLS contrasts (`t(c) %*% beta`).
#' @param H_L,H_U Per-observation Hessians for labeled and unlabeled data (used in sandwich variance for OLS).
#' @param Sigma_YY Covariance of X(Y - X'beta*).
#' @param Sigma_ff_l,Sigma_ff_u Covariances of X(f(X) - X'beta*) for labeled/unlabeled data.
#' @param Sigma_Yf Cross-covariance between labeled and predicted score vectors.
#'
#' @param warn_smallN Warn if unlabeled N is small enough that finite-sample effects dominate.
#' @param smallN_threshold Threshold for triggering small-N warning (default = 500).
#' @param mode Behavior if required n_star exceeds N.
#'   Either `"error"` (default) to throw an error, or `"cap"` to cap n_star = N and return achieved power.
#'
#' @return Integer labeled sample size (`n_star`).
#'
#' @examples
#' # PPI++ mean estimator using direct moments
#' n_required_ppi_pp(
#'   delta = 0.2, N = 4000, alpha = 0.05, power = 0.8,
#'   sigma_y2 = 1.0, sigma_f2 = 0.4, cov_y_f = 0.18,
#'   var_f = 0.4, var_res = 0.6
#' )
#'
#' # Using predictive metrics
#' metrics_pp <- list(type = "continuous", mse = 0.6, var_y = 1.0, cov_y_f = 0.18, m_obs = 300)
#' n_required_ppi_pp(
#'   delta = 0.2, N = 4000, alpha = 0.05, power = 0.8,
#'   metrics = metrics_pp, metric_type = "continuous"
#' )
#'
#' # OLS example
#' c_vec <- c(0, 1, 0)
#' H_pop <- diag(3)
#' Sigma_YY_pop <- diag(3)
#' Sigma_ff_pop <- 0.5 * diag(3)
#' Sigma_Yf_pop <- 0.3 * diag(3)
#' n_required_ppi_pp(
#'   delta = 0.15, N = 2000, alpha = 0.05, power = 0.8,
#'   type = "ols",
#'   c = c_vec, H_L = H_pop, H_U = H_pop,
#'   Sigma_YY = Sigma_YY_pop, Sigma_ff_l = Sigma_ff_pop,
#'   Sigma_ff_u = Sigma_ff_pop, Sigma_Yf = Sigma_Yf_pop
#' )
#'
#' @export
n_required_ppi_pp <- function(delta,
                              N,
                              alpha = 0.05,
                              power = 0.80,
                              type = c("mean", "ols"),
                              # --- mean inputs ---
                              sigma_y2 = NULL,
                              sigma_f2 = NULL,
                              cov_y_f = NULL,
                              var_f = NULL,
                              var_res = NULL,
                              metrics = NULL,
                              metric_type = NULL,
                              m_labeled = NULL,
                              correction = TRUE,
                              # --- OLS inputs ---
                              c = NULL,
                              H_L = NULL,
                              H_U = NULL,
                              Sigma_YY = NULL,
                              Sigma_ff_l = NULL,
                              Sigma_ff_u = NULL,
                              Sigma_Yf = NULL,
                              warn_smallN = TRUE,
                              smallN_threshold = 500,
                              mode = c("error", "cap")) {
  type <- match.arg(type)
  mode <- match.arg(mode)

  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta == 0)
    stop("delta must be a non-zero finite numeric scalar.", call. = FALSE)
  if (!is.numeric(N) || length(N) != 1L || N <= 0)
    stop("N must be a positive numeric scalar.", call. = FALSE)

  if (warn_smallN && N < smallN_threshold)
    warning("Unlabeled N is quite small; finite-N term may dominate.", call. = FALSE)

  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- stats::qnorm(power)
  c0 <- z_alpha + z_beta
  if (!is.finite(c0) || c0 <= 0)
    stop("Unable to compute z-scores for the requested alpha/power.", call. = FALSE)

  if (type == "mean") {
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
      stop("sigma_f2 must be positive", call. = FALSE)
    }

    S2 <- (delta / c0)^2
    term_common <- sigma_y2 - (cov_y_f^2 / sigma_f2)
    discrim <- (sigma_y2 - S2 * N)^2 + 4 * S2 * N * term_common
    if (discrim < 0)
      stop("Infeasible: power target cannot be met with supplied moments.")

    numerator <- sigma_y2 - S2 * N + sqrt(discrim)
    denom <- 2 * S2
    n_star <- ceiling(numerator / denom)

    achieved_fun <- power_ppi_pp_mean
    achieved_args <- list(
      delta = delta, N = N, n = N, alpha = alpha,
      sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f,
      var_f = moments$var_f, var_res = moments$var_res
    )
  }

  if (type == "ols") {
    if (is.null(c) || is.null(H_L) || is.null(H_U) ||
        is.null(Sigma_YY) || is.null(Sigma_ff_l) ||
        is.null(Sigma_ff_u) || is.null(Sigma_Yf)) {
      stop("OLS mode requires c, H_L, H_U, Sigma_YY, Sigma_ff_l, Sigma_ff_u, and Sigma_Yf.")
    }

    lambda <- 1
    H_mix <- (1 - lambda) * H_L + lambda * H_U
    bread_inv <- solve(H_mix)

    coef_var <- function(n) {
      middle <- Sigma_YY / n +
        lambda^2 * (Sigma_ff_u / N + Sigma_ff_l / n) -
        2 * lambda * Sigma_Yf / n
      as.numeric(t(c) %*% bread_inv %*% middle %*% bread_inv %*% c)
    }

    target <- (delta / c0)^2
    objective <- function(n) coef_var(n) - target

    root <- tryCatch(stats::uniroot(objective, c(2, 1e6))$root,
                 error = function(e) NA_real_)
    if (is.na(root))
      stop("Infeasible: unable to find n satisfying desired power.", call. = FALSE)

    n_star <- ceiling(root)

    achieved_fun <- power_ppi_pp_ols
    achieved_args <- list(
      delta = delta, N = N, n = N, alpha = alpha,
      contrast = c,
      H_L = H_L, H_U = H_U,
      Sigma_YY = Sigma_YY,
      Sigma_ff_l = Sigma_ff_l,
      Sigma_ff_u = Sigma_ff_u,
      Sigma_Yf = Sigma_Yf,
      lambda = lambda
    )
  }

  if (!is.finite(n_star) || n_star <= 0)
    stop("Infeasible: required n is non-positive or undefined.", call. = FALSE)

  if (n_star <= N)
    return(as.integer(n_star))

  if (mode == "error") {
    stop(
      sprintf("Infeasible: required n = %d exceeds N = %d.", n_star, as.integer(N)),
      call. = FALSE
    )
  }

  # CAP + WARNING
  n_capped <- as.integer(N)
  achieved <- do.call(achieved_fun, achieved_args)
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

n_required_ppi_pp_ols_closed <- function(
  delta, N, alpha = 0.05, power = 0.90,
  c, H, Sigma_YY, Sigma_ff, Sigma_Yf
) {
  stopifnot(is.numeric(delta), length(delta) == 1, is.finite(delta), delta != 0)
  stopifnot(is.numeric(N), length(N) == 1, N > 0)
  c  <- as.numeric(c)
  H  <- as.matrix(H)
  Sy <- as.matrix(Sigma_YY)
  Sf <- as.matrix(Sigma_ff)
  Syf<- as.matrix(Sigma_Yf)
  p <- length(c)
  if (!all(dim(H) == p)) stop("dim(H) must match length(c)")
  if (!all(dim(Sy) == p) || !all(dim(Sf) == p) || !all(dim(Syf) == p))
    stop("All Sigma_* must be p x p to match c")

  # z-constants and S^2
  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- stats::qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2
  if (!is.finite(S2) || S2 <= 0) stop("Invalid alpha/power/delta.")

  Hinv <- solve(H)
  A <- as.numeric(t(c) %*% Hinv %*% Sf  %*% Hinv %*% c)   # >= 0
  B <- as.numeric(t(c) %*% Hinv %*% Syf %*% Hinv %*% c)
  C <- as.numeric(t(c) %*% Hinv %*% Sy  %*% Hinv %*% c)   # > 0

  if (A <= 0 || C <= 0) stop("Require A>0 and C>0 for a meaningful contrast.")
  K <- B^2 / A                                           # <= C by Cauchy–Schwarz

  # Discriminant and closed-form root
  a <- S2 * N
  b <- S2 * N - C
  d <- (S2 * N + C)^2 - 4 * S2 * N * K
  if (!is.finite(d) || d < 0) {
    stop("Infeasible: power target cannot be met with supplied moments and N.")
  }

  r_star <- (C - S2 * N + sqrt(d)) / (2 * S2 * N)
  if (!is.finite(r_star) || r_star <= 0) {
    stop("Infeasible: required r <= 0. Increase N or |delta|, or relax alpha/power.")
  }

  n_star <- ceiling(r_star * N)
  as.integer(n_star)
}