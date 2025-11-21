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

#' Required labeled sample size for PPI++ (mean estimator or OLS contrast)
#'
#' @description
#' Computes the minimum number of labeled observations \eqn{n} required to achieve
#' a desired power for Prediction-Powered Plus (PPI++) inference.
#'
#' The function supports:
#'
#' * **Mean estimation** (scalar parameter)
#' * **OLS linear contrast estimation** \eqn{a^\top \beta}
#'
#' and three modes for the augmentation weight \eqn{\lambda}:
#'
#' * `"vanilla"` — \eqn{\lambda = 1} (standard PPI)
#' * `"oracle"`  — \eqn{\lambda = \lambda^\star} (variance-minimizing; ratio \eqn{r = n/N} updated iteratively)
#' * `"user"`    — \eqn{\lambda = \lambda_{\text{user}}} (user-specified)
#'
#' The unlabeled sample size \eqn{N} contributes predictions \eqn{f(X)} that reduce
#' the variance of the estimator. The function uses closed-form inversions when
#' available, or numerical root-finding otherwise, to compute the minimal
#' \eqn{n} satisfying the desired power constraint.
#' 
#' @param delta Numeric scalar. Effect size being tested.  
#'   * **Mean case:** \eqn{\theta - \theta_0}.  
#'   * **OLS case:** \eqn{a^\top \beta - \theta_0}.  
#'   Must be nonzero.
#' @param N Positive scalar. Size of the unlabeled dataset used to compute
#'   \eqn{\bar f_U}. Must satisfy \code{N > 0}.
#' @param alpha Significance level for a two-sided Wald test. Default: \code{0.05}.
#' @param power Desired statistical power, e.g., \code{0.80}.
#' @param type Character string specifying the estimation problem:
#'   * `"mean"` — mean estimator \eqn{\hat{\theta}_\lambda}.  
#'   * `"ols"`  — linear contrast \eqn{a^\top \hat{\beta}_\lambda}.
#' @param lambda_mode Character string specifying λ:
#'   * `"vanilla"` — λ = 1.  
#'   * `"oracle"`  — λ = λ\*(n), updated using r = n/N.  
#'   * `"user"`    — λ fixed at user-supplied \code{lambda_user}.
#' @param lambda_user Numeric. The λ value to use when \code{lambda_mode = "user"}. Otherwise ignored.
#' 
#' @param sigma_y2 Variance of the labeled outcomes \eqn{Y}. Optional if using
#'   \code{metrics}.
#' @param sigma_f2 Variance of the predictions \eqn{f(X)}. Optional if using
#'   \code{metrics}.
#' @param cov_y_f Covariance between \eqn{Y} and \eqn{f(X)}.
#' @param var_f Alias for \code{sigma_f2}. Provided for flexibility.
#' @param var_res Residual variance \eqn{\Var(Y - f(X))}. Used if available.
#' @param metrics Optional list containing predictive performance summaries such as:
#'   * mean squared error (MSE)  
#'   * \eqn{R^2}
#'   * empirical covariance estimates  
#'   Used as an alternative interface when \code{sigma_*} are not provided.
#'
#' @param metric_type Character. One of `"continuous"` or `"binary"`.
#'   Determines how prediction metrics are translated to required moments.
#' @param m_labeled Number of labeled samples used to estimate prediction metrics.
#'   Needed for small-sample corrections.
#' @param correction Logical. If TRUE (default), applies finite-sample corrections
#'   to adjust \eqn{\sigma_f^2} and \eqn{\Cov(Y,f)}.
#'  
#' @param c Numeric vector defining the contrast \eqn{a}. Must have length p.
#'
#' @param H_L Per-observation Hessians for labeled design matrices.
#' @param H_U Per-observation Hessians for labeled and unlabeled design
#'   matrices. For OLS, \eqn{H = X^\top X / n}. These are used in the sandwich
#'   variance formula.
#'
#' @param Sigma_YY Sandwich covariance term \eqn{\Sigma_{YY}}.
#' @param Sigma_ff_l Sandwich covariance terms \eqn{\Sigma_{ff}}
#'   computed on labeled samples.
#' @param Sigma_ff_u Sandwich covariance terms \eqn{\Sigma_{ff}}
#'   computed on unlabeled samples.
#' @param Sigma_Yf Cross-covariance sandwich term \eqn{\Sigma_{Yf}}.
#'
#' @param warn_smallN Logical. If TRUE, warns when \code{N < smallN_threshold}.
#'
#' @param smallN_threshold Numeric threshold at which to trigger a warning
#'   (default: 500).
#'
#' @param mode How to behave if the required n exceeds N:
#'   * `"error"` — throw an error (default)
#'   * `"cap"`   — return n = N and attach attribute `"achieved_power"`.
#'
#' @return
#' Integer labeled sample size `n_star`.  
#' When `mode = "cap"` and `n_star > N`, returns `N` and attaches
#' `attr(n_star, "achieved_power")`.
#' 
#' @export
#' 
#' @importFrom stats uniroot
#' 
#' @examples
#' 
#' metrics <- list(
#'   type = "continuous",
#'   mse = 0.6,
#'   var_y = 1.0,
#'   cov_y_f = 0.18,
#'   m_obs = 300
#' )
#'
#' n_required_ppi_pp(
#'   delta = 0.25, N = 5000,
#'   type = "mean",
#'   metrics = metrics,
#'   metric_type = "continuous",
#'   lambda_mode = "vanilla"
#' )
#'
#' n_required_ppi_pp(
#'   delta = 0.25, N = 5000,
#'   type = "mean",
#'   sigma_y2 = 1.0,
#'   sigma_f2 = 0.4,
#'   cov_y_f  = 0.18,
#'   var_f    = 0.4,          # must supply
#'   var_res  = 1.0 - 0.4,    # = 0.6
#'   lambda_mode = "oracle"
#' )
#' 
#' n_required_ppi_pp(
#'   delta = 0.25, N = 5000,
#'   type = "mean",
#'   sigma_y2 = 1.0,
#'   sigma_f2 = 0.4,
#'   cov_y_f  = 0.18,
#'   var_f = 0.4,
#'   var_res = 0.6,
#'   lambda_mode = "user",
#'   lambda_user = 0.5
#' )
#' 
#' p <- 3
#' cvec <- c(1, 0, -1)
#' Hpop <- diag(p)
#' SYY  <- diag(p)
#' Sff  <- 0.5 * diag(p)
#' SYf  <- 0.3 * diag(p)
#'  
#' ## Vanilla λ=1
#' n_required_ppi_pp(
#'   delta = 0.15, N = 2000,
#'   type = "ols",
#'   c = cvec,
#'   H_L = Hpop, H_U = Hpop,
#'   Sigma_YY = SYY,
#'   Sigma_ff_l = Sff, Sigma_ff_u = Sff,
#'   Sigma_Yf   = SYf,
#'   lambda_mode = "vanilla"
#' )
#' 
#' ## Oracle λ*
#' n_required_ppi_pp(
#'   delta = 0.15, N = 2000,
#'   type = "ols",
#'   c = cvec,
#'   H_L = Hpop, H_U = Hpop,
#'   Sigma_YY = SYY,
#'   Sigma_ff_l = Sff, Sigma_ff_u = Sff,
#'   Sigma_Yf   = SYf,
#'   lambda_mode = "oracle"
#' )

n_required_ppi_pp <- function(
  delta, 
  N,
  alpha = 0.05, 
  power = 0.80,
  type = c("mean", "ols"),
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  # mean inputs:
  sigma_y2 = NULL, 
  sigma_f2 = NULL, 
  cov_y_f = NULL,
  var_f = NULL, 
  var_res = NULL, 
  metrics = NULL,
  metric_type = NULL, 
  m_labeled = NULL, 
  correction = TRUE,  # Correct for variance or not
  # OLS inputs:
  c = NULL, 
  H_L = NULL, 
  H_U = NULL,
  Sigma_YY = NULL, 
  Sigma_ff_l = NULL,
  Sigma_ff_u = NULL, 
  Sigma_Yf = NULL,
  warn_smallN = TRUE, 
  smallN_threshold = 500,
  mode = c("error", "cap")
) {
  type <- match.arg(type)
  lambda_mode <- match.arg(lambda_mode)
  mode <- match.arg(mode)

  ## Basic validation
  if (!is.numeric(delta) || length(delta) != 1 || delta == 0)
    stop("delta must be a non-zero numeric scalar.")
  if (!is.numeric(N) || N <= 0)
    stop("N must be a positive scalar.")

  if (warn_smallN && N < smallN_threshold)
    warning("N is small; finite-N effects may dominate.")

  ## z-values
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2

  ####  MEAN ESTIMATION CASE

  if (type == "mean") {

    ## Extract or compute predictive moments
    mm <- resolve_ppi_pp_moments(
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f  = cov_y_f,
      var_f    = var_f,
      var_res  = var_res,
      metrics  = metrics,
      metric_type = metric_type,
      m_labeled = m_labeled,
      correction = correction
    )
    sy2 <- mm$sigma_y2
    sf2 <- mm$sigma_f2
    cyf <- mm$cov_y_f

    if (!is.finite(sf2) || sf2 <= 0) {
  stop("sigma_f2 must be positive", call. = FALSE)
    }

    # Optional consistency check:
    if (cyf^2 > sy2 * sf2 + 1e-8) {
      stop("Infeasible: moments imply |Cov(Y,f)| > sqrt(Var(Y)Var(f)).", call. = FALSE)
    }

    ## Choose lambda
    r <- NA
    if (lambda_mode == "vanilla") {
      var_fun <- function(n) {
        sy2/n + sf2/N + sf2/n - 2*(cyf/n)
      }
    } else if (lambda_mode == "user") {
      if (is.null(lambda_user))
        stop("lambda_user must be provided when lambda_mode='user'")
      lambda <- lambda_user

      var_fun <- function(n) {
        lambda^2 * (sf2/N + sf2/n) + sy2/n - 2 * lambda * (cyf/n)
      }

    } else if (lambda_mode == "oracle") {

      ## exact closed-form oracle variance (Section 3)
      ## Var(theta_hat_{λ*})
      var_fun <- function(n) {
        sy2/n - (cyf^2 / sf2) * (N / (n*(n+N)))
      }
    }

    ## Solve var_fun(n) ≤ S2
    obj <- function(n) var_fun(n) - S2
    root <- tryCatch(
      uniroot(obj, c(2, 1e6))$root,
      error = function(e) NA_real_
    )

    if (is.na(root)) {
      if (mode == "error") {
        stop("Infeasible: cannot find n achieving required power.")
      } else {
        # mode == "cap": cap at n = N and report achieved power
        n_capped <- as.integer(N)
        varN <- var_fun(N)
        achieved <- 1 - pnorm(z_alpha - delta / sqrt(varN))
        warning(
          sprintf(
            "Required n exceeds search interval; capping to n = N. Achieved power = %.4f (target: %.4f).",
            achieved, power
          ),
          call. = FALSE
        )
        attr(n_capped, "achieved_power") <- achieved
        return(n_capped)
      }
    }

    ## Root found, return n*
    n_star <- as.integer(ceiling(root))

    ## If feasible
    if (n_star <= N) return(n_star)

    ## n_star > N
    if (mode == "error") {
      stop(sprintf("Required n=%d exceeds N=%d.", n_star, N), call. = FALSE)
    }
    ## mode = "cap"
    n_capped <- as.integer(N)
    varN <- var_fun(N)
    achieved <- 1 - pnorm(z_alpha - delta / sqrt(varN))
    attr(n_capped, "achieved_power") <- achieved
    return(n_capped)

  }

  ####  OLS CONTRAST CASE

  if (type == "ols") {
    ## Check inputs
    if (is.null(c) || is.null(H_L) || is.null(H_U) ||
        is.null(Sigma_YY) || is.null(Sigma_ff_l) ||
        is.null(Sigma_ff_u) || is.null(Sigma_Yf))
      stop("OLS requires c, H_L, H_U, Sigma_YY, Sigma_ff_l, Sigma_ff_u, Sigma_Yf.")

    H_L <- as.matrix(H_L)
    H_U <- as.matrix(H_U)
    c   <- as.numeric(c)

    # Objects for population level Hessians & contrasts

    H_pop <- H_U
    Hinv_pop <- solve(H_pop)

    A_pop <- as.numeric(t(c) %*% Hinv_pop %*% Sigma_ff_u %*% Hinv_pop %*% c)
    C_pop <- as.numeric(t(c) %*% Hinv_pop %*% Sigma_Yf   %*% Hinv_pop %*% c)

    ## Choose lambda
    if (lambda_mode == "vanilla") {
      lambda <- 1
    } else if (lambda_mode == "user") {
      if (is.null(lambda_user))
        stop("lambda_user must be provided for lambda_mode='user'")
      lambda <- lambda_user
    }

    if (lambda_mode %in% c("vanilla", "user")) {
      lambda <- if (lambda_mode == "vanilla") 1 else lambda_user

      var_fun <- function(n) {
        H_mix <- (1 - lambda)*H_L + lambda*H_U
        Hinv  <- solve(H_mix)
        middle <- Sigma_YY/n +
          lambda^2 * (Sigma_ff_u/N + Sigma_ff_l/n) -
          2 * lambda * (Sigma_Yf/n)
        as.numeric(t(c) %*% Hinv %*% middle %*% Hinv %*% c)
      }

    } else {

      var_fun <- function(n) {
        r <- n/N
        lambda_star <- C_pop / ((1+r)*A_pop)
        lambda_star <- max(0, min(1, lambda_star))

        # Hessian for variance at λ*
        H_mix <- (1 - lambda_star)*H_L + lambda_star*H_U
        Hinv  <- solve(H_mix)

        middle <- Sigma_YY/n +
          lambda_star^2 * (Sigma_ff_u/N + Sigma_ff_l/n) -
          2 * lambda_star * (Sigma_Yf/n)

        out <- as.numeric(t(c) %*% Hinv %*% middle %*% Hinv %*% c)
        out
      }
    }


    obj <- function(n) var_fun(n) - S2
    root <- tryCatch(uniroot(obj, interval = c(2, 1e6))$root,
                    error = function(e) NA_real_)

    # if infeasible

    if (is.na(root)) {
      if (mode == "error") {
        stop("Infeasible: cannot find n achieving required power.")
      } else {
        n_capped <- as.integer(N)
        vN <- var_fun(N)
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
      stop(sprintf("Required n=%d exceeds N=%d.", n_star, N))
    }

    # mode = "cap" 
    n_capped <- as.integer(N)
    vN <- var_fun(N)
    achieved <- 1 - pnorm(z_alpha - delta / sqrt(vN))
    attr(n_capped, "achieved_power") <- achieved
    return(n_capped)
  }
}