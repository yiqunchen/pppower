#' Required labeled sample size for PPI++ (mean or regression/contrast inference)
#'
#' @description
#' Computes the minimum labeled sample size \eqn{n} required to achieve a desired
#' statistical power for Prediction-Powered Plus (PPI++) inference.
#'
#' This function supports:
#'
#' **1. Mean estimation**
#'
#' Estimating a scalar parameter \eqn{\theta = \mathbb{E}[Y]} using
#' prediction-powered estimators with augmentation weight \eqn{\lambda}.
#'
#' **2. Regression contrast estimation**
#'
#' Estimating a linear contrast \eqn{a^\top \beta^\star} in either:
#'   * **OLS** (ordinary least squares)
#'   * **GLM** with canonical link (e.g., logistic or Poisson regression)
#'
#' via sandwich/Fisher blocks computed on labeled and unlabeled samples.
#'
#' The augmentation weight \eqn{\lambda} may be:
#'
#' * `"vanilla"` — \eqn{\lambda = 1} (standard PPI)
#' * `"user"` — user-specified \eqn{\lambda}
#' * `"oracle"` — variance-minimizing \eqn{\lambda^\star(n)} that depends on \eqn{r = n/N}
#'
#' The unlabeled sample size \eqn{N} contributes predictions \eqn{f(X)} that lower
#' the estimator variance.
#' 
#' @param delta Numeric scalar. Effect size being tested.  
#'   * **Mean case:** difference from the null, \eqn{\theta - \theta_0}.  
#'   * **Regression case:** \eqn{a^\top \beta - \theta_0}.  
#'   Must be nonzero.
#' @param N Positive scalar. Number of unlabeled observations providing predictions.
#' @param alpha Significance level (two-sided). Default: `0.05`.
#' @param power Desired statistical power. Default: `0.80`.
#' @param type `"mean"` or `"regression"`.  
#'   `"regression"` covers both **OLS** and **canonical GLM** contrasts.
#' @param lambda_mode `"vanilla"`, `"oracle"`, or `"user"`.  
#' @param lambda_user Scalar specifying \eqn{\lambda} when
#'   `lambda_mode = "user"`.
#' @param sigma_y2 Variance of labeled outcomes \eqn{\Var(Y)}.
#' @param sigma_f2 Variance of predictions \eqn{\Var(f(X))}.
#' @param cov_y_f Covariance \eqn{\Cov(Y, f(X))}.
#' @param var_f Optional alias for `sigma_f2`.
#' @param var_res Residual variance \eqn{\Var(Y - f(X))}.
#' @param metrics List of predictive metrics such as:
#'   * `mse`  
#'   * `r2`  
#'   * `var_y`  
#'   * `p_y` (binary)  
#'   * `bias`  
#'   Required when direct moments are not supplied.
#'
#' @param metric_type `"continuous"` or `"binary"`.  
#'   Determines how `metrics` is interpreted.
#' @param m_labeled Number of labeled samples used to compute prediction metrics.
#' @param correction Logical; apply finite-sample adjustment when deriving metrics.
#' @param c Contrast vector \eqn{a}.
#' @param H_L Labeled Hessian/Fisher block.
#' @param H_U Unlabeled Hessian/Fisher block.
#' @param Sigma_YY Labeled covariance block \eqn{\Sigma_{YY}}.
#' @param Sigma_ff_l Labeled covariance block \eqn{\Sigma_{ff}}.
#' @param Sigma_ff_u Unlabeled covariance block \eqn{\Sigma_{ff}}.
#' @param Sigma_Yf Cross-covariance \eqn{\Sigma_{Yf}}.
#' @param warn_smallN Logical. Warn when `N < smallN_threshold`.  
#' @param smallN_threshold Numeric. Default `500`.
#' @param mode `"error"` or `"cap"`.  
#'   `"cap"` returns `n = N` with an `"achieved_power"` attribute.
#'
#'
#' @return Integer required sample size `n_star`.  
#' If capped, returns `N` with attribute `"achieved_power"`.
#'
#' @export
#' @importFrom stats glm.fit uniroot
#'
#' @examples
#' 
#' ## Example 1: Using prediction metrics
#' metrics <- list(
#'   type   = "continuous",
#'   mse    = 0.6,
#'   var_y  = 1.0,
#'   m_obs  = 300
#' )
#'
#' n_required_ppi_pp(
#'   delta       = 0.25,
#'   N           = 5000,
#'   type        = "mean",
#'   metrics     = metrics,
#'   metric_type = "continuous",
#'   lambda_mode = "vanilla"
#' )
#'
#'
#' ## Example 2: Direct-moment supply
#' n_required_ppi_pp(
#'   delta      = 0.25,
#'   N          = 5000,
#'   type       = "mean",
#'   sigma_y2   = 1.0,
#'   sigma_f2   = 0.4,
#'   cov_y_f    = 0.18,
#'   lambda_mode = "oracle"
#' )
#'
#'
#' ## Generate synthetic blocks for a 3-D contrast
#' p <- 3
#' cvec <- c(1, 0, -1)
#' Hpop <- diag(p)
#' SYY  <- diag(p)
#' Sff  <- 0.5 * diag(p)
#' SYf  <- 0.3 * diag(p)
#'
#' ## Vanilla λ = 1
#' n_required_ppi_pp(
#'   delta       = 0.15,
#'   N           = 2000,
#'   type        = "regression",
#'   c           = cvec,
#'   H_L         = Hpop,
#'   H_U         = Hpop,
#'   Sigma_YY    = SYY,
#'   Sigma_ff_l  = Sff,
#'   Sigma_ff_u  = Sff,
#'   Sigma_Yf    = SYf,
#'   lambda_mode = "vanilla"
#' )
#'
#' ## Oracle λ*
#' n_required_ppi_pp(
#'   delta       = 0.15,
#'   N           = 2000,
#'   type        = "regression",
#'   c           = cvec,
#'   H_L         = Hpop,
#'   H_U         = Hpop,
#'   Sigma_YY    = SYY,
#'   Sigma_ff_l  = Sff,
#'   Sigma_ff_u  = Sff,
#'   Sigma_Yf    = SYf,
#'   lambda_mode = "oracle"
#' )

n_required_ppi_pp <- function(
  delta, 
  N,
  alpha = 0.05, 
  power = 0.80,
  type = c("mean", "regression"),   # "regression" = OLS or GLM contrast
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  # mean inputs (either direct moments or metrics-based):
  sigma_y2 = NULL, 
  sigma_f2 = NULL, 
  cov_y_f  = NULL,
  var_f    = NULL, 
  var_res  = NULL, 
  metrics  = NULL,
  metric_type = NULL, 
  m_labeled = NULL, 
  correction = TRUE,  # finite-sample correction when deriving from metrics
  # regression inputs (Hessians/Fisher + Σ blocks):
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

  ## Basic validation
  if (!is.numeric(delta) || length(delta) != 1L || delta == 0)
    stop("delta must be a non-zero numeric scalar.")
  if (!is.numeric(N) || length(N) != 1L || N <= 0)
    stop("N must be a positive scalar.")

  if (warn_smallN && N < smallN_threshold)
    warning("N is small; finite-N effects may dominate.")

  ## z-values and target variance S^2
  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  S2      <- (delta / (z_alpha + z_beta))^2

  ####  MEAN ESTIMATION CASE (PPI / PPI++)
  if (type == "mean") {

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

      ## B. Metrics-based mode: derive moments from metrics/var_f/var_res
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

      if (!is.finite(sf2) || sf2 <= 0) {
        stop("Var(f) must be positive.", call. = FALSE)
      }

      ## Var(Y) from metrics: either var_y (continuous) or p_y (binary)
      if (!is.null(metrics$var_y)) {
        sy2 <- metrics$var_y
      } else if (!is.null(metrics$p_y)) {
        sy2 <- metrics$p_y * (1 - metrics$p_y)
      } else {
        stop("For mean model: metrics must include var_y (continuous) or p_y (binary).",
             call. = FALSE)
      }

      if (!is.finite(sy2) || sy2 < 0) {
        stop("Var(Y) must be non-negative.", call. = FALSE)
      }

      ## Universal identity: Var(Y-f) = Var(Y) + Var(f) - 2 Cov(Y,f)
      cyf <- (sy2 + sf2 - vres) / 2

      if (cyf^2 > sy2 * sf2 + 1e-8) {
        stop("Infeasible: |Cov(Y,f)| > sqrt(Var(Y)Var(f)).", call. = FALSE)
      }
    }

    ## C. Define variance(n) under lambda_mode
    if (lambda_mode == "vanilla") {

      ## λ = 1
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

  ####  REGRESSION CONTRAST CASE (OLS / GLM via blocks)
  if (type == "regression") {

    ## These blocks are typically produced by compute_ppi_blocks(),
    ## which itself uses compute_hessian_fisher() and compute_sigma_blocks().
    ## Here we only assume they are given.
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

      ## Oracle λ*(a) and associated variance
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

  stop("Unknown type: expected 'mean' or 'regression'.", call. = FALSE)
}