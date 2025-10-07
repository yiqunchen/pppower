#' Required labeled sample size for desired PPI power
#'
#' @description
#' Computes the minimum labeled sample size \eqn{n} needed to attain a target
#' two-sided power when using prediction-powered inference, covering both the
#' mean estimator and linear-regression (OLS) contrasts.
#'
#' @param delta Effect size \eqn{\theta - \theta_0} for the estimand of interest.
#' @param N Unlabeled sample size.
#' @param alpha Two-sided significance level.
#' @param power Target power (\eqn{1-\beta}).
#' @param type Either \code{"mean"} for the PP mean estimator, \code{"ols"} for a linear
#'   contrast in PP-OLS, or \code{"custom"} to supply pre-computed variance pieces.
#' @param var_f Variance of \eqn{f(X)}. Required when \code{type = "mean"}.
#' @param var_res Variance of \eqn{Y - f(X)}. Required when \code{type = "mean"}.
#' @param V_u,V_l Sandwich covariance matrices for the unlabeled and labeled
#'   PP-OLS pieces. Required when \code{type = "ols"}.
#' @param c Numeric vector specifying the linear contrast \eqn{c^\top \beta}
#'   for PP-OLS. Required when \code{type = "ols"}.
#' @param var_unlabeled,var_labeled Optional scalar variance components when
#'   \code{type = "custom"}. These correspond to \eqn{\mathrm{Var}_u} and
#'   \eqn{\mathrm{Var}_l} so that the standard error is
#'   \eqn{\sqrt{\mathrm{Var}_u/N  \mathrm{Var}_l/n}}.
#' @param warn_smallN Logical; warn if \code{N} is below \code{smallN_threshold}.
#' @param smallN_threshold Threshold for triggering the small-\code{N} warning.
#' @param mode Either \code{"error"} (default) to throw when the required \eqn{n}
#'   exceeds \code{N}, or \code{"cap"} to return \code{n = N} with the achieved power
#'   attached as an attribute.
#'
#' @return Integer labeled sample size meeting the target power, possibly with an
#'   \code{"achieved_power"} attribute when \code{mode = "cap"}.
#'
#' @examples
#' # Mean estimator
#' n_required_PP(
#'   delta = 0.1, N = 5000, var_f = 0.2, var_res = 0.5,
#'   alpha = 0.05, power = 0.8, type = "mean"
#' )
#'
#' # PP-OLS contrast
#' V_u <- diag(c(0.8, 0.6)); V_l <- diag(c(0.5, 0.7))
#' n_required_PP(
#'   delta = 0.15, N = 4000, V_u = V_u, V_l = V_l, c = c(1, 0),
#'   alpha = 0.05, power = 0.85, type = "ols"
#' )
#' @export
n_required_PP <- function(delta,
                          N,
                          alpha = 0.05,
                          power = 0.80,
                          type = c("mean", "ols", "custom"),
                          var_f = NULL,
                          var_res = NULL,
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

  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta)) {
    stop("delta must be a finite numeric scalar.", call. = FALSE)
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

  if (warn_smallN && is.finite(N) && N < smallN_threshold) {
    warning(
      "Unlabeled N is quite small. The finite-N term may dominate the variance.",
      call. = FALSE
    )
  }

  components <- switch(
    type,
    mean = {
      if (is.null(var_f) || is.null(var_res)) {
        stop("var_f and var_res must be supplied when type = 'mean'.", call. = FALSE)
      }
      if (var_f < 0 || var_res <= 0) {
        stop("var_f must be >= 0 and var_res must be > 0.", call. = FALSE)
      }
      list(
        var_u = as.numeric(var_f),
        var_l = as.numeric(var_res),
        power_fn = function(n) power_ppi_mean(delta, var_f, var_res, N, n, alpha)
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
      if (var_u < 0) {
        stop("c' V_u c must be non-negative.", call. = FALSE)
      }
      if (var_l <= 0) {
        stop("c' V_l c must be positive.", call. = FALSE)
      }
      list(
        var_u = var_u,
        var_l = var_l,
        power_fn = function(n) power_ppi_ols(delta, V_u, V_l, N, n, c_vec, alpha)
      )
    },
    custom = {
      if (is.null(var_unlabeled) || is.null(var_labeled)) {
        stop("var_unlabeled and var_labeled must be supplied when type = 'custom'.",
             call. = FALSE)
      }
      var_u <- as.numeric(var_unlabeled)
      var_l <- as.numeric(var_labeled)
      if (var_u < 0) {
        stop("var_unlabeled must be non-negative.", call. = FALSE)
      }
      if (var_l <= 0) {
        stop("var_labeled must be positive.", call. = FALSE)
      }
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
