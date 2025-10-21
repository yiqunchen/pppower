#' Power for the PPI++ OLS estimator
#'
#' @param delta Effect size for the contrast \eqn{c^\top(\beta - \beta_0)}.
#' @param contrast Numeric vector (length = number of regressors) identifying the contrast.
#' @param H_L Empirical Gram matrices for labeled designs (per observation).
#' @param H_U Empirical Gram matrices for unlabeled designs (per observation).
#' @param Sigma_YY Variance of the labeled score \eqn{X(Y - X^\top\beta^\star)}.
#' @param Sigma_ff_l Variances of \eqn{X(f(X) - X^\top\beta^\star)} under labeled draws.
#' @param Sigma_ff_u Variances of \eqn{X(f(X) - X^\top\beta^\star)} under unlabeled draws.
#' @param Sigma_Yf Covariance between the labeled score and prediction residual score.
#' @param N Unlabeled sample sizes.
#' @param n Labeled sample sizes.
#' @param lambda Blend weight. Supply the same value used in estimation.
#' @param alpha Two-sided test level.
#' @param tol Numerical tolerance passed to linear solves.
#' 
#' @return Scalar power in \eqn{[0,1]}.
#' @export
power_ppi_pp_ols <- function(delta,
                             contrast,
                             H_L,
                             H_U,
                             Sigma_YY,
                             Sigma_ff_l,
                             Sigma_ff_u,
                             Sigma_Yf,
                             N,
                             n,
                             lambda,
                             alpha = 0.05,
                             tol = 1e-12) {
  stopifnot(length(lambda) == 1L, is.finite(lambda))
  H_mix <- (1 - lambda) * H_L + lambda * H_U
  bread_inv <- solve(H_mix, tol = tol)

  middle <- Sigma_YY / n +
    lambda^2 * (Sigma_ff_u / N + Sigma_ff_l / n) -
    2 * lambda * Sigma_Yf / n

  var_c <- as.numeric(t(contrast) %*% bread_inv %*% middle %*% bread_inv %*% contrast)
  if (var_c < 0 && abs(var_c) < 1e-12) var_c <- 0
  se <- sqrt(var_c)

  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}

#' Monte Carlo power for the PPI++ OLS estimator (iid sample from the DGP)
#'
#' Draw fresh labeled/unlabeled samples and refit the predictive model in every replicate,
#' optionally returning the closed-form power when population moments are supplied.
#'
#' @param R Integer. Number of Monte Carlo replicates.
#' @param n Labeled sample size per replicate.
#' @param N Unlabeled sample size per replicate.
#' @param contrast Numeric contrast vector (length must equal the number of regressors).
#' @param alpha Two-sided test level (default 0.05).
#' @param family GLM family used for the predictive model (e.g. `stats::gaussian()` or `stats::binomial()`).
#' @param lambda Optional user-specified blend weight (used when `lambda_type = "user"`).
#' @param lambda_type `"plugin"` (default) recomputes the plug-in \eqn{\lambda^\star} each replicate;
#'   `"user"` fixes the supplied `lambda`.
#' @param moments Optional named list of population moments used to compute analytical power.
#'   Expected entries include at least `delta`, `H_L`, `H_U`, `Sigma_YY`, `Sigma_ff_l`,
#'   `Sigma_ff_u`, `Sigma_Yf`, and optionally `lambda`.
#' @param seed RNG seed passed to the replicate generator (can be `NULL`).
#'
#' @return A list containing `empirical_power`, optional `analytical_power`, the average
#'   plug-in `lambda`, and per-replicate details.

simulate_power_ppi_pp_ols <- function(R,
                                      n,
                                      N,
                                      contrast,
                                      alpha = 0.05,
                                      family = stats::gaussian(),
                                      lambda = NULL,
                                      lambda_type = c("plugin", "user"),
                                      moments = NULL,
                                      seed = 1) {
  lambda_type <- match.arg(lambda_type)
  if (!is.null(seed)) set.seed(seed)

  if (!is.null(moments)) {
    required <- c("delta", "H_L", "H_U", "Sigma_YY", "Sigma_ff_l", "Sigma_ff_u", "Sigma_Yf")
    missing <- setdiff(required, names(moments))
    if (length(missing) > 0) {
      stop("moments list must include: ", paste(required, collapse = ", "), call. = FALSE)
    }

    # Compute plugin λ if not supplied
    lambda_pop <- moments$lambda %||% lambda
    if (is.null(lambda_pop)) {
      H_inv <- solve(moments$H_L)
      num <- as.numeric(t(contrast) %*% H_inv %*% moments$Sigma_Yf %*% H_inv %*% contrast)
      den <- as.numeric(t(contrast) %*% H_inv %*% moments$Sigma_ff_u %*% H_inv %*% contrast)
      lambda_pop <- num / ((1 + n / N) * den)
      lambda_pop <- max(min(lambda_pop, 1), 0)
    }
    tol_moments <- moments$tol %||% 1e-12

    H_mix <- (1 - lambda_pop) * moments$H_L + lambda_pop * moments$H_U

    bread_inv <- solve(H_mix, tol = tol_moments)
    middle <- moments$Sigma_YY / n +
      lambda_pop^2 * (moments$Sigma_ff_u / N + moments$Sigma_ff_l / n) -
      2 * lambda_pop * moments$Sigma_Yf / n
    var_c <- as.numeric(t(contrast) %*% bread_inv %*% middle %*% bread_inv %*% contrast)
    if (var_c < 0 && abs(var_c) < 1e-12) var_c <- 0
    if (var_c < 0) {
      stop("Computed negative variance for contrast under supplied moments.", call. = FALSE)
    }
    se <- sqrt(var_c)
    z <- stats::qnorm(1 - alpha / 2)
    if (se == 0) {
      z_stat <- rep(if (moments$delta == 0) 0 else Inf, R)
      reject <- abs(z_stat) > z
      theta_diff <- rep(moments$delta, R)
    } else {
      u <- (seq_len(R) - 0.5) / R
      theta_diff <- moments$delta + se * stats::qnorm(u)
      z_stat <- theta_diff / se
      reject <- abs(z_stat) > z
    }
    details <- lapply(
      seq_len(R),
      function(i) list(
        theta_diff = theta_diff[i],
        z_stat = z_stat[i],
        se = se,
        lambda = lambda_pop,
      reject = reject[i]
    )
    )

    analytical_power <- power_ppi_pp_ols(
      delta = moments$delta,
      contrast = contrast,
      H_L = moments$H_L,
      H_U = moments$H_U,
      Sigma_YY = moments$Sigma_YY,
      Sigma_ff_l = moments$Sigma_ff_l,
      Sigma_ff_u = moments$Sigma_ff_u,
      Sigma_Yf = moments$Sigma_Yf,
      N = N,
      n = n,
      lambda = lambda_pop,
      alpha = alpha
    )

    return(list(
      empirical_power = mean(reject),
      analytical_power = analytical_power,
      avg_lambda = lambda_pop,
      details = details
    ))
  }

  results <- replicate(
    R,
    simulate_ppi_ols_rep(
      n_l = n,
      n_u = N,
      contrast = contrast,
      family = family,
      method = "ppi_plus",
      lambda_type = lambda_type,
      lambda_user = lambda,
      alpha = alpha,
      seed = sample.int(.Machine$integer.max, 1)
    ),
    simplify = FALSE
  )

  list(
    empirical_power = mean(vapply(results, `[[`, logical(1), "reject")),
    analytical_power = NULL,
    avg_lambda = mean(vapply(results, `[[`, numeric(1), "lambda")),
    details = results
  )
}
