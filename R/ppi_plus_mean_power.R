ppi_pp_lambda_star <- function(n, N, cov_y_f, sigma_f2) {
  if (!is.numeric(sigma_f2) || sigma_f2 <= 0) stop("sigma_f2 must be > 0.")
  r <- n / N
  cov_y_f / ((1 + r) * sigma_f2)
}

ppi_pp_variance <- function(n,
                            N,
                            sigma_y2,
                            sigma_f2,
                            cov_y_f,
                            lambda = NULL,
                            lambda_type = c("oracle", "user")) {
  lambda_type <- match.arg(lambda_type)
  if (lambda_type == "oracle") {
    lambda <- ppi_pp_lambda_star(n, N, cov_y_f, sigma_f2)
  } else {
    if (is.null(lambda)) stop("Provide lambda when lambda_type = 'user'.")
    lambda <- as.numeric(lambda)
  }
  term_res <- sigma_y2 / n
  term_pred <- lambda^2 * sigma_f2 * (1 / N + 1 / n)
  term_cross <- -2 * lambda * cov_y_f / n
  max(term_res + term_pred + term_cross, 0)
}

resolve_ppi_pp_moments <- function(sigma_y2 = NULL,
                                   sigma_f2 = NULL,
                                   cov_y_f = NULL,
                                   var_f = NULL,
                                   var_res = NULL,
                                   metrics = NULL,
                                   metric_type = NULL,
                                   m_labeled = NULL,
                                   correction = TRUE) {
  comps <- resolve_ppi_variances(
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )
  sigma_f2 <- sigma_f2 %||% comps$var_f
  sigma_y2 <- sigma_y2 %||% metrics$var_y %||% (sigma_f2 + comps$var_res)

  cov_y_f <- cov_y_f %||% metrics$cov_y_f

  # Didn't supply cov_y_f, supplied metric_type
  if (is.null(cov_y_f) && !is.null(metric_type)) {
    metric_type_clean <- tolower(metric_type)

    if (metric_type_clean == "classification") {
      if (is.null(m_labeled) && !is.null(metrics$m_obs)) m_labeled <- metrics$m_obs
      if (is.null(m_labeled)) {
        stop("Need m_labeled (or metrics$m_obs) to recover covariances.")
      }

      p_y <- metrics$p_y
      if (is.null(p_y)) stop("classification metrics require p_y (prevalence).")

      if (!is.null(metrics$tp)) {
        tp_rate <- metrics$tp / m_labeled
        p_hat <- metrics$p_hat %||% ((metrics$tp + metrics$fp) / m_labeled)
      } else if (!is.null(metrics$precision) && !is.null(metrics$recall)) {
        tp_rate <- metrics$recall * p_y
        p_hat <- metrics$p_hat %||% (tp_rate + tp_rate * (1 / metrics$precision - 1))
      } else if (!is.null(metrics$accuracy) && !is.null(metrics$p_hat)) {
        p_hat <- metrics$p_hat
        acc <- metrics$accuracy
        tp_rate <- (p_hat + p_y + acc - 1) / 2
      } else {
        stop("classification metrics require either (tp, fp), precision/recall, or accuracy + p_hat.")
      }

      if (is.null(p_hat)) stop("Unable to infer p_hat for classification metrics.")
      cov_y_f <- tp_rate - p_y * p_hat

    } else if (metric_type_clean == "prob") {
      cov_y_f <- sigma_f2  # cross-fitted plug-in; override by passing metrics$cov_y_f if available
    }
  }

  if (is.null(cov_y_f)) stop("Unable to infer cov_y_f; pass it directly or via metrics.")

  list(
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_res = comps$var_res,
    var_f = sigma_f2
  )
}

#' Power for the PPI++ mean estimator
#'
#' @description
#' Computes two-sided normal-theory power for the PPI++ mean estimator. You can
#' supply raw variance pieces (`sigma_y2`, `sigma_f2`, `cov_y_f`, `var_f`,
#' `var_res`) or reuse the metrics interface (`metrics`, `metric_type`,
#' `m_labeled`, `correction`) to back out the required inputs.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param alpha Two-sided significance level (default 0.05).
#' @param sigma_y2 Optional outcome variance; overrides anything implied by `metrics`.
#' @param sigma_f2 Optional prediction variance; overrides anything implied by `metrics`.
#' @param cov_y_f Optional covariance \eqn{\Cov(Y, f(X))}. When supplied (directly
#'   or via `metrics$cov_y_f`) the PPI++ power is returned in addition to the PP
#'   quantities.
#' @param var_f Optional variance of \eqn{f(X)}.
#' @param var_res Optional residual variance of residual \eqn{Y - f(X)}.
#' @param metrics Optional list of predictive-performance summaries.
#' @param metric_type Character string describing the metric bundle (e.g.,
#'   `"continuous"`, `"hard"`, `"prob"`).
#' @param m_labeled Labeled sample size associated with the metrics (defaults to `n`).
#' @param correction Logical; apply finite-sample variance corrections (default `TRUE`).
#' @param lambda Optional user-specified blend weight.
#' @param lambda_type `"oracle"` (default) uses the closed-form blend; `"user"`
#'   evaluates power at the supplied `lambda`.
#'
#' @return Scalar power in \[0, 1\].
#' @export
#' 
power_ppi_pp_mean <- function(delta,
                              N,
                              n,
                              alpha = 0.05,
                              sigma_y2 = NULL,
                              sigma_f2 = NULL,
                              cov_y_f = NULL,
                              var_f = NULL,
                              var_res = NULL,
                              metrics = NULL,
                              metric_type = NULL,
                              m_labeled = n,
                              correction = TRUE,
                              lambda = NULL,
                              lambda_type = c("oracle", "user")) {
  lambda_type <- match.arg(lambda_type)
  moments <- resolve_ppi_pp_moments(
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )
  v <- ppi_pp_variance(
    n = n,
    N = N,
    sigma_y2 = moments$sigma_y2,
    sigma_f2 = moments$sigma_f2,
    cov_y_f = moments$cov_y_f,
    lambda = lambda,
    lambda_type = lambda_type
  )
  se <- sqrt(v)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}

#' Monte Carlo power for the PPI++ mean estimator
#'
#' @param R Integer, number of Monte Carlo draws.
#' @param n Labeled sample size.
#' @param N Unlabeled sample size.
#' @param alpha Two-sided test level.
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param sigma_y2 Total variance of \eqn{Y}.
#' @param sigma_f2 Variance of \eqn{f(X)}.
#' @param cov_y_f Covariance between \eqn{Y} and \eqn{f(X)}.
#' @param var_f Variance of \eqn{f(X)} (optional, overrides `sigma_f2`).
#' @param var_res Residual variance \eqn{Y - f(X)}.
#' @param lambda Optional user-specified blend weight.
#' @param lambda_type One of `"oracle"`, `"plugin"`, or `"user"`.
#' @param moments Optional list with entries for population moments.
#' @param seed RNG seed.
#' @param family A `stats::family()` object (`gaussian()` or `binomial()`).
#' @param theta0 Null hypothesis value \eqn{\theta_0}.
#'
#' @return A list with empirical and analytical power, average standard error, average lambda, and simulation details.
#' @export
simulate_power_ppiplus_mean <- function(R,
                                        n,
                                        N,
                                        alpha = 0.05,
                                        delta = NULL,
                                        sigma_y2 = NULL,
                                        sigma_f2 = NULL,
                                        cov_y_f = NULL,
                                        var_f = NULL,
                                        var_res = NULL,
                                        lambda = NULL,
                                        lambda_type = c("oracle", "plugin", "user"),
                                        moments = NULL,
                                        seed = 1,
                                        family = stats::gaussian(),
                                        theta0 = 0) {
  lambda_type <- match.arg(lambda_type)
  if (!is.null(seed)) set.seed(seed)

  # Resolve population moments (allow override via `moments`) 
  if (!is.null(moments)) {
    delta     <- moments$delta     %||% delta
    sigma_y2  <- moments$sigma_y2  %||% sigma_y2
    sigma_f2  <- moments$sigma_f2  %||% sigma_f2
    cov_y_f   <- moments$cov_y_f   %||% cov_y_f
    var_f     <- moments$var_f     %||% var_f
    var_res   <- moments$var_res   %||% var_res
    lambda    <- moments$lambda    %||% lambda
  }

  if (is.null(delta) || is.null(var_f)) {
    stop("simulate_power_ppiplus_mean requires `delta` and `var_f`.", call. = FALSE)
  }

  theta_true <- theta0 + delta

  if (identical(family$family, "binomial")) {
    if (theta_true <= 0 || theta_true >= 1) {
      stop("For binomial mean simulation, theta_true must lie in (0,1).", call. = FALSE)
    }
    # For Bernoulli with random p=f, the natural identities:
    # Var(Y) = Var(f) + E[f(1-f)], Cov(Y,f) = Var(f)
    var_res  <- var_res  %||% (theta_true * (1 - theta_true) - var_f)
    sigma_y2 <- sigma_y2 %||% (var_f + var_res)
    sigma_f2 <- sigma_f2 %||% var_f
    cov_y_f  <- cov_y_f  %||% var_f
  } else if (identical(family$family, "gaussian")) {
    if (is.null(var_res)) {
      stop("Gaussian mean simulation requires `var_res`.", call. = FALSE)
    }
    sigma_y2 <- sigma_y2 %||% (var_f + var_res)
    sigma_f2 <- sigma_f2 %||% var_f
    cov_y_f  <- cov_y_f  %||% var_f
  } else {
    stop("Unsupported family; use stats::gaussian() or stats::binomial().", call. = FALSE)
  }

  # Choose lambda used for both analytical and MC parts 
  lambda_eff <- if (lambda_type == "user") {
    if (is.null(lambda) || !is.finite(lambda)) {
      stop("Provide finite `lambda` when lambda_type = 'user'.", call. = FALSE)
    }
    lambda
  } else {
    ppi_pp_lambda_star(n = n, N = N, cov_y_f = cov_y_f, sigma_f2 = sigma_f2)
  }

  # Analytical power 
  var_ppiplus <- sigma_y2 / n +
    lambda_eff^2 * sigma_f2 * (1 / N + 1 / n) -
    2 * lambda_eff * cov_y_f / n

  var_ppiplus <- if (var_ppiplus < 0 && abs(var_ppiplus) < 1e-15) 0 else var_ppiplus
  if (var_ppiplus < 0) {
    stop("Analytical variance computed negative; check moments.", call. = FALSE)
  }

  se_ppiplus <- sqrt(var_ppiplus)
  z_alpha <- stats::qnorm(1 - alpha / 2)

  if (se_ppiplus == 0) {
    # Degenerate: reject only if delta != 0
    analytical_power <- if (delta == 0) alpha else 1.0
  } else {
    mu_ppiplus <- abs(delta) / se_ppiplus
    analytical_power <- 1 - stats::pnorm(z_alpha - mu_ppiplus) + stats::pnorm(-z_alpha - mu_ppiplus)
  }

  # Monte Carlo simulation 
  seeds <- sample.int(.Machine$integer.max, R)
  reps <- lapply(
    seq_len(R),
    function(i) simulate_ppi_mean_rep(
      n_l = n,
      n_u = N,
      family = family,
      method = "ppi_plus",
      lambda_type = lambda_type,
      lambda_user = if (lambda_type == "user") lambda_eff else NULL,
      alpha = alpha,
      theta0 = theta0,
      delta = delta,                     # DGP mean shift
      var_f = var_f,                     # DGP variance of f
      var_res = var_res,                 # DGP residual variance (Gaussian) or Bernoulli implied
      seed = seeds[i]
    )
  )

  empirical <- mean(vapply(reps, function(x) x$reject_ppplus, logical(1)))
  avg_SE <- mean(vapply(reps, function(x) x$se_ppplus,    numeric(1)))
  lambda_vals <- vapply(reps, function(x) x$lambda,        numeric(1))
  avg_lambda <- mean(lambda_vals, na.rm = TRUE)
  if (!is.finite(avg_lambda)) avg_lambda <- NA_real_

  list(
    empirical_power = empirical,
    analytical_power = analytical_power,
    avg_SE = avg_SE,
    avg_lambda = avg_lambda,
    delta_used = delta,
    details = reps
  )
}