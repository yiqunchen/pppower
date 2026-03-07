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

# Takes in type of metric and statistics of the metrics, then output variances and covariances
# for downstream functions.
resolve_ppi_pp_moments <- function(sigma_y2 = NULL,
                                   sigma_f2 = NULL,
                                   cov_y_f = NULL,
                                   var_f = NULL,
                                   var_res = NULL,
                                   metrics = NULL,
                                   metric_type = NULL,
                                   m_labeled = NULL,
                                   correction = TRUE) {
  # Fast path: if all three direct moments are provided, skip resolve_ppi_variances
  if (!is.null(sigma_y2) && !is.null(sigma_f2) && !is.null(cov_y_f)) {
    return(list(
      sigma_y2 = as.numeric(sigma_y2),
      sigma_f2 = as.numeric(sigma_f2),
      cov_y_f = as.numeric(cov_y_f),
      var_res = NA_real_,
      var_f = as.numeric(sigma_f2)
    ))
  }

  # Fast path: if var_f, var_res, and cov_y_f are provided
  if (!is.null(var_f) && !is.null(var_res) && !is.null(cov_y_f)) {
    sf2 <- as.numeric(var_f)
    sy2 <- sigma_y2 %||% (sf2 + as.numeric(var_res))
    return(list(
      sigma_y2 = as.numeric(sy2),
      sigma_f2 = sf2,
      cov_y_f = as.numeric(cov_y_f),
      var_res = as.numeric(var_res),
      var_f = sf2
    ))
  }

  comps <- resolve_ppi_variances(
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )
  sigma_f2 <- sigma_f2 %||% comps$var_f

  # For binary outcomes (classification), Var(Y) = p_y * (1 - p_y)
  # The formula sigma_y2 = var_f + var_res only holds when Y and f are independent
  if (!is.null(metric_type) && tolower(metric_type) == "classification" && !is.null(metrics$p_y)) {
    sigma_y2 <- sigma_y2 %||% metrics$var_y %||% (metrics$p_y * (1 - metrics$p_y))
  } else {
    sigma_y2 <- sigma_y2 %||% metrics$var_y %||% (sigma_f2 + comps$var_res)
  }

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
      } else if (!is.null(metrics$sensitivity) && !is.null(metrics$specificity)) {
        # Sensitivity/specificity path
        tp_rate <- metrics$sensitivity * p_y
        fp_rate <- (1 - metrics$specificity) * (1 - p_y)
        p_hat <- metrics$p_hat %||% (tp_rate + fp_rate)
      } else if (!is.null(metrics$precision) && !is.null(metrics$recall)) {
        tp_rate <- metrics$recall * p_y
        p_hat <- metrics$p_hat %||% (tp_rate + tp_rate * (1 / metrics$precision - 1))
      } else if (!is.null(metrics$accuracy) && !is.null(metrics$p_hat)) {
        p_hat <- metrics$p_hat
        acc <- metrics$accuracy
        tp_rate <- (p_hat + p_y + acc - 1) / 2
      } else {
        stop("classification metrics require either (tp, fp), sensitivity/specificity, precision/recall, or accuracy + p_hat.")
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

# ============================================================================
# Unified power_ppi_mean(): pwr-style "leave one NULL" function
# ============================================================================

#' Power or required sample size for the PPI++ mean estimator
#'
#' @description
#' Unified interface for power analysis of the PPI++ mean estimator.
#' Follows the \pkg{pwr} convention: supply **either** `n` or `power` (but not
#' both). The function computes whichever is left `NULL`.
#'
#' * `power = NULL` (default): compute power given `n`.
#' * `n = NULL`: solve for the minimum labeled sample size achieving `power`.
#'
#' You can supply raw variance pieces (`sigma_y2`, `sigma_f2`, `cov_y_f`,
#' `var_f`, `var_res`) or reuse the metrics interface (`metrics`,
#' `metric_type`, `m_labeled`, `correction`) to back out the required inputs.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size. Set to `NULL` to solve for the required `n`.
#' @param power Target power. Set to `NULL` (default) to compute power from `n`.
#' @param alpha Two-sided significance level (default 0.05).
#' @param sigma_y2 Optional outcome variance; overrides anything implied by `metrics`.
#' @param sigma_f2 Optional prediction variance; overrides anything implied by `metrics`.
#' @param cov_y_f Optional covariance \eqn{\mathrm{Cov}(Y, f(X))}. When supplied (directly
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
#' @param lambda_mode Alias for `lambda_type` used in sample-size mode;
#'   one of `"vanilla"`, `"oracle"`, or `"user"`.
#' @param lambda_user Scalar specifying \eqn{\lambda} when
#'   `lambda_mode = "user"` (sample-size mode only).
#' @param mode `"error"` (default) or `"cap"`.
#'   `"cap"` returns `n = N` with an `"achieved_power"` attribute when the
#'   required sample size exceeds `N`. Only used when solving for `n`.
#' @param warn_smallN Logical. Warn when `N < smallN_threshold` (default `TRUE`).
#' @param smallN_threshold Numeric. Default `500`.
#'
#' @return When computing power: scalar in \[0, 1\].
#'   When computing sample size: integer `n` required, possibly with
#'   `"achieved_power"` attribute if capped.
#' @export
#'
#' @examples
#' # Compute power given n:
#' power_ppi_mean(delta = 0.3, N = 5000, n = 200,
#'                sigma_y2 = 1, sigma_f2 = 0.5, cov_y_f = 0.4)
#'
#' # Compute required n given target power:
#' power_ppi_mean(delta = 0.3, N = 5000, power = 0.8,
#'                sigma_y2 = 1, sigma_f2 = 0.5, cov_y_f = 0.4)
#'
power_ppi_mean <- function(delta,
                           N,
                           n = NULL,
                           power = NULL,
                           alpha = 0.05,
                           sigma_y2 = NULL,
                           sigma_f2 = NULL,
                           cov_y_f = NULL,
                           var_f = NULL,
                           var_res = NULL,
                           metrics = NULL,
                           metric_type = NULL,
                           m_labeled = NULL,
                           correction = TRUE,
                           lambda = NULL,
                           lambda_type = c("oracle", "user"),
                           lambda_mode = NULL,
                           lambda_user = NULL,
                           mode = c("error", "cap"),
                           warn_smallN = TRUE,
                           smallN_threshold = 500) {
  # Validate: exactly one of n, power must be NULL
  if (is.null(n) && is.null(power)) {
    stop("Exactly one of 'n' and 'power' must be NULL. Both are NULL.", call. = FALSE)
  }
  if (!is.null(n) && !is.null(power)) {
    stop("Exactly one of 'n' and 'power' must be NULL. Neither is NULL.", call. = FALSE)
  }

  if (!is.null(n)) {
    # ---- Compute power given n ----
    lambda_type <- match.arg(lambda_type)
    if (is.null(m_labeled)) m_labeled <- n
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
    return(1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu))

  } else {
    # ---- Solve for n given target power ----
    # Determine lambda_mode from lambda_type or explicit argument
    if (!is.null(lambda_mode)) {
      lm <- match.arg(lambda_mode, c("vanilla", "oracle", "user"))
    } else if (!is.null(lambda) || (!missing(lambda_type) && match.arg(lambda_type) == "user")) {
      lm <- "user"
    } else {
      lm <- "oracle"
    }

    n_required_ppi_pp_mean_internal(
      delta = delta,
      N = N,
      alpha = alpha,
      power = power,
      lambda_mode = lm,
      lambda_user = lambda_user %||% lambda,
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f,
      var_f = var_f,
      var_res = var_res,
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = m_labeled,
      correction = correction,
      warn_smallN = warn_smallN,
      smallN_threshold = smallN_threshold,
      mode = match.arg(mode)
    )
  }
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
simulate_ppi_mean <- function(R,
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
    stop("simulate_ppi_mean requires `delta` and `var_f`.", call. = FALSE)
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

#' Power for PPI++ paired t-test (within-subject differences)
#'
#' @description
#' Computes power for a paired t-test using PPI++. This is equivalent to
#' mean estimation on the differences \eqn{D_i = Y_i^A - Y_i^B} with auxiliary
#' predictions \eqn{f^D_i = f^A(X_i) - f^B(X_i)}.
#'
#' @param delta Effect size: the true mean difference \eqn{E[Y^A - Y^B]}.
#' @param N Number of unlabeled pairs with predictions.
#' @param n Number of labeled pairs.
#' @param power Target power. Set to `NULL` (default) to compute power from `n`.
#'   Set to a value in (0, 1) to solve for `n` instead.
#' @param alpha Two-sided significance level (default 0.05).
#' @param sigma_D2 Variance of the differences \eqn{D = Y^A - Y^B}.
#' @param rho_D Correlation between \eqn{D} and \eqn{f^D}.
#' @param sigma_fD2 Optional variance of predicted differences. If NULL,
#'   assumed equal to \code{sigma_D2}.
#'
#' @return When `power = NULL`: scalar power in \eqn{[0, 1]}.
#'   When `n = NULL`: required number of labeled pairs (integer).
#'
#' @details
#' The paired t-test reduces to a one-sample problem on the differences.
#' The PPI++ variance for paired data is:
#' \deqn{\mathrm{Var}(\hat\Delta_{\lambda^*}) \approx \frac{\sigma_D^2(1 - \rho_D^2)}{n}}
#' where \eqn{\rho_D} is the correlation between the observed difference and
#' the predicted difference.
#'
#' @examples
#' # Compute power
#' power_ppi_paired(
#'   delta = 0.3,
#'   N = 1000,
#'   n = 50,
#'   sigma_D2 = 1,
#'   rho_D = 0.7
#' )
#'
#' # Compute required sample size
#' power_ppi_paired(
#'   delta = 0.3,
#'   N = 1000,
#'   power = 0.80,
#'   sigma_D2 = 1,
#'   rho_D = 0.7
#' )
#'
#' @export
power_ppi_paired <- function(delta,
                              N,
                              n = NULL,
                              power = NULL,
                              alpha = 0.05,
                              sigma_D2,
                              rho_D,
                              sigma_fD2 = NULL) {
  if (is.null(sigma_fD2)) {
    sigma_fD2 <- sigma_D2
  }

  # Covariance between D and f^D
  cov_D_fD <- rho_D * sqrt(sigma_D2 * sigma_fD2)

  # Delegate to unified power_ppi_mean
  power_ppi_mean(
    delta = delta,
    N = N,
    n = n,
    power = power,
    alpha = alpha,
    sigma_y2 = sigma_D2,
    sigma_f2 = sigma_fD2,
    cov_y_f = cov_D_fD
  )
}

#' Power for PPI++ paired test with binary outcomes
#'
#' @description
#' Computes power for a paired test when outcomes are binary and predictions are
#' derived from a classifier with sensitivity and specificity. The paired
#' differences \eqn{D = Y^A - Y^B} and prediction differences \eqn{f^D = f^A - f^B}
#' are treated via a CLT approximation using population moments implied by a
#' joint Bernoulli model with within-pair correlation.
#'
#' @param delta Effect size: mean difference \eqn{E[Y^A - Y^B]}. If NULL, uses
#'   \eqn{p_A - p_B}.
#' @param N Number of unlabeled pairs with predictions.
#' @param n Number of labeled pairs.
#' @param alpha Two-sided significance level (default 0.05).
#' @param p_A Marginal probability \eqn{P(Y^A=1)}.
#' @param p_B Marginal probability \eqn{P(Y^B=1)}.
#' @param rho_within Pearson correlation between \eqn{Y^A} and \eqn{Y^B}.
#' @param sens Classifier sensitivity.
#' @param spec Classifier specificity.
#'
#' @return Scalar power in \eqn{[0, 1]}.
#'
#' @details
#' The joint probability \eqn{P(Y^A=1, Y^B=1)} is set to
#' \deqn{p_{11} = p_A p_B + \rho_{\text{within}} \sqrt{p_A(1-p_A)p_B(1-p_B)}}
#' and clipped to the feasible range. Predictions are generated conditionally
#' on \eqn{Y} with sensitivity/specificity, and population moments are computed
#' by marginalizing over the joint distribution of \eqn{(Y^A, Y^B)}.
#'
#' @examples
#' power_ppi_paired_binary(
#'   delta = 0.08,
#'   N = 5000,
#'   n = 500,
#'   p_A = 0.34,
#'   p_B = 0.26,
#'   rho_within = 0.3,
#'   sens = 0.85,
#'   spec = 0.85
#' )
#'
#' @export
power_ppi_paired_binary <- function(delta,
                                     N,
                                     n,
                                     alpha = 0.05,
                                     p_A,
                                     p_B,
                                     rho_within,
                                     sens,
                                     spec) {
  if (is.null(delta)) {
    delta <- p_A - p_B
  } else if (abs(delta - (p_A - p_B)) > 1e-8) {
    warning("delta does not match p_A - p_B; using provided delta.")
  }

  denom <- sqrt(p_A * (1 - p_A) * p_B * (1 - p_B))
  if (denom == 0) {
    p11 <- min(p_A, p_B)
  } else {
    p11 <- p_A * p_B + rho_within * denom
  }
  lower <- max(0, p_A + p_B - 1)
  upper <- min(p_A, p_B)
  p11 <- min(max(p11, lower), upper)
  p10 <- p_A - p11
  p01 <- p_B - p11
  p00 <- 1 - p11 - p10 - p01
  probs <- c(p11, p10, p01, p00)
  if (any(probs < -1e-10)) {
    stop("Invalid joint probabilities for (Y_A, Y_B).")
  }
  probs <- pmax(probs, 0)
  probs <- probs / sum(probs)
  p11 <- probs[1]; p10 <- probs[2]; p01 <- probs[3]; p00 <- probs[4]

  mean_D <- p10 - p01
  mean_D2 <- p10 + p01
  sigma_D2 <- mean_D2 - mean_D^2

  q1 <- sens
  q0 <- 1 - spec
  calc_fD_moments <- function(yA, yB) {
    qA <- if (yA == 1) q1 else q0
    qB <- if (yB == 1) q1 else q0
    mean_fD <- qA - qB
    var_fD <- qA * (1 - qA) + qB * (1 - qB)
    list(mean = mean_fD, mean2 = var_fD + mean_fD^2)
  }

  m11 <- calc_fD_moments(1, 1)
  m10 <- calc_fD_moments(1, 0)
  m01 <- calc_fD_moments(0, 1)
  m00 <- calc_fD_moments(0, 0)

  mean_fD <- p11 * m11$mean + p10 * m10$mean + p01 * m01$mean + p00 * m00$mean
  mean_fD2 <- p11 * m11$mean2 + p10 * m10$mean2 + p01 * m01$mean2 + p00 * m00$mean2
  sigma_fD2 <- mean_fD2 - mean_fD^2

  mean_DfD <- p10 * (1 * m10$mean) + p01 * (-1 * m01$mean)
  cov_D_fD <- mean_DfD - mean_D * mean_fD

  power_ppi_mean(
    delta = delta,
    N = N,
    n = n,
    alpha = alpha,
    sigma_y2 = sigma_D2,
    sigma_f2 = sigma_fD2,
    cov_y_f = cov_D_fD
  )
}

#' Power for PPI++ two-sample t-test (difference in means)
#'
#' @description
#' Computes power for a two-sample t-test using PPI++ where we estimate
#' the difference in means between two independent groups.
#'
#' @param delta Effect size: the true difference in means \eqn{\mu_A - \mu_B}.
#' @param n_A Labeled sample size for group A.
#' @param n_B Labeled sample size for group B.
#' @param N_A Unlabeled sample size for group A.
#' @param N_B Unlabeled sample size for group B.
#' @param alpha Two-sided significance level (default 0.05).
#' @param sigma_y2_A Variance of outcomes in group A.
#' @param sigma_f2_A Variance of predictions in group A.
#' @param cov_yf_A Covariance between Y and f in group A.
#' @param sigma_y2_B Variance of outcomes in group B.
#' @param sigma_f2_B Variance of predictions in group B.
#' @param cov_yf_B Covariance between Y and f in group B.
#'
#' @return Scalar power in \eqn{[0, 1]}.
#'
#' @details
#' The variance of the difference estimator is the sum of the PPI++ variances
#' from each group (assuming independence):
#' \deqn{\mathrm{Var}(\hat\mu_A - \hat\mu_B) = \mathrm{Var}_{PPI++}(\hat\mu_A) + \mathrm{Var}_{PPI++}(\hat\mu_B)}
#'
#' @examples
#' # Two-sample test with equal groups
#' power_ppi_ttest(
#'   delta = 0.3,
#'   n_A = 100, n_B = 100,
#'   N_A = 5000, N_B = 5000,
#'   sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = 0.7,
#'   sigma_y2_B = 1, sigma_f2_B = 1, cov_yf_B = 0.7
#' )
#'
#' @export
power_ppi_ttest <- function(delta,
                             n_A, n_B,
                             N_A, N_B,
                             alpha = 0.05,
                             sigma_y2_A, sigma_f2_A, cov_yf_A,
                             sigma_y2_B, sigma_f2_B, cov_yf_B) {
  # Oracle lambda for each group
  lambda_A <- cov_yf_A / ((1 + n_A / N_A) * sigma_f2_A)
  lambda_B <- cov_yf_B / ((1 + n_B / N_B) * sigma_f2_B)

  # PPI++ variance for each group
  var_A <- sigma_y2_A / n_A +
    lambda_A^2 * sigma_f2_A * (1 / N_A + 1 / n_A) -
    2 * lambda_A * cov_yf_A / n_A

  var_B <- sigma_y2_B / n_B +
    lambda_B^2 * sigma_f2_B * (1 / N_B + 1 / n_B) -
    2 * lambda_B * cov_yf_B / n_B

  # Total variance of difference (independent groups)
  var_total <- max(var_A, 0) + max(var_B, 0)
  se <- sqrt(var_total)

  # Power calculation
  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu <- abs(delta) / se
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}

#' Power for PPI++ two-sample test with binary outcomes (sens/spec input)
#'
#' @description
#' Computes PPI++ power for a two-sample difference in proportions when
#' predictions come from a binary classifier summarized by sensitivity and
#' specificity. This helper converts the classification metrics into the
#' required moment inputs and delegates to `power_ppi_ttest()`.
#'
#' @param p_A,p_B Group-specific prevalences \eqn{P(Y=1)}.
#' @param n_A,n_B Labeled sample sizes per group.
#' @param N_A,N_B Unlabeled sample sizes per group.
#' @param sens_A,spec_A Sensitivity and specificity for group A predictions.
#' @param sens_B,spec_B Sensitivity and specificity for group B predictions
#'   (defaults match group A).
#' @param delta Optional effect size \eqn{p_A - p_B}. If `NULL`, it is
#'   computed from `p_A` and `p_B`.
#' @param alpha Two-sided significance level (default 0.05).
#'
#' @return Scalar power in \[0, 1\] with a `"moments"` attribute containing
#'   the per-group moment lists.
#' @export
power_ppi_ttest_binary <- function(p_A,
                                    p_B,
                                    n_A,
                                    n_B,
                                    N_A,
                                    N_B,
                                    sens_A,
                                    spec_A,
                                    sens_B = sens_A,
                                    spec_B = spec_A,
                                    delta = NULL,
                                    alpha = 0.05) {
  if (is.null(delta)) {
    delta <- p_A - p_B
  } else if (abs(delta - (p_A - p_B)) > 1e-8) {
    warning("delta does not equal p_A - p_B; using supplied delta.")
  }

  moments_A <- binary_moments_from_sens_spec(p = p_A, sens = sens_A, spec = spec_A)
  moments_B <- binary_moments_from_sens_spec(p = p_B, sens = sens_B, spec = spec_B)

  power <- power_ppi_ttest(
    delta = delta,
    n_A = n_A, n_B = n_B,
    N_A = N_A, N_B = N_B,
    alpha = alpha,
    sigma_y2_A = moments_A$sigma_y2,
    sigma_f2_A = moments_A$sigma_f2,
    cov_yf_A = moments_A$cov_y_f,
    sigma_y2_B = moments_B$sigma_y2,
    sigma_f2_B = moments_B$sigma_f2,
    cov_yf_B = moments_B$cov_y_f
  )
  attr(power, "moments") <- list(A = moments_A, B = moments_B)
  power
}

#' Monte Carlo power for binary two-sample PPI++ test
#'
#' @param R Number of Monte Carlo replicates.
#' @param p_A,p_B Group-specific prevalences \eqn{P(Y=1)}.
#' @param sens_A,spec_A Sensitivity/specificity for group A predictions.
#' @param sens_B,spec_B Sensitivity/specificity for group B predictions
#'   (defaults match group A).
#' @param n_A,n_B Labeled sample sizes.
#' @param N_A,N_B Unlabeled sample sizes.
#' @param alpha Two-sided significance level.
#' @param delta Optional effect size (defaults to `p_A - p_B`).
#' @param seed RNG seed (set to `NULL` to avoid seeding).
#'
#' @return List with empirical and analytical power plus simulation details.
#' @export
simulate_ppi_ttest_binary <- function(R,
                                      p_A,
                                      p_B,
                                      sens_A,
                                      spec_A,
                                      n_A,
                                      n_B,
                                      N_A,
                                      N_B,
                                      sens_B = sens_A,
                                      spec_B = spec_A,
                                      alpha = 0.05,
                                      delta = NULL,
                                      seed = 1) {
  if (!is.null(seed)) set.seed(seed)

  moments_A <- binary_moments_from_sens_spec(p = p_A, sens = sens_A, spec = spec_A)
  moments_B <- binary_moments_from_sens_spec(p = p_B, sens = sens_B, spec = spec_B)

  if (is.null(delta)) {
    delta <- p_A - p_B
  } else if (abs(delta - (p_A - p_B)) > 1e-8) {
    warning("delta does not equal p_A - p_B; using supplied delta.")
  }

  lambda_A <- ppi_pp_lambda_star(n = n_A, N = N_A,
                                 cov_y_f = moments_A$cov_y_f,
                                 sigma_f2 = moments_A$sigma_f2)
  lambda_B <- ppi_pp_lambda_star(n = n_B, N = N_B,
                                 cov_y_f = moments_B$cov_y_f,
                                 sigma_f2 = moments_B$sigma_f2)

  analytical <- power_ppi_ttest(
    delta = delta,
    n_A = n_A, n_B = n_B,
    N_A = N_A, N_B = N_B,
    alpha = alpha,
    sigma_y2_A = moments_A$sigma_y2,
    sigma_f2_A = moments_A$sigma_f2,
    cov_yf_A = moments_A$cov_y_f,
    sigma_y2_B = moments_B$sigma_y2,
    sigma_f2_B = moments_B$sigma_f2,
    cov_yf_B = moments_B$cov_y_f
  )

  seeds <- sample.int(.Machine$integer.max, R)
  rejects <- vapply(
    seq_len(R),
    function(i) {
      if (!is.null(seed)) set.seed(seeds[i])

      Y_A <- stats::rbinom(n_A, 1, p_A)
      f_A_L <- ifelse(
        Y_A == 1,
        stats::rbinom(n_A, 1, sens_A),
        stats::rbinom(n_A, 1, 1 - spec_A)
      )
      Y_A_U <- stats::rbinom(N_A, 1, p_A)
      f_A_U <- ifelse(
        Y_A_U == 1,
        stats::rbinom(N_A, 1, sens_A),
        stats::rbinom(N_A, 1, 1 - spec_A)
      )

      Y_B <- stats::rbinom(n_B, 1, p_B)
      f_B_L <- ifelse(
        Y_B == 1,
        stats::rbinom(n_B, 1, sens_B),
        stats::rbinom(n_B, 1, 1 - spec_B)
      )
      Y_B_U <- stats::rbinom(N_B, 1, p_B)
      f_B_U <- ifelse(
        Y_B_U == 1,
        stats::rbinom(N_B, 1, sens_B),
        stats::rbinom(N_B, 1, 1 - spec_B)
      )

      res <- ppi_ttest(
        Y_A, f_A_L, f_A_U,
        Y_B, f_B_L, f_B_U,
        delta0 = 0,
        alpha = alpha,
        lambda_A = lambda_A,
        lambda_B = lambda_B
      )
      res$reject
    },
    logical(1)
  )

  empirical <- mean(rejects)
  mc_se <- sqrt(analytical * (1 - analytical) / R)

  list(
    empirical_power = empirical,
    analytical_power = analytical,
    mc_se = mc_se,
    lambda_A = lambda_A,
    lambda_B = lambda_B,
    moments_A = moments_A,
    moments_B = moments_B,
    delta_used = delta,
    R = R
  )
}


# ============================================================================
# Deprecated shims (old names → new names)
# ============================================================================

#' @rdname deprecated
#' @export
power_ppi_pp_mean <- function(delta, N, n, ...) {
  .Deprecated("power_ppi_mean", package = "pppower")
  power_ppi_mean(delta = delta, N = N, n = n, ...)
}

#' @rdname deprecated
#' @export
power_ppi_pp_paired <- function(...) {
  .Deprecated("power_ppi_paired", package = "pppower")
  power_ppi_paired(...)
}

#' @rdname deprecated
#' @export
power_ppi_pp_paired_binary <- function(...) {
  .Deprecated("power_ppi_paired_binary", package = "pppower")
  power_ppi_paired_binary(...)
}

#' @rdname deprecated
#' @export
n_required_ppi_pp_paired <- function(delta,
                                      N,
                                      power = 0.80,
                                      alpha = 0.05,
                                      sigma_D2,
                                      rho_D,
                                      sigma_fD2 = NULL) {
  .Deprecated("power_ppi_paired", package = "pppower",
              msg = "n_required_ppi_pp_paired() is deprecated. Use power_ppi_paired(..., n = NULL, power = 0.80) instead.")
  power_ppi_paired(
    delta = delta,
    N = N,
    n = NULL,
    power = power,
    alpha = alpha,
    sigma_D2 = sigma_D2,
    rho_D = rho_D,
    sigma_fD2 = sigma_fD2
  )
}

#' @rdname deprecated
#' @export
power_ppi_pp_ttest <- function(...) {
  .Deprecated("power_ppi_ttest", package = "pppower")
  power_ppi_ttest(...)
}

#' @rdname deprecated
#' @export
power_ppi_pp_ttest_binary <- function(...) {
  .Deprecated("power_ppi_ttest_binary", package = "pppower")
  power_ppi_ttest_binary(...)
}

#' @rdname deprecated
#' @export
simulate_power_ppiplus_mean <- function(...) {
  .Deprecated("simulate_ppi_mean", package = "pppower")
  simulate_ppi_mean(...)
}

#' @rdname deprecated
#' @export
simulate_power_ppi_pp_ttest_binary <- function(...) {
  .Deprecated("simulate_ppi_ttest_binary", package = "pppower")
  simulate_ppi_ttest_binary(...)
}
