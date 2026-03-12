# ============================================================================
# Power analysis for 2x2 contingency tables with PPI++ surrogates
# ============================================================================

#' Power or required sample size for a 2x2 table with PPI++ surrogates
#'
#' @description
#' Computes power or solves for the required labeled sample size for testing
#' an odds ratio (OR) or relative risk (RR) from a 2x2 contingency table
#' when a binary PPI++ surrogate is available.
#'
#' Follows the \pkg{pwr} convention: exactly one of `n` or `power` must be
#' `NULL`.
#'
#' **Odds ratio** (default): Uses logistic regression with design matrix
#' \eqn{[1, X]} where \eqn{X \in \{0,1\}} is the group indicator. The effect
#' is \eqn{\beta_1 = \log(\text{OR})}. Analytical variance blocks are
#' computed from the cell probabilities and surrogate performance, then
#' passed to [power_ppi_regression()].
#'
#' **Relative risk**: Uses the delta method on per-group PPI++ mean estimates.
#' The test statistic is a Wald test for \eqn{\log(\hat p_1 / \hat p_0)}.
#'
#' @param p_exp Numeric scalar. Probability of outcome in the exposed group,
#'   \eqn{P(Y = 1 \mid X = 1)}.
#' @param p_ctrl Numeric scalar. Probability of outcome in the control group,
#'   \eqn{P(Y = 1 \mid X = 0)}.
#' @param N Positive integer. Total number of unlabeled observations.
#' @param n Total number of labeled observations. Set to `NULL` to solve for
#'   the required `n`.
#' @param power Target power. Set to `NULL` to compute power from `n`.
#' @param alpha Significance level (default `0.05`).
#' @param prev_exp Numeric scalar. Prevalence of exposure, \eqn{P(X = 1)}.
#'   Default `0.5` (balanced groups).
#' @param sens Sensitivity of the binary surrogate classifier.
#' @param spec Specificity of the binary surrogate classifier.
#' @param effect_measure Character. `"OR"` (odds ratio, default) or `"RR"`
#'   (relative risk).
#' @param lambda_mode `"oracle"` (default), `"vanilla"`, or `"user"`.
#' @param lambda_user Scalar lambda when `lambda_mode = "user"`.
#' @param warn_smallN Logical. Warn when `N` is small.
#' @param ... Additional arguments (currently unused).
#'
#' @return When computing power: scalar in \[0, 1\].
#'   When computing sample size: integer `n` required.
#' @export
#'
#' @examples
#' # OR test: exposed group has higher risk
#' power_ppi_2x2(
#'   p_exp = 0.40, p_ctrl = 0.20,
#'   N = 500, n = 80, alpha = 0.05,
#'   sens = 0.85, spec = 0.85,
#'   effect_measure = "OR"
#' )
#'
#' # Solve for required n
#' power_ppi_2x2(
#'   p_exp = 0.40, p_ctrl = 0.20,
#'   N = 500, power = 0.80, alpha = 0.05,
#'   sens = 0.85, spec = 0.85,
#'   effect_measure = "OR"
#' )
#'
#' # RR test
#' power_ppi_2x2(
#'   p_exp = 0.40, p_ctrl = 0.20,
#'   N = 500, n = 80, alpha = 0.05,
#'   sens = 0.85, spec = 0.85,
#'   effect_measure = "RR"
#' )
power_ppi_2x2 <- function(
  p_exp,
  p_ctrl,
  N,
  n = NULL,
  power = NULL,

  alpha = 0.05,
  prev_exp = 0.5,
  sens,
  spec,
  effect_measure = c("OR", "RR"),
  lambda_mode = "oracle",
  lambda_user = NULL,
  warn_smallN = TRUE,
  ...
) {
  effect_measure <- match.arg(effect_measure)

  # Validate inputs
  if (is.null(n) && is.null(power))
    stop("Exactly one of 'n' and 'power' must be NULL.", call. = FALSE)
  if (!is.null(n) && !is.null(power))
    stop("Exactly one of 'n' and 'power' must be NULL.", call. = FALSE)

  for (v in c(p_exp, p_ctrl, prev_exp, sens, spec)) {
    if (!is.numeric(v) || length(v) != 1 || v <= 0 || v >= 1)
      stop("p_exp, p_ctrl, prev_exp, sens, spec must be scalars in (0,1).",
           call. = FALSE)
  }
  if (!is.numeric(N) || N <= 0) stop("N must be positive.", call. = FALSE)

  if (effect_measure == "OR") {
    .ppi_2x2_or(p_exp, p_ctrl, N, n, power, alpha, prev_exp, sens, spec,
                lambda_mode, lambda_user, warn_smallN)
  } else {
    .ppi_2x2_rr(p_exp, p_ctrl, N, n, power, alpha, prev_exp, sens, spec,
                lambda_mode, lambda_user)
  }
}


# --------------------------------------------------------------------------
# Internal: OR via logistic regression blocks
# --------------------------------------------------------------------------
.ppi_2x2_or <- function(p_exp, p_ctrl, N, n, power, alpha, prev_exp,
                         sens, spec, lambda_mode, lambda_user, warn_smallN) {
  blk <- .ppi_2x2_logistic_blocks(p_exp, p_ctrl, prev_exp, sens, spec)

  delta <- log(blk$OR)
  cvec  <- c(0, 1)  # contrast for beta_1

  power_ppi_regression(
    delta = delta,
    N = N,
    n = n,
    power = power,
    alpha = alpha,
    lambda_mode = lambda_mode,
    lambda_user = lambda_user,
    c = cvec,
    H_L = blk$H,
    H_U = blk$H,
    Sigma_YY   = blk$Sigma_YY,
    Sigma_ff_l = blk$Sigma_ff,
    Sigma_ff_u = blk$Sigma_ff,
    Sigma_Yf   = blk$Sigma_Yf,
    warn_smallN = warn_smallN
  )
}


# --------------------------------------------------------------------------
# Internal: Analytical logistic regression blocks for the 2x2 table
# --------------------------------------------------------------------------
#
# Model: logit P(Y=1|X) = beta_0 + beta_1 * X,  X in {0,1}
# Design matrix columns: (1, X).
#
# Population quantities for X ~ Bernoulli(pi), pi = prev_exp:
#   P(X=0) = 1-pi,  P(Y=1|X=0) = p_ctrl,  P(Y=1|X=1) = p_exp
#
# Surrogate: f|Y ~ Bernoulli(sens) if Y=1, Bernoulli(1-spec) if Y=0
#   => p_f_x = sens * p_x + (1-spec)*(1-p_x)  for x in {0,1}
#
.ppi_2x2_logistic_blocks <- function(p_exp, p_ctrl, prev_exp, sens, spec) {
  pi <- prev_exp
  p1 <- p_exp
  p0 <- p_ctrl

  # Logistic coefficients
  beta0 <- log(p0 / (1 - p0))
  beta1 <- log(p1 / (1 - p1)) - beta0
  OR <- exp(beta1)

  # Weights: w_x = p_x * (1-p_x)

  w0 <- p0 * (1 - p0)
  w1 <- p1 * (1 - p1)

  # Fisher information H = E[w(X) * x x^T]
  # x = (1, X)^T
  # E[w(X) * (1,X)^T (1,X)] = (1-pi)*w0*[[1,0],[0,0]] + pi*w1*[[1,1],[1,1]]
  H <- matrix(0, 2, 2)
  H[1, 1] <- (1 - pi) * w0 + pi * w1
  H[1, 2] <- pi * w1
  H[2, 1] <- pi * w1
  H[2, 2] <- pi * w1

  # Sigma_YY = H for correctly specified logistic
  Sigma_YY <- H

  # Surrogate prevalence per stratum
  p_f0 <- sens * p0 + (1 - spec) * (1 - p0)
  p_f1 <- sens * p1 + (1 - spec) * (1 - p1)

  # Sigma_ff: E[(f - mu(X))^2 * x x^T]
  # For binary f, Var(f - mu|X=x) = E[(f - p_x)^2 | X=x]
  #   = p_f_x * (1-p_x)^2 + (1-p_f_x) * p_x^2
  v_f0 <- p_f0 * (1 - p0)^2 + (1 - p_f0) * p0^2
  v_f1 <- p_f1 * (1 - p1)^2 + (1 - p_f1) * p1^2

  Sigma_ff <- matrix(0, 2, 2)
  Sigma_ff[1, 1] <- (1 - pi) * v_f0 + pi * v_f1
  Sigma_ff[1, 2] <- pi * v_f1
  Sigma_ff[2, 1] <- pi * v_f1
  Sigma_ff[2, 2] <- pi * v_f1

  # Sigma_Yf: E[(Y - mu(X))(f - mu(X)) * x x^T]
  # Cov(Y - p_x, f - p_x | X = x) = E[Y*f|X=x] - p_x * p_f_x
  #   Y*f = 1 iff Y=1 and f=1, P(Y=1,f=1|X=x) = p_x * sens
  #   So cov = p_x * sens - p_x * p_f_x = p_x * (sens - p_f_x)
  c_Yf0 <- p0 * (sens - p_f0)
  c_Yf1 <- p1 * (sens - p_f1)

  Sigma_Yf <- matrix(0, 2, 2)
  Sigma_Yf[1, 1] <- (1 - pi) * c_Yf0 + pi * c_Yf1
  Sigma_Yf[1, 2] <- pi * c_Yf1
  Sigma_Yf[2, 1] <- pi * c_Yf1
  Sigma_Yf[2, 2] <- pi * c_Yf1

  list(H = H, Sigma_YY = Sigma_YY, Sigma_ff = Sigma_ff, Sigma_Yf = Sigma_Yf,
       beta = c(beta0, beta1), OR = OR)
}


# --------------------------------------------------------------------------
# Internal: RR (relative risk) via delta method on per-group PPI++ means
# --------------------------------------------------------------------------
.ppi_2x2_rr <- function(p_exp, p_ctrl, N, n, power, alpha, prev_exp,
                         sens, spec, lambda_mode, lambda_user) {
  pi <- prev_exp
  p1 <- p_exp
  p0 <- p_ctrl
  RR <- p1 / p0
  delta <- log(RR)

  if (abs(delta) < 1e-15)
    stop("delta = log(RR) is zero; no effect to detect.", call. = FALSE)

  # Per-group sample sizes
  # N_1 = pi*N, N_0 = (1-pi)*N; similarly for n
  # Variance of log(RR_hat) ~ Var(p1_hat)/p1^2 + Var(p0_hat)/p0^2

  # Binary moments for each group
  mom1 <- binary_moments_from_sens_spec(p = p1, sens = sens, spec = spec)
  mom0 <- binary_moments_from_sens_spec(p = p0, sens = sens, spec = spec)

  # PPI++ oracle variance for each group's mean estimator
  # Var_ppi(p_hat_x) = sigma_y2_x/n_x - cov_yf_x^2/(sigma_f2_x * n_x * (1 + n_x/N_x))
  #                  = (sigma_y2_x - cov_yf_x^2/sigma_f2_x * N_x/(n_x+N_x)) / n_x
  # For oracle lambda mode.

  ppi_var_fun <- function(n_total) {
    n1 <- pi * n_total
    n0 <- (1 - pi) * n_total
    N1 <- pi * N
    N0 <- (1 - pi) * N

    if (lambda_mode == "oracle") {
      # Oracle variance for group x: sigma_y2/n_x - cov^2/(sf2 * n_x*(1+n_x/N_x))
      var_ppi_1 <- mom1$sigma_y2 / n1 -
        mom1$cov_y_f^2 / (mom1$sigma_f2 * n1 * (1 + n1 / N1))
      var_ppi_0 <- mom0$sigma_y2 / n0 -
        mom0$cov_y_f^2 / (mom0$sigma_f2 * n0 * (1 + n0 / N0))
    } else if (lambda_mode == "vanilla") {
      var_ppi_1 <- mom1$sigma_y2 / n1 + mom1$sigma_f2 / N1 +
        mom1$sigma_f2 / n1 - 2 * mom1$cov_y_f / n1
      var_ppi_0 <- mom0$sigma_y2 / n0 + mom0$sigma_f2 / N0 +
        mom0$sigma_f2 / n0 - 2 * mom0$cov_y_f / n0
    } else {
      lam <- lambda_user
      var_ppi_1 <- lam^2 * (mom1$sigma_f2 / N1 + mom1$sigma_f2 / n1) +
        mom1$sigma_y2 / n1 - 2 * lam * mom1$cov_y_f / n1
      var_ppi_0 <- lam^2 * (mom0$sigma_f2 / N0 + mom0$sigma_f2 / n0) +
        mom0$sigma_y2 / n0 - 2 * lam * mom0$cov_y_f / n0
    }

    # Delta method: Var(log(RR_hat)) ~ Var(p1_hat)/p1^2 + Var(p0_hat)/p0^2
    var_ppi_1 / p1^2 + var_ppi_0 / p0^2
  }

  z_alpha <- stats::qnorm(1 - alpha / 2)

  if (!is.null(n)) {
    # Compute power given n
    v <- ppi_var_fun(n)
    se <- sqrt(max(v, 0))
    if (se == 0) return(if (delta == 0) alpha else 1)
    mu <- abs(delta) / se
    return(1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu))
  }

  # Solve for n given power
  z_beta <- stats::qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2

  obj <- function(n_total) ppi_var_fun(n_total) - S2
  root <- tryCatch(
    stats::uniroot(obj, interval = c(4, 1e6))$root,
    error = function(e) NA_real_
  )

  if (is.na(root))
    stop("Infeasible: cannot find n achieving required power.", call. = FALSE)

  as.integer(ceiling(root))
}
