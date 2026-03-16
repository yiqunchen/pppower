#' Monte Carlo vs. analytical power for PPI / PPI++ mean estimation
#'
#' @description
#' `r lifecycle::badge("soft-deprecated")`
#'
#' Prefer [power_ppi_mean()] for PPI++ power calculations. This function is
#' retained for backward compatibility with `power_curve_mean()` and
#' `type1_error_curve_mean()`.
#'
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param alpha Two-sided significance level.
#' @param R Number of Monte Carlo draws (default 100000).
#' @param var_f Variance of \eqn{f(X)}.
#' @param var_res Variance of residuals \eqn{Y - f(X)}.
#' @param sigma_y2 Optional outcome variance (needed for PP++ if `cov_y_f` is supplied).
#' @param sigma_f2 Optional prediction variance (needed for PP++ if `cov_y_f` is supplied).
#' @param cov_y_f Optional covariance
#'   \eqn{\operatorname{Cov}(Y, f(X))}; when present, PP++ power is returned.
#' @param metrics Optional predictive summaries to recover missing moments.
#' @param metric_type Character string describing `metrics`.
#' @param m_labeled Labeled sample size associated with `metrics` (defaults to `n`).
#' @param correction Apply finite-sample corrections when recovering moments.
#'
#' @return Named numeric vector with `Exact_PP`, `Empirical_PP`,
#'   and when possible `Exact_PPplus`, `Empirical_PPplus`.
#' @export
simulate_power <- function(delta,
                           N,
                           n,
                           alpha = 0.05,
                           R = 100000,
                           var_f = NULL,
                           var_res = NULL,
                           sigma_y2 = NULL,
                           sigma_f2 = NULL,
                           cov_y_f = NULL,
                           metrics = NULL,
                           metric_type = NULL,
                           m_labeled = n,
                           correction = TRUE) {

  comps <- resolve_ppi_variances(
    var_f = var_f,
    var_res = var_res,
    metrics = metrics,
    metric_type = metric_type,
    m_labeled = m_labeled,
    correction = correction
  )
  var_f <- comps$var_f
  var_res <- comps$var_res

  se_pp <- sqrt(var_f / N + var_res / n)
  z_alpha <- stats::qnorm(1 - alpha / 2)
  mu_pp <- if (se_pp > 0) abs(delta) / se_pp else Inf

  if (se_pp == 0) {
    draws_pp <- rep(if (delta == 0) 0 else sign(delta) * Inf, R)
  } else {
    u_pp <- (seq_len(R) - 0.5) / R
    draws_pp <- stats::qnorm(u_pp, mean = mu_pp, sd = 1)
  }
  empirical_pp <- mean(abs(draws_pp) > z_alpha)
  exact_pp <- 1 - stats::pnorm(z_alpha - mu_pp) + stats::pnorm(-z_alpha - mu_pp)
  out <- c(Exact_PP = exact_pp, Empirical_PP = empirical_pp)

  cov_y_f_input <- cov_y_f %||% metrics$cov_y_f
  if (!is.null(cov_y_f_input)) {
    ppplus <- resolve_ppi_pp_moments(
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f_input,
      var_f = var_f,
      var_res = var_res,
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = m_labeled,
      correction = correction
    )

    var_ppplus <- ppi_pp_variance(
      n = n,
      N = N,
      sigma_y2 = ppplus$sigma_y2,
      sigma_f2 = ppplus$sigma_f2,
      cov_y_f = ppplus$cov_y_f,
      lambda_type = "oracle"
    )
    se_ppplus <- sqrt(var_ppplus)
    mu_ppplus <- if (se_ppplus > 0) abs(delta) / se_ppplus else Inf

    if (se_ppplus == 0) {
      draws_ppplus <- rep(if (delta == 0) 0 else sign(delta) * Inf, R)
    } else {
      u_ppplus <- (seq_len(R) - 0.5) / R
      draws_ppplus <- stats::qnorm(u_ppplus, mean = mu_ppplus, sd = 1)
    }
    empirical_ppplus <- mean(abs(draws_ppplus) > z_alpha)
    exact_ppplus <- 1 - stats::pnorm(z_alpha - mu_ppplus) + stats::pnorm(-z_alpha - mu_ppplus)
    out <- c(out, Exact_PPplus = exact_ppplus, Empirical_PPplus = empirical_ppplus)
  }

  out
}



#' Simulate cross-fitted predictive data
#'
#' @description
#' Generates a synthetic regression/classification dataset with covariates, outcomes,
#' and out-of-fold predictions from a K-fold cross-fitted GLM.
#'
#' @param n Number of observations to simulate.
#' @param p Number of covariates (columns in the design matrix).
#' @param family GLM family (defaults to `stats::binomial()`; `stats::gaussian()` gives a Gaussian response).
#' @param K Number of folds used for cross-fitting.
#' @param seed RNG seed for reproducibility.
#'
#' @return A data frame containing the outcome `y`, covariate columns `x1, ..., xp`, and
#' `fhat_cf`, the cross-fitted predictions. Attributes `family` and `K` record the generating family and fold count.
#' 
#' @export
simulate_crossfit_data <- function(n = 20000, p = 5,
                                   family = stats::binomial(),
                                   K = 5, seed = 1) {
  set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)

  beta <- seq(0.3, length.out = p, by = 0.5)

  if (identical(family$family, "binomial")) {
    # Logistic DGP
    eta <- drop(X %*% beta)
    prob <- stats::plogis(eta)
    y <- stats::rbinom(n, size = 1, prob = prob)
  } else {
    # Gaussian DGP
    mu <- drop(X %*% beta)
    y <- mu + stats::rnorm(n, sd = 2.0)
  }

  # Cross-fitted predictions (out-of-fold)
  fhat_cf <- crossfit_glm(X, y, K = K, family = family, seed = seed)

  # Return tidy data.frame
  df <- data.frame(y = y, X, fhat_cf = fhat_cf)
  attr(df, "family") <- family$family
  attr(df, "K") <- K
  df
}

draw_population_sample <- function(n, p, family = stats::gaussian(), seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X <- matrix(stats::rnorm(n * p), n, p)
  colnames(X) <- paste0("x", seq_len(p))
  beta <- seq(0.3, length.out = p, by = 0.5)
  if (identical(family$family, "binomial")) {
    eta <- drop(X %*% beta)
    prob <- stats::plogis(eta)
    y <- stats::rbinom(n, size = 1, prob = prob)
  } else {
    mu <- drop(X %*% beta)
    y <- mu + stats::rnorm(n, sd = 2.0)
  }
  list(X = X, y = y, family = family)
}

# Helper to draw labeled & unlabeled sample from the data generating distribution (OLS)
draw_labeled_unlabeled <- function(n_l, n_u, p, family = stats::gaussian(), seed = NULL) {
  labeled <- draw_population_sample(n_l, p, family = family, seed = seed)
  unlabeled <- draw_population_sample(n_u, p, family = family,
                                      seed = if (is.null(seed)) NULL else seed + 1)
  list(labeled = labeled, unlabeled = unlabeled)
}


crossfit_rectifier <- function(X_l, y_l, family, K = 2, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  folds <- kfold_split(nrow(X_l), K = K, seed = sample.int(.Machine$integer.max, 1))
  f_hat <- numeric(nrow(X_l))
  models <- vector("list", K)

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(nrow(X_l)), te)
    fit <- stats::glm(y ~ ., data = data.frame(y = y_l, X_l)[tr, , drop = FALSE],
                      family = family)
    f_hat[te] <- stats::predict(fit, newdata = data.frame(X_l)[te, , drop = FALSE],
                                type = if (identical(family$family, "binomial")) "response" else "response")
    models[[k]] <- fit
  }

  list(pred = f_hat, models = models, folds = folds)
}

crossfit_unlabeled <- function(X_u, pseudo_labels, family, folds, models_labeled) {
  # pseudo_labels should already be aligned with rows of X_u
  f_hat <- numeric(nrow(X_u))
  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(nrow(X_u)), te)
    fit <- stats::glm(y ~ ., data = data.frame(y = pseudo_labels, X_u)[tr, , drop = FALSE],
                      family = family)
    f_hat[te] <- stats::predict(fit, newdata = data.frame(X_u)[te, , drop = FALSE],
                                type = if (identical(family$family, "binomial")) "response" else "response")
  }
  f_hat
}

# Draw a population sample for mean estimation
draw_population_mean_sample <- function(n,
                                        theta_true,
                                        var_f,
                                        var_res = NULL,
                                        family = stats::gaussian(),
                                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  if (identical(family$family, "gaussian")) {
    if (is.null(var_res)) stop("Gaussian mean simulation requires var_res.")
    # Predictor component: f(X) centered at theta_true
    f <- stats::rnorm(n, mean = theta_true, sd = sqrt(var_f))
    # Residual noise
    eps <- stats::rnorm(n, mean = 0, sd = sqrt(var_res))
    # Observed outcome
    y <- f + eps

  } else if (identical(family$family, "binomial")) {
    if (theta_true <= 0 || theta_true >= 1) {
      stop("For binomial simulation, theta_true must lie in (0,1).")
    }
    if (var_f >= theta_true * (1 - theta_true)) {
      stop("var_f must be smaller than theta_true * (1 - theta_true).")
    }

    # Shape parameters of Beta prior for f(X)
    k_param <- theta_true * (1 - theta_true) / var_f - 1
    a_param <- theta_true * k_param
    b_param <- (1 - theta_true) * k_param

    f <- stats::rbeta(n, a_param, b_param)             # predictor distribution
    y <- stats::rbinom(n, size = 1, prob = f)          # binary outcome

    if (is.null(var_res)) {
      var_res <- (a_param * b_param) / ((a_param + b_param)^2 * (a_param + b_param + 1))
    }
  } else {
    stop("Unsupported family.")
  }

  list(
    y = y,
    f = f,
    family = family,
    theta_true = theta_true,
    var_f = var_f,
    var_res = var_res
  )
}

# Draw labeled & unlabeled samples for mean estimation
draw_labeled_unlabeled_mean <- function(n_l,
                                        n_u,
                                        theta0,
                                        delta,
                                        var_f,
                                        var_res = NULL,
                                        family = stats::gaussian(),
                                        seed = NULL) {
  theta_true <- theta0 + delta
  labeled <- draw_population_mean_sample(
    n = n_l,
    theta_true = theta_true,
    var_f = var_f,
    var_res = var_res,
    family = family,
    seed = seed
  )
  unlabeled <- draw_population_mean_sample(
    n = n_u,
    theta_true = theta_true,
    var_f = var_f,
    var_res = var_res,
    family = family,
    seed = if (is.null(seed)) NULL else seed + 1
  )
  list(labeled = labeled, unlabeled = unlabeled)
}


simulate_ppi_mean_rep <- function(n_l,
                                  n_u,
                                  family = stats::gaussian(),
                                  method = c("ppi", "ppi_plus"),
                                  lambda_type = c("plugin", "oracle", "user"),
                                  lambda_user = NULL,
                                  alpha = 0.05,
                                  theta0 = 0,
                                  delta = NULL,
                                  var_f = NULL,
                                  var_res = NULL,
                                  seed = NULL) {
  method <- match.arg(method)
  lambda_type <- match.arg(lambda_type)
  if (!is.null(seed)) set.seed(seed)

  # required population inputs for the DGP
  if (is.null(delta) || is.null(var_f)) {
    stop("Must supply `delta` and `var_f` for mean-based simulation.", call. = FALSE)
  }

  #  draw labeled & unlabeled from the MEAN DGP
  smp <- draw_labeled_unlabeled_mean(
    n_l = n_l, n_u = n_u,
    theta0 = theta0, delta = delta,
    var_f = var_f, var_res = var_res,
    family = family, seed = seed
  )

  y_l <- smp$labeled$y
  f_l <- smp$labeled$f
  y_u <- smp$unlabeled$y
  f_u <- smp$unlabeled$f
  theta_true <- smp$labeled$theta_true  # theta0 + delta

  # Vanilla PPI (mean): A_N + B_n
  A_N <- mean(f_u)
  B_n <- mean(y_l - f_l)
  theta_hat_pp <- A_N + B_n

  # variance components (natural estimators)
  # use unlabeled f for var_f; labeled residuals for var_res
  var_f_hat   <- stats::var(f_u)
  var_res_hat <- stats::var(y_l - f_l)
  se_pp <- sqrt(var_f_hat / n_u + var_res_hat / n_l)

  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_pp <- (theta_hat_pp - theta0) / se_pp
  reject_pp <- abs(z_pp) > z_alpha

  out <- list(
    method = method,
    theta_hat = theta_hat_pp,
    se = se_pp,
    reject = reject_pp,
    components = list(
      A_N = A_N,
      B_n = B_n,
      var_f_hat = var_f_hat,
      var_res_hat = var_res_hat,
      theta_true = theta_true
    )
  )

  # Mean estimation of PPI++
  if (method == "ppi_plus") {
    # empirical cov(Y,f) from labeled
    cov_yf_hat <- stats::cov(y_l, f_l, use = "complete.obs")
    # empirical var(f) from pooled f (stabilizes)
    sigma_f2_mix_hat <- stats::var(c(f_l, f_u), na.rm = TRUE)

    lambda <-
      if (lambda_type == "user") {
        if (is.null(lambda_user) || !is.finite(lambda_user)) {
          stop("Provide a finite scalar lambda_user when lambda_type='user'.", call. = FALSE)
        }
        lambda_user
      } else if (lambda_type == "plugin") {
        cov_yf_hat / ((1 + n_l / n_u) * sigma_f2_mix_hat)
      } else { # "oracle" (simulation-only; uses unlabeled Y_u which is unavailable in practice)
        cov_yf_oracle <- stats::cov(c(y_l, y_u), c(f_l, f_u), use = "complete.obs")
        sigma_f2_oracle <- stats::var(c(f_l, f_u), na.rm = TRUE)
        cov_yf_oracle / ((1 + n_l / n_u) * sigma_f2_oracle)
      }

    ybar_l <- mean(y_l)
    fbar_l <- mean(f_l)
    theta_hat_ppplus <- ybar_l + lambda * (A_N - fbar_l)

    # variance pieces for PPI++ (estimated from the draws)
    sigma_y2_hat <- stats::var(c(y_l, y_u))
    sigma_f2_hat <- stats::var(c(f_l, f_u))

    var_ppplus <- sigma_y2_hat / n_l +
      lambda^2 * sigma_f2_hat * (1 / n_u + 1 / n_l) -
      2 * lambda * cov_yf_hat / n_l
    var_ppplus <- max(var_ppplus, 0)
    se_ppplus <- sqrt(var_ppplus)

    z_ppplus <- (theta_hat_ppplus - theta0) / se_ppplus
    reject_ppplus <- abs(z_ppplus) > z_alpha

    out <- c(
      out,
      list(
        theta_hat_ppplus = theta_hat_ppplus,
        se_ppplus = se_ppplus,
        lambda = lambda,
        reject_ppplus = reject_ppplus,
        components_ppplus = list(
          sigma_y2_hat = sigma_y2_hat,
          sigma_f2_hat = sigma_f2_hat,
          cov_yf_hat = cov_yf_hat
        )
      )
    )
  }

  out
}
