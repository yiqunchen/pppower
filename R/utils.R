#' @importFrom stats lm sd pnorm cor
NULL

`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}


# stable least-squares via lm.fit
ols_fit <- function(X, y) {
  # X: n x p (already includes intercept column if desired)
  # y: n vector
  fit <- stats::lm.fit(X, y)
  # residuals and QR object come for free
  list(coef = fit$coefficients,
       resid = fit$residuals,
       qr    = fit$qr)
}

# OLS sandwich estimator
meat_matrix <- function(X, resid) {
  # returns M = E[(eps^2) x x^T] estimated by sample mean
  # M_hat = (1/n) sum_i (resid_i^2) x_i x_i^T
  n <- nrow(X)
  # vectorized: crossprod with residual weights
  WX <- X * as.numeric(resid)
  crossprod(WX) / n
}

bread_inv <- function(X) {
  # (X'X / n)^{-1}
  n <- nrow(X)
  solve(crossprod(X) / n, tol = 1e-12)
}

# k-fold splitting
kfold_split <- function(n, K = 5, seed = 1) {
  set.seed(seed)
  idx <- sample.int(n)
  folds <- split(idx, rep(1:K, length.out = n))
  lapply(folds, sort)
}

# helper to derive variances from mse and r2 (continuous case)
derive_vars_continuous <- function(mse,
                                   r2        = NULL,
                                   var_y     = NULL,
                                   bias      = 0,
                                   m_obs     = NULL,
                                   correction = TRUE) {
  stopifnot(is.numeric(mse), mse >= 0)
  if (is.null(m_obs)) stop("Need the sample size m_obs.")
  adj_factor <- if (correction) m_obs / (m_obs - 1) else 1
  var_res <- adj_factor * (mse - bias^2)
  if (!is.null(var_y)) {
    var_f <- var_y - var_res
  } else if (!is.null(r2)) {
    if (!(0 <= r2 && r2 <= 1)) stop("r2 must be in [0, 1].")
    var_y <- (m_obs * mse) / ((1 - r2) * (m_obs - 1))
    var_f <- r2 * var_y
  } else {
    stop("Provide either var_y or r2.")
  }
  list(var_f = var_f, var_res = var_res, var_y = var_y)
}

## Brier change to Y and Y_hat
derive_vars_binary <- function(metric_type = c("classification", "prob"),
                               stats,
                               m_obs,
                               correction = TRUE) {
  metric_type <- match.arg(metric_type)

  if (!is.numeric(m_obs) || length(m_obs) != 1L || m_obs < 2) {
    stop("m_obs must be a scalar >= 2.")
  }
  adj <- if (isTRUE(correction)) m_obs / (m_obs - 1) else 1

  if (metric_type == "classification") {
    # Allow: (1) confusion matrix, (2) sensitivity/specificity, or (3) precision/recall
    if (!is.null(stats$tp) && !is.null(stats$fp) && !is.null(stats$fn)) {
      # Path 1: Confusion matrix counts
      tp <- stats$tp
      fp <- stats$fp
      fn <- stats$fn
      tn <- stats$tn %||% (m_obs - tp - fp - fn)
      if (any(c(tp, fp, fn, tn) < 0) || abs(tp + fp + fn + tn - m_obs) > 1e-8) {
        stop("Confusion-matrix counts must be non-negative and sum to m_obs.")
      }
      p_y   <- (tp + fn) / m_obs
      p_hat <- (tp + fp) / m_obs
      err_rate <- (fp + fn) / m_obs

    } else if (!is.null(stats$sensitivity) && !is.null(stats$specificity)) {
      # Path 2: Sensitivity/specificity + prevalence
      sensitivity <- stats$sensitivity
      specificity <- stats$specificity
      p_y <- stats$p_y
      if (any(!is.numeric(c(sensitivity, specificity, p_y))) ||
          sensitivity < 0 || sensitivity > 1 ||
          specificity < 0 || specificity > 1 ||
          p_y < 0 || p_y > 1) {
        stop("Need sensitivity in [0,1], specificity in [0,1], and p_y (prevalence) in [0,1].")
      }
      # Derive confusion matrix from sensitivity/specificity
      tp <- sensitivity * p_y * m_obs
      fn <- (1 - sensitivity) * p_y * m_obs
      tn <- specificity * (1 - p_y) * m_obs
      fp <- (1 - specificity) * (1 - p_y) * m_obs
      p_hat <- (tp + fp) / m_obs
      err_rate <- (fp + fn) / m_obs

    } else {
      # Path 3: Precision/recall + prevalence
      precision <- stats$precision
      recall    <- stats$recall
      p_y       <- stats$p_y
      if (any(!is.numeric(c(precision, recall, p_y))) ||
          precision <= 0 || precision > 1 ||
          recall < 0 || recall > 1 ||
          p_y < 0 || p_y > 1) {
        stop("Need precision in (0,1], recall in [0,1], and p_y in [0,1] when supplying precision/recall.")
      }
      tp <- recall * p_y * m_obs
      fp <- tp * (1 / precision - 1)
      fn <- p_y * m_obs - tp
      if (min(tp, fp, fn) < 0) stop("Inconsistent precision/recall/p_y triple.")
      p_hat   <- (tp + fp) / m_obs
      err_rate <- (fp + fn) / m_obs
    }

    bias <- if (!is.null(stats$bias)) stats$bias else p_y - p_hat
    var_res <- adj * (err_rate - bias^2)
    var_res <- max(var_res, 0)
    var_f   <- adj * p_hat * (1 - p_hat)
    var_f   <- max(var_f, 0)

    return(list(
      var_f = var_f,
      var_res = var_res,
      p_hat = p_hat,
      bias = bias,
      tp = tp,
      fp = fp,
      fn = fn
    ))
  }

  # metric_type == "prob"
  brier <- stats$brier
  if (!is.numeric(brier) || brier < 0) {
    stop("For metric_type = 'prob', provide a non-negative Brier score.")
  }

  bias <- if (!is.null(stats$bias)) stats$bias else 0
  var_res <- adj * (brier - bias^2)
  var_res <- max(var_res, 0)

  if (!is.null(stats$var_y)) {
    var_y <- stats$var_y
    if (!is.numeric(var_y) || var_y < var_res) {
      stop("var_y must be numeric and at least as large as var_res.")
    }
    var_f <- var_y - var_res
  } else if (!is.null(stats$r2)) {
    r2 <- stats$r2
    if (!is.numeric(r2) || r2 < 0 || r2 >= 1) {
      stop("r2 must lie in [0, 1) when supplied.")
    }
    var_f <- (r2 / (1 - r2)) * var_res
    var_y <- var_f + var_res
  } else if (!is.null(stats$p_y)) {
    p_y <- stats$p_y
    if (!is.numeric(p_y) || p_y < 0 || p_y > 1) {
      stop("p_y must lie in [0, 1].")
    }
    var_y <- p_y * (1 - p_y)
    if (var_res > var_y) {
      stop("var_res exceeds var_y implied by p_y; check inputs.")
    }
    var_f <- var_y - var_res
  } else {
    stop("Supply at least one of var_y, r2, or p_y for metric_type = 'prob'.")
  }

  p_hat <- stats$p_hat %||% NA_real_

  list(
    var_f = max(var_f, 0),
    var_res = var_res,
    var_y = var_y,
    p_hat = p_hat,
    bias = bias
  )
}

# brier score for binary case
brier_score <- function(y, p) {
  stopifnot(length(y) == length(p))
  mean((y - p)^2)
}

#' Moments for binary outcome/prediction pairs from sensitivity/specificity
#'
#' @param p Prevalence \eqn{P(Y = 1)}.
#' @param sens Classifier sensitivity \eqn{P(\hat Y = 1 \mid Y = 1)}.
#' @param spec Classifier specificity \eqn{P(\hat Y = 0 \mid Y = 0)}.
#'
#' @return List with `sigma_y2`, `sigma_f2`, `cov_y_f`, and `p_hat`.
#' @export
binary_moments_from_sens_spec <- function(p, sens, spec) {
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
    stop("p must be a scalar prevalence in (0, 1).", call. = FALSE)
  }
  if (any(!is.finite(c(sens, spec))) || any(c(sens, spec) < 0) || any(c(sens, spec) > 1)) {
    stop("sens and spec must lie in [0, 1].", call. = FALSE)
  }

  p_hat <- sens * p + (1 - spec) * (1 - p)

  sigma_y2 <- p * (1 - p)
  sigma_f2 <- p_hat * (1 - p_hat)

  # Cov(Y, f) = P(Y=1, f=1) - P(Y=1)P(f=1) = p(1-p)(sens + spec - 1)
  cov_y_f <- p * (1 - p) * (sens + spec - 1)

  list(
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    p_hat = p_hat
  )
}

#' Resolve PP variance components from metrics
#'
#' @param var_f Optional numeric scalar supplying \eqn{\operatorname{Var}(f)} directly.
#' @param var_res Optional numeric scalar supplying \eqn{\operatorname{Var}(Y-f)} directly.
#' @param metrics Optional named list of predictive-performance summaries. The required
#'   fields depend on `metric_type`.
#' @param metric_type Character string identifying the metric bundle. Supported values are:
#'   `"continuous"` (regression-style metrics), `"prob"` (binary probabilistic metrics such
#'   as the Brier score), and `"classification"` (binary classification metrics such as a
#'   confusion matrix or precision/recall).
#' @param m_labeled Labeled sample size associated with `metrics`; defaults to
#'   `metrics$m_obs` when present.
#' @param correction Logical; apply the finite-sample adjustment when deriving moments
#'   from metrics (default `TRUE`).
#'
#' @return List with numeric elements `var_f` and `var_res`.
#' @keywords internal
resolve_ppi_variances <- function(var_f = NULL,
                                  var_res = NULL,
                                  metrics = NULL,
                                  metric_type = NULL,
                                  m_labeled = NULL,
                                  correction = TRUE) {
  if (!is.null(var_f) && !is.null(var_res)) {
    return(list(var_f = as.numeric(var_f), var_res = as.numeric(var_res)))
  }

  if (is.null(metrics)) {
    stop("Supply either (var_f, var_res) or a metrics list.", call. = FALSE)
  }

  # allow the metrics list itself to carry the type or sample size
  if (is.null(metric_type) && !is.null(metrics$type)) {
    metric_type <- metrics$type
  }
  if (is.null(m_labeled) && !is.null(metrics$m_obs)) {
    m_labeled <- metrics$m_obs
  }
  if (is.null(metric_type)) {
    stop("metric_type must be provided (or stored in metrics$type).", call. = FALSE)
  }
  if (is.null(m_labeled)) {
    stop("m_labeled must be provided (or stored in metrics$m_obs).", call. = FALSE)
  }

  stats <- metrics
  stats$type <- NULL
  stats$m_obs <- NULL

  metric_type_clean <- match.arg(tolower(metric_type),
                                 c("classification", "prob", "continuous"))

  if (metric_type_clean == "continuous") {
    vars <- derive_vars_continuous(
      mse        = stats$mse,
      r2         = stats$r2,
      var_y      = stats$var_y,
      bias       = stats$bias %||% 0,
      m_obs      = m_labeled,
      correction = correction
    )
  } else {
    vars <- derive_vars_binary(
      metric_type = metric_type_clean,
      stats = stats,
      m_obs = m_labeled,
      correction = correction
    )
  }

  list(var_f = vars$var_f, var_res = vars$var_res)
}

compute_hessian_fisher <- function(
  model_type = c("ols", "glm"),
  X,
  Y = NULL,         # optional for OLS, unused for GLM fisher
  beta = NULL,      # required for GLM
  family = NULL     # required for GLM
) {
  model_type <- match.arg(model_type)
  n <- nrow(X)

  if (model_type == "ols") {

    # H = (1/n) Xᵀ X
    H <- crossprod(X) / n
    return(H)
  }

  # GLM case (canonical link)
  if (is.null(beta))
    stop("GLM requires 'beta' (coefficient vector).")

  if (is.null(family))
    stop("GLM requires 'family' (e.g., 'binomial','poisson','gaussian').")

  eta <- as.numeric(X %*% beta)

  W <- switch(family,
              "binomial" = {
                mu <- 1 / (1 + exp(-eta))
                mu * (1 - mu)   # canonical logit variance
              },
              "poisson" = {
                mu <- exp(eta)
                mu               # canonical log variance
              },
              "gaussian" = rep(1, n),   # canonical identity
              stop("Unsupported GLM family: ", family)
  )

  Wmat <- diag(W, n, n)

  J <- t(X) %*% Wmat %*% X / n
  return(J)
}

compute_sigma_blocks <- function(
  X_l, Y_l, f_l,
  X_u, f_u,
  model_type = c("ols","glm"),
  beta = NULL,
  family = NULL
) {
  model_type <- match.arg(model_type)

  n_l <- nrow(X_l)
  n_u <- nrow(X_u)

  # Compute means μ_l, μ_u for GLM or OLS fitted value
  if (model_type == "ols") {

    mu_l <- as.numeric(X_l %*% beta)
    mu_u <- as.numeric(X_u %*% beta)

  } else {

    eta_l <- as.numeric(X_l %*% beta)
    eta_u <- as.numeric(X_u %*% beta)

    mu_l <- switch(family,
                   "binomial" = 1 / (1 + exp(-eta_l)),
                   "poisson"  = exp(eta_l),
                   "gaussian" = eta_l,
                   stop("Unsupported GLM family: ", family)
    )

    mu_u <- switch(family,
                   "binomial" = 1 / (1 + exp(-eta_u)),
                   "poisson"  = exp(eta_u),
                   "gaussian" = eta_u,
                   stop("Unsupported GLM family: ", family)
    )
  }

  # Score residuals
  res_Y_l <- Y_l - mu_l
  res_f_l <- f_l - mu_l
  res_f_u <- f_u - mu_u

  # Σ blocks
  Sigma_YY   <- crossprod(X_l * res_Y_l) / n_l
  Sigma_ff_l <- crossprod(X_l * res_f_l) / n_l
  Sigma_ff_u <- crossprod(X_u * res_f_u) / n_u
  Sigma_Yf   <- crossprod(X_l * res_Y_l, X_l * res_f_l) / n_l

  list(
    Sigma_YY   = Sigma_YY,
    Sigma_ff_l = Sigma_ff_l,
    Sigma_ff_u = Sigma_ff_u,
    Sigma_Yf   = Sigma_Yf
  )
}

#' Construct Hessian/Fisher and Covariance Blocks for PPI/PPI++ Regression
#'
#' @description
#' Computes all required Hessian/Fisher information matrices and covariance
#' blocks used by Prediction-Powered Inference (PPI) and PPI++ for
#' regression-based estimands (OLS or GLM contrasts).
#'
#' This helper wraps:
#'   * \code{compute_hessian_fisher} — model Hessian or Fisher information  
#'   * \code{compute_sigma_blocks}   --- Sigma-block covariance components
#'
#' and produces a unified object that can be passed directly to
#' [power_ppi_regression()] for labeled sample size calculations.
#'
#' @param model_type Character string: `"ols"` or `"glm"`.
#'   Determines whether the estimator is linear regression (OLS) or a GLM.
#'
#' @param X_l Matrix of labeled covariates (n x p).
#' @param Y_l Vector of labeled responses (required for OLS and GLM).
#' @param f_l Vector of model predictions on labeled data.
#'
#' @param X_u Matrix of unlabeled covariates (N x p).
#' @param f_u Vector of model predictions on unlabeled data.
#'
#' @param beta Numeric vector of regression coefficients used for Hessian/Fisher
#'   evaluation. Required for GLM; optional for OLS.
#'
#' @param family GLM family (`"binomial"` or `"gaussian"`). Only used for
#'   `model_type = "glm"`.
#'
#' @return
#' A named list of matrices:
#'
#' \describe{
#'   \item{H_L}{Labeled-data Hessian / Fisher information (p x p)}
#'   \item{H_U}{Unlabeled-data Hessian / Fisher information (p x p)}
#'   \item{Sigma_YY}{Covariance of labeled score (`Y - X beta`) (p x p)}
#'   \item{Sigma_ff_l}{Covariance of prediction score on labeled data (p x p)}
#'   \item{Sigma_ff_u}{Covariance of prediction score on unlabeled data (p x p)}
#'   \item{Sigma_Yf}{Cross-covariance between labeled scores and prediction scores (p x p)}
#' }
#'
#' These matrices are exactly the inputs needed for
#' [power_ppi_regression()] in `type = "regression"` mode.
#'
#' @examples
#' set.seed(1)
#' p <- 3
#' n_l <- 400
#' n_u <- 2000
#'
#' X_l <- matrix(rnorm(n_l * p), n_l, p)
#' X_u <- matrix(rnorm(n_u * p), n_u, p)
#' beta <- c(1, -0.5, 0.3)
#'
#' # Labeled response and predictions
#' Y_l <- drop(X_l %*% beta + rnorm(n_l))
#' f_l <- drop(X_l %*% beta + rnorm(n_l, sd = 0.3))
#' f_u <- drop(X_u %*% beta + rnorm(n_u, sd = 0.3))
#'
#' blocks <- compute_ppi_blocks(
#'   model_type = "ols",
#'   X_l = X_l, Y_l = Y_l, f_l = f_l,
#'   X_u = X_u, f_u = f_u,
#'   beta = beta
#' )
#'
#' # Use in sample size solver:
#' c_vec <- c(1, 0, 0)
#' delta <- as.numeric(t(c_vec) %*% beta)
#'
#' power_ppi_regression(
#'   delta = delta, N = n_u, power = 0.80,
#'   lambda_mode = "oracle",
#'   c = c_vec,
#'   H_L = blocks$H_L, H_U = blocks$H_U,
#'   Sigma_YY = blocks$Sigma_YY,
#'   Sigma_ff_l = blocks$Sigma_ff_l,
#'   Sigma_ff_u = blocks$Sigma_ff_u,
#'   Sigma_Yf = blocks$Sigma_Yf
#' )
#'
#' @seealso \link{power_ppi_regression}
#'
#' @export
compute_ppi_blocks <- function(
  model_type = c("ols","glm"),
  X_l, Y_l, f_l,
  X_u, f_u,
  beta = NULL,
  family = NULL
) {
  model_type <- match.arg(model_type)

  # Hessians / Fisher information
  H_L <- compute_hessian_fisher(
    model_type = model_type,
    X = X_l,
    Y = Y_l,
    beta = beta,
    family = family
  )

  H_U <- compute_hessian_fisher(
    model_type = model_type,
    X = X_u,
    Y = NULL,
    beta = beta,
    family = family
  )

  # Covariance Σ-blocks
  sig <- compute_sigma_blocks(
    X_l = X_l, Y_l = Y_l, f_l = f_l,
    X_u = X_u, f_u = f_u,
    model_type = model_type,
    beta = beta,
    family = family
  )

  # Output list
  list(
    H_L         = H_L,
    H_U         = H_U,
    Sigma_YY    = sig$Sigma_YY,
    Sigma_ff_l  = sig$Sigma_ff_l,
    Sigma_ff_u  = sig$Sigma_ff_u,
    Sigma_Yf    = sig$Sigma_Yf
  )
}

#' Generate One Labeled + Unlabeled Sample
#'
#' @keywords internal
#' 
#' @param n_labeled Number of labeled samples
#' @param n_unlabeled Number of unlabeled samples
#' @param X_sampler_L Function sampling labeled covariates
#' @param X_sampler_U Function sampling unlabeled covariates (defaults to X_sampler_L)
#' @param f_generator Function generating f(X)
#' @param eps_sampler Error sampler, default N(0,1)
#' @param delta Mean shift applied to f
#'
#' @return A list with labeled and unlabeled data and true mean.
#'
#' @importFrom stats rnorm
simulate_one_draw <- function(n_labeled, 
                              n_unlabeled,
                              X_sampler_L, 
                              X_sampler_U = NULL,
                              f_generator,
                              eps_sampler = function(n) rnorm(n, 0, 1),
                              delta = 0) {
  if (is.null(X_sampler_U)) X_sampler_U <- X_sampler_L
 
  X_L <- X_sampler_L(n_labeled)
  X_U <- X_sampler_U(n_unlabeled)

  f_L_base <- f_generator(X_L)
  f_U_base <- f_generator(X_U)

  # mean shift by delta
  f_L <- f_L_base + delta
  f_U <- f_U_base + delta

  y_L <- f_L + eps_sampler(n_labeled)

  list(
    labeled   = list(X = X_L, y = y_L, f = f_L),
    unlabeled = list(X = X_U, f = f_U),
    theta_true = mean(f_L),   # population mean of Y is ~ mean(f) here
    delta = delta
  )
}

#' Internal helper for model prediction
#'
#' @keywords internal
#'
#' @param fit Fitted model object
#' @param newdata New data frame to predict on
#'
#' @return Numeric vector of predictions.
#'
#' @importFrom stats predict
.predict_any <- function(fit, newdata) {
  if (inherits(fit, "ranger")) {
    predict(fit, data = as.data.frame(newdata))$predictions
  } else {
    stats::predict(fit, newdata = as.data.frame(newdata))
  }
}

#' Fit predictive model and return predictions for Mean Estimation for PPI / PPI++ 
#'
#' @description
#' Internal unified interface for fitting predictive models used in the
#' PPI and PPI++ estimators.  
#' Supports correctly specified, misspecified, and incorrectly specified
#' linear models, as well as random forests.  
#' Returns fitted model object and predictions on the unlabeled covariates.
#'
#' @keywords internal
#'
#' @param model_type Type of model ("glm_correct", "glm_mis", "glm_wrong", "rf")
#' @param X_L Labeled covariates
#' @param y_L Labeled outcomes
#' @param X_U Unlabeled covariates
#' @param mtry RF mtry value
#' @param rf_engine Random forest engine ("ranger" or "randomForest")
#' @param rf_trees Number of trees
#' @param rf_min_node_size Minimum node size for ranger
#' @param rf_num_threads Threads (ranger)
#' @param rf_seed Random seed
#'
#' @return List containing model fit object and predictions on X_U.
#'
#' @importFrom stats glm predict
fit_predict_model <- function(model_type, X_L, y_L, X_U,
                              mtry = NULL,
                              rf_engine = c("ranger", "randomForest"),
                              rf_trees = 200,
                              rf_min_node_size = 5,
                              rf_num_threads = NULL,
                              rf_seed = NULL) {
  model_type <- match.arg(model_type, c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  rf_engine  <- match.arg(rf_engine, c("ranger", "randomForest"))

  ## ensure names for formula use
  if (is.matrix(X_L)) colnames(X_L) <- paste0("x", seq_len(ncol(X_L)))
  if (is.matrix(X_U)) colnames(X_U) <- paste0("x", seq_len(ncol(X_U)))

  dat_L <- data.frame(y = y_L, X_L)

  if (model_type == "glm_correct") {
    fit <- stats::glm(y ~ x1 + I(x1 * x2) + I(x2^3), data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))

  } else if (model_type == "glm_mis") {
    fit <- stats::glm(y ~ x1 + x2, data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))

  } else if (model_type == "glm_wrong") {
    fit <- stats::glm(y ~ I(x1 * x2), data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))

  } else if (model_type == "rf") {
    mtry_eff  <- if (is.null(mtry)) floor(sqrt(ncol(X_L))) else mtry
    use_rgr   <- (rf_engine == "ranger") && requireNamespace("ranger", quietly = TRUE)

    if (use_rgr) {
      fit <- ranger::ranger(
        y ~ ., data = dat_L,
        num.trees     = rf_trees,
        mtry          = mtry_eff,
        min.node.size = rf_min_node_size,
        importance    = "none",
        num.threads   = rf_num_threads
      )
      fhat_U <- predict(fit, data = as.data.frame(X_U))$predictions

    } else {
      fit <- randomForest::randomForest(y ~ ., data = dat_L, ntree = rf_trees, mtry = mtry_eff)
      fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))
    }
  }

  list(fit = fit, fhat_U = as.numeric(fhat_U))
}
