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
    # Allow either confusion-matrix pieces or precision/recall
    if (!is.null(stats$tp) && !is.null(stats$fp) && !is.null(stats$fn)) {
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
    } else {
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
#' @importFrom ranger ranger
#' @importFrom randomForest randomForest
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
        seed          = rf_seed,
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

#' Fit Predictive Model for OLS PPI / PPI++ (Internal Helper)
#'
#' @description
#' Internal unified interface for fitting predictive models used in the
#' OLS-based PPI and PPI++ estimators.  
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
#' @importFrom ranger ranger
#' @importFrom randomForest randomForest
fit_predict_model_ols <- function(model_type, X_L, y_L, X_U,
                              mtry = NULL,
                              rf_engine = c("ranger", "randomForest"),
                              rf_trees = 200,
                              rf_min_node_size = 5,
                              rf_num_threads = NULL,
                              rf_seed = NULL) {
  
  model_type <- match.arg(model_type, c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  rf_engine  <- match.arg(rf_engine, c("ranger", "randomForest"))
  
  ## Ensure column names for formula models
  if (is.matrix(X_L)) colnames(X_L) <- paste0("x", seq_len(ncol(X_L)))
  if (is.matrix(X_U)) colnames(X_U) <- paste0("x", seq_len(ncol(X_U)))
  
  dat_L <- data.frame(y = y_L, X_L)
  
  ## --- GLM MODELS ---
  if (model_type == "glm_correct") {
    
    ## Correct: true regression is y = β0 + β1 x1 + β2 x2
    #fit <- stats::glm(y ~ x1 + x2, data = dat_L)
    fit <- stats::glm(y ~ x1 + x2 + x3, data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))
    
  } else if (model_type == "glm_mis") {
    
    ## Misspecified: omit x2
    #fit <- stats::glm(y ~ x1, data = dat_L)
    fit <- stats::glm(y ~ x1 + x2, data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))
    
  } else if (model_type == "glm_wrong") {
    
    ## Wrong with interaction terms
    #fit <- stats::glm(y ~ I(x1 * x2), data = dat_L)
    fit <- stats::glm(y ~ x1, data = dat_L)
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))
    
  ## RANDOM FOREST 
  } else if (model_type == "rf") {
    
    mtry_eff <- if (is.null(mtry)) floor(sqrt(ncol(X_L))) else mtry
    use_rgr  <- (rf_engine == "ranger") && requireNamespace("ranger", quietly = TRUE)
    
    if (use_rgr) {
      fit <- ranger::ranger(
        y ~ ., 
        data = dat_L,
        num.trees     = rf_trees,
        mtry          = mtry_eff,
        min.node.size = rf_min_node_size,
        importance    = "none",
        seed          = rf_seed,
        num.threads   = rf_num_threads
      )
      fhat_U <- predict(fit, data = as.data.frame(X_U))$predictions
      
    } else {
      fit <- randomForest::randomForest(
        y ~ ., data = dat_L,
        ntree = rf_trees,
        mtry  = mtry_eff
      )
      fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U))
    }
  }
  
  list(
    fit    = fit,
    fhat_U = as.numeric(fhat_U)
  )
}

#' Fit Predictive Model for GLM PPI / PPI++ (Internal Helper)
#' @description
#' Internal unified interface for GLMs and random forests used in the
#' PPI / PPI++ GLM simulation framework.  
#' This function fits a predictive model on labeled data \code{(X_L, y_L)}
#' and returns fitted objects and predictions on both labeled and
#' unlabeled covariates.
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
#' @return A list with:
#'   \item{fit}{The fitted model object.}
#'   \item{fhat_L}{Predicted conditional means \eqn{\hat\mu_f(X_L)} on labeled data.}
#'   \item{fhat_U}{Predicted conditional means \eqn{\hat\mu_f(X_U)} on unlabeled data.}
#'
#' @importFrom stats glm predict binomial
#' @importFrom ranger ranger
#' @importFrom randomForest randomForest

fit_predict_model_glm <- function(
  model_type, X_L, y_L, X_U,
  mtry = NULL,
  rf_engine = c("ranger", "randomForest"),
  rf_trees = 200,
  rf_min_node_size = 5,
  rf_num_threads = NULL,
  rf_seed = NULL
) {

  model_type <- match.arg(model_type,
                          c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  rf_engine  <- match.arg(rf_engine, c("ranger", "randomForest"))

  ## Ensure column names for formula interface
  if (is.matrix(X_L)) colnames(X_L) <- paste0("x", seq_len(ncol(X_L)))
  if (is.matrix(X_U)) colnames(X_U) <- paste0("x", seq_len(ncol(X_U)))

  dat_L <- data.frame(y = y_L, X_L)

  ## GLM: CORRECT SPECIFICATION
  if (model_type == "glm_correct") {

    fit <- stats::glm(
      y ~ x1 + x2 + x3,
      data = dat_L,
      family = binomial()
    )

    fhat_L <- stats::predict(fit, newdata = dat_L, type = "response")
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U), type = "response")

  ## GLM: MIS-SPECIFIED (omits x2)
  } else if (model_type == "glm_mis") {

    fit <- stats::glm(
      y ~ x1 + x3,
      data = dat_L,
      family = binomial()
    )

    fhat_L <- stats::predict(fit, newdata = dat_L, type = "response")
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U), type = "response")

  ## GLM: WRONG MODEL (probit link + wrong structure)
  } else if (model_type == "glm_wrong") {

    fit <- stats::glm(
      y ~ x1 + I(x2^2),
      data = dat_L,
      family = binomial(link = "probit")
    )

    fhat_L <- stats::predict(fit, newdata = dat_L, type = "response")
    fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U), type = "response")

  ## RANDOM FOREST CLASSIFIER
  } else if (model_type == "rf") {

    ## Convert y to factor for classification
    dat_L_rf <- dat_L
    dat_L_rf$y <- factor(y_L, levels = c(0,1))

    mtry_eff <- if (is.null(mtry)) floor(sqrt(ncol(X_L))) else mtry

    if (rf_engine == "ranger" && requireNamespace("ranger", quietly = TRUE)) {

      fit <- ranger::ranger(
        y ~ ., data = dat_L_rf,
        probability = TRUE,
        num.trees   = rf_trees,
        mtry        = mtry_eff,
        min.node.size = rf_min_node_size,
        seed        = rf_seed,
        num.threads = rf_num_threads
      )

      fhat_L <- predict(fit, data = dat_L_rf)$predictions[, 2]
      fhat_U <- predict(fit, data = as.data.frame(X_U))$predictions[, 2]

    } else {

      fit <- randomForest::randomForest(
        y ~ ., data = dat_L_rf,
        ntree = rf_trees,
        mtry  = mtry_eff
      )

      fhat_L <- stats::predict(fit, newdata = dat_L_rf, type = "prob")[, 2]
      fhat_U <- stats::predict(fit, newdata = as.data.frame(X_U), type = "prob")[, 2]
    }
  }

  ## Output
  list(
    fit    = fit,
    fhat_L = as.numeric(fhat_L),
    fhat_U = as.numeric(fhat_U)
  )
}

simulate_one_draw_contrast <- function(
  n_labeled,
  X_sampler_L,
  f_generator,
  eps_sampler,
  a,
  delta
) {
  X <- X_sampler_L(n_labeled)
  X_mat <- as.matrix(X)

  f_true <- f_generator(X)
  f_delta <- f_true + delta * as.numeric(X_mat %*% a)

  y <- f_delta + eps_sampler(n_labeled)

  list(
    X = X,
    y = y
  )
}

compute_theta0 <- function(a, X_sampler_L, f_generator, n_ref = 5e5) {
  Xref <- X_sampler_L(n_ref)
  Xref_m <- as.matrix(Xref)
  yref <- f_generator(Xref_m)

  H <- crossprod(Xref_m) / n_ref
  G <- crossprod(Xref_m, yref) / n_ref
  beta_star <- solve(H, G)

  as.numeric(sum(a * beta_star))
}