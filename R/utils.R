`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}

# Compute PP estimate + SE from a single draw (superpopulation)
pp_once <- function(df, N, n, use_sample_var = TRUE, seed = NULL) {
  # df must have columns: y, fhat_cf
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(c("y", "fhat_cf") %in% names(df)))

  idx <- sample.int(nrow(df), size = N + n, replace = FALSE)
  idx_u <- idx[1:N]
  idx_l <- idx[(N+1):(N+n)]

  y_u  <- df$y[idx_u]         # not used, unlabeled in PP estimator
  f_u  <- df$fhat_cf[idx_u]
  y_l  <- df$y[idx_l]
  f_l  <- df$fhat_cf[idx_l]

  # PP estimator
  A_N <- mean(f_u)
  B_n <- mean(y_l - f_l)
  theta_hat <- A_N + B_n

  # SE: superpopulation variance estimate
  if (use_sample_var) {
    # use sample variances from this draw
    var_f   <- stats::var(f_u)         # unbiased (ddof=1)
    var_res <- stats::var(y_l - f_l)
  } else {
    # use "population" variances from the whole df (plug-in)
    var_f   <- stats::var(df$fhat_cf)
    var_res <- stats::var(df$y - df$fhat_cf)
  }
  se <- sqrt(var_f / N + var_res / n)

  list(theta_hat = theta_hat, se = se,
       var_f = var_f, var_res = var_res,
       A_N = A_N, B_n = B_n)
}

ppplus_once <- function(df,
                        N,
                        n,
                        lambda = NULL,
                        lambda_type = c("oracle", "plugin", "user"),
                        use_sample_var = TRUE,
                        seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(c("y", "fhat_cf") %in% names(df)))
  lambda_type <- match.arg(lambda_type)

  idx <- sample.int(nrow(df), size = N + n, replace = FALSE)
  idx_u <- idx[seq_len(N)]
  idx_l <- idx[seq_len(n) + N]

  y_l <- df$y[idx_l]
  f_l <- df$fhat_cf[idx_l]
  f_u <- df$fhat_cf[idx_u]

  ybar_l <- mean(y_l)
  fbar_l <- mean(f_l)
  fbar_u <- mean(f_u)

  if (is.null(lambda)) {
    if (lambda_type == "oracle") {
      sigma_f2_pop <- stats::var(df$fhat_cf)
      cov_yf_pop   <- stats::cov(df$y, df$fhat_cf)
      lambda <- cov_yf_pop / ((1 + n / N) * sigma_f2_pop) # keep the same labmda
    } else if (lambda_type == "plugin") {
      sigma_f2_s <- stats::var(c(f_u, f_l))
      cov_yf_s   <- stats::cov(y_l, f_l)
      lambda <- cov_yf_s / ((1 + n / N) * sigma_f2_s) # substitute lambda within each draw
    } else {
      stop("Supply lambda when lambda_type = 'user'.", call. = FALSE)
    }
  }

  if (use_sample_var) {
    sigma_y2 <- stats::var(y_l)
    sigma_f2 <- stats::var(c(f_u, f_l))
    cov_yf   <- stats::cov(y_l, f_l)
  } else {
    sigma_y2 <- stats::var(df$y)
    sigma_f2 <- stats::var(df$fhat_cf)
    cov_yf   <- stats::cov(df$y, df$fhat_cf)
  }

  theta_hat <- ybar_l + lambda * (fbar_u - fbar_l)

  var_hat <- sigma_y2 / n +
    lambda^2 * sigma_f2 * (1 / N + 1 / n) -
    2 * lambda * cov_yf / n
  se <- sqrt(pmax(var_hat, 0))

  list(
    theta_hat = theta_hat,
    se = se,
    lambda = lambda,
    lambda_type = lambda_type,
    ybar_l = ybar_l,
    fbar_l = fbar_l,
    fbar_u = fbar_u,
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_yf = cov_yf
  )
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

derive_vars_binary <- function(metric_type = c("hard", "prob", "precision_recall"),
                               stats,
                               m_obs,
                               correction = TRUE) {
  metric_type <- match.arg(metric_type)

  if (!is.numeric(m_obs) || length(m_obs) != 1L || m_obs < 2) {
    stop("m_obs must be a scalar >= 2.")
  }
  adj <- if (isTRUE(correction)) m_obs / (m_obs - 1) else 1

  if (metric_type == "hard") {
    acc   <- stats$accuracy
    p_y   <- stats$p_y
    p_hat <- stats$p_hat

    if (any(!is.numeric(c(acc, p_y, p_hat))) ||
        any(c(acc, p_y, p_hat) < 0) ||
        any(c(acc, p_y, p_hat) > 1)) {
      stop("For metric_type = 'hard', accuracy, p_y, and p_hat must lie in [0, 1].")
    }

    err_rate <- 1 - acc
    bias     <- if (!is.null(stats$bias)) stats$bias else p_y - p_hat
    var_res  <- adj * (err_rate - bias^2)
    var_res  <- max(var_res, 0)  # guard against numerical negatives
    var_f    <- adj * p_hat * (1 - p_hat)
    var_f    <- max(var_f, 0)

    return(list(var_f = var_f,
                var_res = var_res,
                p_hat = p_hat,
                bias = bias))
  }

  if (metric_type == "precision_recall") {
    precision <- stats$precision
    recall    <- stats$recall
    p_y       <- stats$p_y

    if (any(!is.numeric(c(precision, recall, p_y))) ||
        precision <= 0 || precision > 1 ||
        recall < 0 || recall > 1 ||
        p_y < 0 || p_y > 1) {
      stop("For metric_type = 'precision_recall', precision in (0,1], recall in [0,1], p_y in [0,1].")
    }

    tp     <- recall * p_y * m_obs
    fp     <- tp * (1 / precision - 1)
    fn     <- p_y * m_obs - tp
    p_hat  <- (tp + fp) / m_obs
    err_rate <- (fp + fn) / m_obs
    bias     <- if (!is.null(stats$bias)) stats$bias else p_y - p_hat

    var_res <- adj * (err_rate - bias^2)
    var_res <- max(var_res, 0)
    var_f   <- adj * p_hat * (1 - p_hat)
    var_f   <- max(var_f, 0)

    return(list(var_f = var_f,
                var_res = var_res,
                p_hat = p_hat,
                bias = bias,
                tp = tp,
                fp = fp,
                fn = fn))
  }

  # metric_type == "prob"
  brier <- stats$brier
  if (!is.numeric(brier) || brier < 0) {
    stop("For metric_type = 'prob', provide a non-negative Brier score.")
  }

  bias <- if (!is.null(stats$bias)) stats$bias else 0
  var_res <- adj * (brier - bias^2)
  var_res <- max(var_res, 0)

  # sanity check
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

  if (!is.null(stats$p_hat)) {
    p_hat <- stats$p_hat
  } else {
    p_hat <- NA_real_  # mean predicted probability not identifiable from Brier alone
  }

  list(var_f = max(var_f, 0),
       var_res = var_res,
       var_y = var_y,
       p_hat = p_hat,
       bias = bias)
}

# brier score for binary case
brier_score <- function(y, p) {
  stopifnot(length(y) == length(p))
  mean((y - p)^2)
}

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

  if (metric_type %in% c("hard", "prob", "precision_recall")) {
    vars <- derive_vars_binary(
      metric_type = metric_type,
      stats = stats,
      m_obs = m_labeled,
      correction = correction
    )
    return(list(var_f = vars$var_f, var_res = vars$var_res))
  }

  if (metric_type %in% c("continuous", "regression")) {
    vars <- derive_vars_continuous(
      mse        = stats$mse,
      r2         = stats$r2,
      var_y      = stats$var_y,
      bias       = stats$bias %||% 0,
      m_obs      = m_labeled,
      correction = correction
    )
    return(list(var_f = vars$var_f, var_res = vars$var_res))
  }

  stop(sprintf("Unsupported metric_type: %s", metric_type), call. = FALSE)
}
