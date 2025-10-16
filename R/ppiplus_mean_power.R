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
    
    if (is.null(m_labeled) && !is.null(metrics$m_obs)) m_labeled <- metrics$m_obs
    if (is.null(m_labeled)) stop("Need m_labeled (or metrics$m_obs) to recover covariances.")

    if (metric_type == "hard") {
      acc <- metrics$accuracy; p_y <- metrics$p_y; p_hat <- metrics$p_hat
      tp_rate <- (p_hat + p_y + acc - 1) / 2
      cov_y_f <- tp_rate - p_y * p_hat
    } 
    
    else if (metric_type == "precision_recall") {
      tp <- metrics$tp %||% (metrics$recall * metrics$p_y * m_labeled)
      p_y <- metrics$p_y
      p_hat <- metrics$p_hat %||% ((tp + metrics$fp)/m_labeled)
      cov_y_f <- (tp / m_labeled) - p_y * p_hat
    }
    
    else if (metric_type == "prob") {
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

simulate_power_ppiplus_mean <- function(df,
                                        N,
                                        n,
                                        theta0,
                                        alpha = 0.05,
                                        R = 2000,
                                        lambda = NULL,
                                        lambda_type = c("oracle", "plugin", "user"),
                                        use_sample_var = TRUE,
                                        seed = 1) {
  stopifnot(all(c("y", "fhat_cf") %in% names(df)))
  lambda_type <- match.arg(lambda_type)
  set.seed(seed)

  theta_true <- mean(df$y)
  delta <- theta_true - theta0

  sigma_y2 <- stats::var(df$y)
  sigma_f2 <- stats::var(df$fhat_cf)
  cov_yf   <- stats::cov(df$y, df$fhat_cf)

  lambda_star <- cov_yf / ((1 + n / N) * sigma_f2)
  var_oracle <- sigma_y2 / n - (cov_yf^2 / sigma_f2) * (N / (n * (n + N)))
  se_oracle <- sqrt(max(var_oracle, 0))
  z_alpha <- stats::qnorm(1 - alpha / 2)

  theta_hats <- numeric(R)
  se_hats    <- numeric(R)
  rej        <- logical(R)
  lambda_draw <- numeric(R)

  for (r in seq_len(R)) {
    draw <- ppplus_once(
      df = df,
      N = N,
      n = n,
      lambda = lambda,
      lambda_type = lambda_type,
      use_sample_var = use_sample_var,
      seed = NULL
    )
    theta_hats[r]  <- draw$theta_hat
    se_hats[r]     <- draw$se
    lambda_draw[r] <- draw$lambda
    rej[r] <- abs((draw$theta_hat - theta0) / draw$se) > z_alpha
  }

  list(
    theta_true        = theta_true,
    theta0            = theta0,
    delta             = delta,
    empirical_power   = mean(rej),
    analytical_power  = 1 - stats::pnorm(z_alpha - abs(delta) / se_oracle) +
                        stats::pnorm(-z_alpha - abs(delta) / se_oracle),
    avg_SE            = mean(se_hats),
    lambda_star       = lambda_star,
    lambda_draw       = lambda_draw,
    draws             = data.frame(theta_hat = theta_hats, se = se_hats, reject = rej)
  )
}
