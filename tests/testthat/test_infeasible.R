test_that("n_required_pp correctly detects infeasible configurations", {
  # Case where unlabeled term dominates: var_f / N >= delta^2 / c0^2, also for small N
  delta  <- 0.1
  var_f  <- 2.0
  var_res <- 1.0
  N <- 100
  alpha <- 0.05
  power <- 0.8

  expect_error(
    n_required_pp(
      delta = delta,
      N = N,
      alpha = alpha,
      power = power,
      type = "mean",
      var_f = var_f,
      var_res = var_res
    ),
    regexp = "Infeasible"
  )

  delta2 <- 0.5
  expect_warning(
    n_required_pp(
      delta = delta2,
      N = N,
      alpha = alpha,
      power = power,
      type = "mean",
      var_f = var_f,
      var_res = var_res,
      mode = "cap"
    ),
    regexp = "small|Unlabeled N is too small"
  )

  # Feasible configuration should not throw
  delta  <- 0.5
  var_f  <- 0.1
  var_res <- 1.0
  N <- 1000

  expect_no_error(
    n_required_pp(
      delta = delta,
      N = N,
      alpha = alpha,
      power = power,
      type = "mean",
      var_f = var_f,
      var_res = var_res
    )
  )

  N_big <- 1e6
  expect_no_error(
    n_required_pp(
      delta = delta,
      N = N_big,
      alpha = alpha,
      power = power,
      type = "mean",
      var_f = var_f,
      var_res = var_res
    )
  )
})


test_that("n_required_pp returns expected sample size for the mean estimator", {
  delta <- 0.2
  N <- 1000
  alpha <- 0.05
  power <- 0.8
  var_f <- 0.3
  var_res <- 1.1

  c0 <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  denom <- (delta^2 / c0^2) - (var_f / N)
  expected_n <- as.integer(ceiling(var_res / denom))

  result <- n_required_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    type = "mean",
    var_f = var_f,
    var_res = var_res,
    warn_smallN = FALSE
  )

  expect_identical(result, expected_n)
})

test_that("n_required_pp matches analytical solution for PP-OLS contrasts", {
  delta <- 0.15
  N <- 4000
  alpha <- 0.05
  power <- 0.9
  V_u <- diag(c(0.8, 0.6))
  V_l <- diag(c(0.5, 0.7))
  c_vec <- c(1, -0.5)

  var_u <- as.numeric(drop(t(c_vec) %*% V_u %*% c_vec))
  var_l <- as.numeric(drop(t(c_vec) %*% V_l %*% c_vec))
  c0 <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  denom <- (delta^2 / c0^2) - (var_u / N)
  expected_n <- as.integer(ceiling(var_l / denom))

  result <- n_required_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    type = "ols",
    V_u = V_u,
    V_l = V_l,
    c = c_vec,
    warn_smallN = FALSE
  )

  expect_identical(result, expected_n)
})

test_that("n_required_pp caps infeasible requests and records achieved power", {
  delta <- 0.05
  N <- 120
  alpha <- 0.05
  power <- 0.9
  var_unlabeled <- 0.01
  var_labeled <- 2

  expect_warning(
    capped <- n_required_pp(
      delta = delta,
      N = N,
      alpha = alpha,
      power = power,
      type = "custom",
      var_unlabeled = var_unlabeled,
      var_labeled = var_labeled,
      mode = "cap",
      warn_smallN = FALSE
    ),
    regexp = "Capping to n = N"
  )

  expect_identical(as.integer(capped), as.integer(N))
  achieved <- attr(capped, "achieved_power")
  expect_true(is.numeric(achieved) && length(achieved) == 1L)
  expect_lt(achieved, power)
  expect_gt(achieved, 0)
})

test_that("n_required_pp validates required inputs by type", {
  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_res = 1
    ),
    regexp = "var_f"
  )

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_f = -0.1,
      var_res = 1
    ),
    regexp = "var_f must be >= 0"
  )

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_f = 0.1,
      var_res = 0
    ),
    regexp = "var_res must be > 0"
  )

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "ols",
      V_u = matrix(1),
      c = c(1, 0)
    ),
    regexp = "V_u, V_l, and c must be supplied"
  )

  V_u <- matrix(0.1, nrow = 2, ncol = 2)
  V_l <- diag(c(0.5, 0.7))

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "ols",
      V_u = V_u,
      V_l = V_l,
      c = c(1, NA)
    ),
    regexp = "Contrast c must be a non-empty numeric vector without NA"
  )

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "custom",
      var_unlabeled = -1,
      var_labeled = 1
    ),
    regexp = "var_unlabeled must be non-negative"
  )

  expect_error(
    n_required_pp(
      delta = 0.2,
      N = 600,
      type = "custom",
      var_unlabeled = 1,
      var_labeled = 0
    ),
    regexp = "var_labeled must be positive"
  )
})


test_that("mean-mode: oracle lambda reproduces quadratic closed-form", {
  delta <- 0.25
  N <- 3000
  alpha <- 0.05
  power <- 0.8
  sigma_y2 <- 1.1
  sigma_f2 <- 0.4
  cov_y_f  <- 0.2
  var_f <- sigma_f2
  var_res <- sigma_y2 - sigma_f2

  z_alpha <- stats::qnorm(1 - alpha / 2)
  z_beta  <- stats::qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2

  discrim <- (sigma_y2 - S2 * N)^2 +
    4 * S2 * N * (sigma_y2 - (cov_y_f^2 / sigma_f2))

  expected_n <- ceiling((sigma_y2 - S2 * N + sqrt(discrim)) / (2 * S2))

  result <- n_required_ppi_pp(
    delta = delta,
    N     = N,
    alpha = alpha,
    power = power,
    type  = "mean",
    lambda_mode = "oracle",
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f,
    var_f    = var_f,
    var_res  = var_res,
    warn_smallN = FALSE
  )

  expect_identical(result, as.integer(expected_n))
})


test_that("mean-mode: metrics input matches direct moment input", {
  delta <- 0.18
  N <- 2500
  alpha <- 0.05
  power <- 0.85
  sigma_y2 <- 0.9
  sigma_f2 <- 0.35
  cov_y_f  <- 0.12
  var_res  <- sigma_y2 - sigma_f2

  metrics <- list(
    type   = "continuous",
    mse    = var_res,
    var_y  = sigma_y2,
    cov_y_f = cov_y_f,
    m_obs  = 200
  )

  direct <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    type = "mean",
    lambda_mode = "vanilla",
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = sigma_f2,
    var_res = var_res,
    correction = FALSE,
    warn_smallN = FALSE
  )

  via_metrics <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    type = "mean",
    lambda_mode = "vanilla",
    metrics = metrics,
    metric_type = "continuous",
    correction = FALSE,
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})

test_that("mean-mode: binary metrics (confusion matrix & PR) match direct moments", {
  delta <- 0.22
  N <- 1800
  alpha <- 0.05
  power <- 0.9

  tp <- 96; fp <- 24; fn <- 30; tn <- 90
  m_obs <- tp + fp + fn + tn

  p_y <- (tp + fn) / m_obs
  p_hat <- (tp + fp) / m_obs

  accuracy <- (tp + tn) / m_obs
  adj <- m_obs / (m_obs - 1)
  bias <- p_y - p_hat

  var_res <- adj * ((1 - accuracy) - bias^2)
  var_f   <- adj * p_hat * (1 - p_hat)

  sigma_y2 <- p_y * (1 - p_y)
  sigma_f2 <- var_f
  cov_y_f  <- tp/m_obs - p_y * p_hat

  # direct version
  direct <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean",
    lambda_mode = "vanilla",
    sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f, var_f = var_f, var_res = var_res,
    warn_smallN = FALSE
  )

  # confusion-matrix version
  metrics_conf <- list(
    type = "classification",
    tp = tp, fp = fp, fn = fn, tn = tn,
    p_y = p_y, var_y = sigma_y2, m_obs = m_obs
  )

  via_conf <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean", lambda_mode = "vanilla",
    metrics = metrics_conf, metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_conf, direct)

  # PR version
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)

  metrics_pr <- list(
    type = "classification",
    precision = precision,
    recall = recall,
    p_y = p_y, var_y = sigma_y2, m_obs = m_obs
  )

  via_pr <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean", lambda_mode = "vanilla",
    metrics = metrics_pr, metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_pr, direct)
})

test_that("mean-mode: binary metrics (confusion matrix & PR) match direct moments", {
  delta <- 0.22
  N <- 1800
  alpha <- 0.05
  power <- 0.9

  tp <- 96; fp <- 24; fn <- 30; tn <- 90
  m_obs <- tp + fp + fn + tn

  p_y <- (tp + fn) / m_obs
  p_hat <- (tp + fp) / m_obs

  accuracy <- (tp + tn) / m_obs
  adj <- m_obs / (m_obs - 1)
  bias <- p_y - p_hat

  var_res <- adj * ((1 - accuracy) - bias^2)
  var_f   <- adj * p_hat * (1 - p_hat)

  sigma_y2 <- p_y * (1 - p_y)
  sigma_f2 <- var_f
  cov_y_f  <- tp/m_obs - p_y * p_hat

  # direct version
  direct <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean",
    lambda_mode = "vanilla",
    sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f, var_f = var_f, var_res = var_res,
    warn_smallN = FALSE
  )

  # confusion-matrix version
  metrics_conf <- list(
    type = "classification",
    tp = tp, fp = fp, fn = fn, tn = tn,
    p_y = p_y, var_y = sigma_y2, m_obs = m_obs
  )

  via_conf <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean", lambda_mode = "vanilla",
    metrics = metrics_conf, metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_conf, direct)

  # PR version
  precision <- tp / (tp + fp)
  recall    <- tp / (tp + fn)

  metrics_pr <- list(
    type = "classification",
    precision = precision,
    recall = recall,
    p_y = p_y, var_y = sigma_y2, m_obs = m_obs
  )

  via_pr <- n_required_ppi_pp(
    delta = delta, N = N, alpha = alpha, power = power,
    type = "mean", lambda_mode = "vanilla",
    metrics = metrics_pr, metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_pr, direct)
})

test_that("mean-mode: infeasible case is capped with achieved power", {
  delta <- 0.03
  N <- 150
  alpha <- 0.05
  power <- 0.9
  sigma_y2 <- 1.0
  sigma_f2 <- 0.25
  cov_y_f  <- 0.05
  var_f <- sigma_f2
  var_res <- sigma_y2 - sigma_f2

  expect_warning(
    capped <- n_required_ppi_pp(
      delta = delta, N = N, alpha = alpha, power = power,
      type = "mean", lambda_mode = "vanilla",
      sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f, var_f = var_f, var_res = var_res,
      warn_smallN = FALSE, mode = "cap"
    ),
    regexp = "capping to n = N"
  )

  expect_identical(as.integer(capped), as.integer(N))
  achieved <- attr(capped, "achieved_power")
  expect_true(is.numeric(achieved) && achieved > 0 && achieved < power)
})

test_that("mean-mode: invalid or inconsistent moments produce errors", {
  expect_error(
    n_required_ppi_pp(
      delta = 0.2, N = 500, alpha = 0.05, power = 0.8,
      type = "mean", lambda_mode = "vanilla",
      sigma_y2 = 1.0, sigma_f2 = 0, cov_y_f = 0.1,
      var_f = 0, var_res = 1.0
    ),
    regexp = "sigma_f2 must be positive"
  )

  expect_error(
    n_required_ppi_pp(
      delta = 0.2, N = 500, alpha = 0.05, power = 0.8,
      type = "mean", lambda_mode = "vanilla",
      sigma_y2 = 1.0, sigma_f2 = 0.3,
      cov_y_f = 5, var_f = 0.3, var_res = 0.7
    ),
    regexp = "Infeasible"
  )
})

test_that("OLS vanilla lambda matches numeric root of variance equation", {
  delta <- 0.15
  N <- 3000
  alpha <- 0.05
  power <- 0.8

  c_vec <- c(1, 0, -1)
  H <- diag(3)
  Sigma_YY <- diag(3)
  Sigma_ff <- 0.4 * diag(3)
  Sigma_Yf <- 0.25 * diag(3)

  lambda <- 1

  var_fun <- function(n) {
    H_mix <- (1 - lambda) * H + lambda * H
    Hinv  <- solve(H_mix)

    middle <- Sigma_YY / n +
      lambda^2 * (Sigma_ff / N + Sigma_ff / n) -
      2 * lambda * (Sigma_Yf / n)

    as.numeric(t(c_vec) %*% Hinv %*% middle %*% Hinv %*% c_vec)
  }

  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2

  obj <- function(n) var_fun(n) - S2
  expected_n <- ceiling(uniroot(obj, c(2, 1e6))$root)

  out <- n_required_ppi_pp(
    delta = delta, N = N,
    alpha = alpha, power = power,
    type = "ols",
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sigma_YY,
    Sigma_ff_l = Sigma_ff,
    Sigma_ff_u = Sigma_ff,
    Sigma_Yf = Sigma_Yf,
    correction = FALSE,
    warn_smallN = FALSE
  )

  expect_identical(out, as.integer(expected_n))
})

test_that("OLS user-defined lambda works and matches vanilla when lambda_user = 1", {
  delta <- 0.2
  N <- 1500
  alpha <- 0.05
  power <- 0.8

  c_vec <- c(1, 0, 1)
  H <- diag(3)
  Sigma_YY <- 1.2 * diag(3)
  Sigma_ff <- 0.5 * diag(3)
  Sigma_Yf <- 0.3 * diag(3)

  out_vanilla <- n_required_ppi_pp(
    delta = delta, N = N,
    alpha = alpha, power = power,
    type = "ols",
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sigma_YY,
    Sigma_ff_l = Sigma_ff,
    Sigma_ff_u = Sigma_ff,
    Sigma_Yf = Sigma_Yf,
    warn_smallN = FALSE
  )

  out_user <- n_required_ppi_pp(
    delta = delta, N = N,
    alpha = alpha, power = power,
    type = "ols",
    lambda_mode = "user",
    lambda_user = 1,
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sigma_YY,
    Sigma_ff_l = Sigma_ff,
    Sigma_ff_u = Sigma_ff,
    Sigma_Yf = Sigma_Yf,
    warn_smallN = FALSE
  )

  expect_identical(out_user, out_vanilla)
})

test_that("OLS oracle lambda gives n <= vanilla lambda", {
  delta <- 0.15
  N <- 2000
  alpha <- 0.05
  power <- 0.8

  c_vec <- c(1, 1, 0)
  H <- diag(3)

  # Construct Sigma matrices so that oracle λ* < 1
  Sigma_YY <- diag(3)
  Sigma_ff <- 2 * diag(3)
  Sigma_Yf <- 1 * diag(3)

  n_vanilla <- n_required_ppi_pp(
    delta = delta, N = N,
    alpha = alpha, power = power,
    type = "ols",
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sigma_YY,
    Sigma_ff_l = Sigma_ff,
    Sigma_ff_u = Sigma_ff,
    Sigma_Yf = Sigma_Yf,
    mode = "cap",
    warn_smallN = FALSE
  )

  n_oracle <- n_required_ppi_pp(
    delta = delta, N = N,
    alpha = alpha, power = power,
    type = "ols",
    lambda_mode = "oracle",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sigma_YY,
    Sigma_ff_l = Sigma_ff,
    Sigma_ff_u = Sigma_ff,
    Sigma_Yf = Sigma_Yf,,
    mode = "cap",
    warn_smallN = FALSE
  )

  expect_true(n_oracle <= n_vanilla)
})

test_that("OLS: infeasible request is capped with achieved power", {
  delta <- 0.02        # tiny effect size → infeasible
  N <- 100
  alpha <- 0.05
  power <- 0.9

  c_vec <- c(1,0,0)
  H <- diag(3)
  Sigma_YY <- diag(3)
  Sigma_ff <- diag(3)
  Sigma_Yf <- 0.2 * diag(3)

  expect_warning(
    capped <- n_required_ppi_pp(
      delta = delta, N = N,
      alpha = alpha, power = power,
      type = "ols",
      lambda_mode = "vanilla",
      c = c_vec,
      H_L = H, H_U = H,
      Sigma_YY = Sigma_YY,
      Sigma_ff_l = Sigma_ff,
      Sigma_ff_u = Sigma_ff,
      Sigma_Yf = Sigma_Yf,
      mode = "cap",
      warn_smallN = FALSE
    ),
    regexp = "capping to n = N"
  )

  expect_identical(as.integer(capped), as.integer(N))

  achieved <- attr(capped, "achieved_power")
  expect_true(is.numeric(achieved) && achieved > 0 && achieved < power)
})

test_that("OLS: validation errors on missing or invalid arguments", {
  c_vec <- c(1,0,1)
  H <- diag(3)
  Sigma <- diag(3)

  # Missing required matrices
  expect_error(
    n_required_ppi_pp(
      delta = 0.2, N = 500,
      type = "ols", lambda_mode = "vanilla",
      c = c_vec, H_L = H, H_U = H
    ),
    regexp = "OLS requires"
  )

  # Invalid lambda_mode="user" with missing lambda_user
  expect_error(
    n_required_ppi_pp(
      delta = 0.2, N = 500,
      type = "ols", lambda_mode = "user",
      c = c_vec,
      H_L = H, H_U = H,
      Sigma_YY = Sigma,
      Sigma_ff_l = Sigma,
      Sigma_ff_u = Sigma,
      Sigma_Yf = Sigma
    ),
    regexp = "lambda_user must be provided"
  )
})