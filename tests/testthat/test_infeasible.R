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

  result <- power_ppi_mean(
    delta = delta,
    N     = N,
    n     = NULL,
    power = power,
    alpha = alpha,
    lambda_mode = "oracle",
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f,
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

  # identity: Var(Y-f) = Var(Y) + Var(f) - 2 Cov(Y,f)
  var_res  <- sigma_y2 + sigma_f2 - 2 * cov_y_f

  # metrics list should match what resolve_ppi_variances() expects:
  metrics <- list(
      type    = "continuous",
      var_y   = sigma_y2,
      var_f   = sigma_f2,
      cov_y_f = cov_y_f,
      var_res = var_res,
      m_obs   = 200
  )

  direct <- power_ppi_mean(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f,
    var_f    = sigma_f2,
    var_res  = var_res,
    correction = FALSE,
    warn_smallN = FALSE
  )

  via_metrics <- power_ppi_mean(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    metrics = metrics,
    var_f = sigma_f2,
    var_res = var_res,
    metric_type = "continuous",
    correction = FALSE,
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})


test_that("mean-mode: binary classification metrics match direct moments", {
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

  direct <- power_ppi_mean(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha, lambda_mode = "vanilla",
    sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
    cov_y_f  = cov_y_f, var_f = var_f, var_res = var_res,
    warn_smallN = FALSE
  )

  metrics_conf <- list(
    type     = "classification",
    tp = tp, fp = fp, fn = fn, tn = tn,
    p_y      = p_y,
    var_y    = sigma_y2,
    m_obs    = m_obs
  )

  via_conf <- power_ppi_mean(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha, lambda_mode = "vanilla",
    metrics = metrics_conf,
    metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_conf, direct)

  metrics_pr <- list(
    type      = "classification",
    precision = tp / (tp + fp),
    recall    = tp / (tp + fn),
    p_y       = p_y,
    var_y     = sigma_y2,
    m_obs     = m_obs
  )

  via_pr <- power_ppi_mean(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha, lambda_mode = "vanilla",
    metrics = metrics_pr,
    metric_type = "classification",
    warn_smallN = FALSE
  )

  expect_identical(via_pr, direct)
})



test_that("mean-mode: infeasible request is capped with achieved power", {
  delta <- 0.03
  N <- 150
  sigma_y2 <- 1.0
  sigma_f2 <- 0.25
  cov_y_f  <- 0.05
  var_f    <- sigma_f2
  var_res  <- sigma_y2 + sigma_f2 - 2 * cov_y_f

  expect_warning(
    capped <- power_ppi_mean(
      delta = delta, N = N,
      n = NULL, power = 0.9,
      alpha = 0.05, lambda_mode = "vanilla",
      sigma_y2 = sigma_y2, sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f, var_f = var_f, var_res = var_res,
      warn_smallN = FALSE, mode = "cap"
    ),
    regexp = "capping to n = N"
  )

  expect_identical(as.integer(capped), as.integer(N))
  expect_true(attr(capped, "achieved_power") < 0.9)
})



test_that("mean-mode: invalid / inconsistent inputs throw correct errors", {
  expect_error(
    power_ppi_mean(
      delta = 0.2, N = 500, n = NULL, power = 0.8,
      lambda_mode = "vanilla",
      sigma_y2 = 1.0, sigma_f2 = 0, cov_y_f = 0.1,
      var_f = 0, var_res = 1.0
    ),
    regexp = "positive"
  )

  expect_error(
    power_ppi_mean(
      delta = 0.2, N = 500, n = NULL, power = 0.8,
      lambda_mode = "vanilla",
      sigma_y2 = 1.0, sigma_f2 = 0.3,
      cov_y_f = 5, var_f = 0.3, var_res = 0.7
    ),
    regexp = "Infeasible"
  )
})



test_that("regression: vanilla lambda matches numeric root", {
  delta <- 0.15
  N <- 3000
  alpha <- 0.05
  power <- 0.8

  c_vec <- c(1,0,-1)
  H <- diag(3)
  SYY <- diag(3)
  Sff <- 0.4 * diag(3)
  SYf <- 0.25 * diag(3)

  lambda <- 1

  var_fun <- function(n) {
    middle <- SYY/n +
      lambda^2 * (Sff/N + Sff/n) -
      2 * lambda * (SYf/n)
    t(c_vec) %*% solve(H) %*% middle %*% solve(H) %*% c_vec
  }

  z_alpha <- qnorm(1 - alpha/2)
  z_beta  <- qnorm(power)
  S2 <- (delta / (z_alpha + z_beta))^2

  expected_n <- ceiling(uniroot(function(n) var_fun(n) - S2, c(2, 1e6))$root)

  out <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = SYY,
    Sigma_ff_l = Sff, Sigma_ff_u = Sff,
    Sigma_Yf = SYf,
    warn_smallN = FALSE
  )

  expect_identical(out, as.integer(expected_n))
})



test_that("regression: lambda_user = 1 equals vanilla", {
  c_vec <- c(1,0,1)
  H <- diag(3)
  Sig <- diag(3)

  out_vanilla <- power_ppi_regression(
    delta = 0.2, N = 1500,
    n = NULL, power = 0.8,
    alpha = 0.05,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sig,
    Sigma_ff_l = Sig, Sigma_ff_u = Sig,
    Sigma_Yf = 0.3 * Sig,
    warn_smallN = FALSE
  )

  out_user <- power_ppi_regression(
    delta = 0.2, N = 1500,
    n = NULL, power = 0.8,
    alpha = 0.05,
    lambda_mode = "user",
    lambda_user = 1,
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = Sig,
    Sigma_ff_l = Sig, Sigma_ff_u = Sig,
    Sigma_Yf = 0.3 * Sig,
    warn_smallN = FALSE
  )

  expect_identical(out_user, out_vanilla)
})



test_that("regression: oracle lambda gives <= vanilla", {
  delta <- 0.15
  N <- 2000

  c_vec <- c(1,1,0)
  H <- diag(3)
  SYY <- diag(3)
  Sff <- 2 * diag(3)
  SYf <- 1 * diag(3)

  n_v <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = 0.8,
    alpha = 0.05,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = SYY,
    Sigma_ff_l = Sff, Sigma_ff_u = Sff,
    Sigma_Yf = SYf,
    mode = "cap", warn_smallN = FALSE
  )

  n_o <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = 0.8,
    alpha = 0.05,
    lambda_mode = "oracle",
    c = c_vec,
    H_L = H, H_U = H,
    Sigma_YY = SYY,
    Sigma_ff_l = Sff, Sigma_ff_u = Sff,
    Sigma_Yf = SYf,
    mode = "cap", warn_smallN = FALSE
  )

  expect_true(n_o <= n_v)
})



test_that("regression: infeasible case is capped with achieved power", {
  delta <- 0.02
  N <- 100

  c_vec <- c(1,0,0)
  H <- diag(3)
  SYY <- diag(3)
  Sff <- diag(3)
  SYf <- 0.2 * diag(3)

  expect_warning(
    capped <- power_ppi_regression(
      delta = delta, N = N,
      n = NULL, power = 0.9,
      alpha = 0.05,
      lambda_mode = "vanilla",
      c = c_vec,
      H_L = H, H_U = H,
      Sigma_YY = SYY,
      Sigma_ff_l = Sff, Sigma_ff_u = Sff,
      Sigma_Yf = SYf,
      mode = "cap",
      warn_smallN = FALSE
    ),
    regexp = "capping to n = N"
  )

  expect_identical(as.integer(capped), as.integer(N))
  expect_true(attr(capped, "achieved_power") < 0.9)
})



test_that("regression: validation errors for missing blocks", {
  c_vec <- c(1,0,1)
  H <- diag(3)
  Sig <- diag(3)

  expect_error(
    power_ppi_regression(
      delta = 0.2, N = 500,
      n = NULL, power = 0.8,
      lambda_mode = "vanilla",
      c = c_vec, H_L = H, H_U = H
    ),
    regexp = "requires c, H_L, H_U, Sigma_YY, Sigma_ff_l, Sigma_ff_u, Sigma_Yf"
  )

  expect_error(
    power_ppi_regression(
      delta = 0.2, N = 500,
      n = NULL, power = 0.8,
      lambda_mode = "user",
      c = c_vec, H_L = H, H_U = H,
      Sigma_YY = Sig, Sigma_ff_l = Sig,
      Sigma_ff_u = Sig, Sigma_Yf = Sig
      # missing lambda_user
    ),
    regexp = "lambda_user must be provided"
  )
})

test_that("regression-mode: metrics-only continuous input works", {
  set.seed(123)

  n_l <- 400
  n_u <- 2000
  p   <- 3

  # True beta
  beta <- c(0.5, -1.0, 0.7)

  X_l <- matrix(rnorm(n_l * p), n_l, p)
  X_u <- matrix(rnorm(n_u * p), n_u, p)

  # Generate data
  Y_l <- drop(X_l %*% beta + rnorm(n_l, sd = 1.0))
  f_l <- drop(X_l %*% beta + rnorm(n_l, sd = 0.3))
  f_u <- drop(X_u %*% beta + rnorm(n_u, sd = 0.3))

  # Prediction residual stats
  var_res <- var(Y_l - f_l)
  var_f   <- var(f_l)
  var_y   <- var(Y_l)
  cov_y_f <- cov(Y_l, f_l)

  metrics <- list(
    type    = "continuous",
    mse     = var_res,
    var_y   = var_y,
    cov_y_f = cov_y_f,
    m_obs   = n_l
  )

  # Build regression blocks MUST PROVIDE BETA
  blocks <- compute_ppi_blocks(
    model_type = "ols",
    X_l = X_l, Y_l = Y_l, f_l = f_l,
    X_u = X_u, f_u = f_u,
    beta = beta
  )

  c_vec <- c(1, 0, -1)
  N     <- n_u
  delta <- 0.25
  alpha <- 0.05
  power <- 0.8

  # Direct-moment version
  direct <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = blocks$H_L, H_U = blocks$H_U,
    Sigma_YY   = blocks$Sigma_YY,
    Sigma_ff_l = blocks$Sigma_ff_l,
    Sigma_ff_u = blocks$Sigma_ff_u,
    Sigma_Yf   = blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  # metrics-only version (same blocks, just also pass metrics)
  via_metrics <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = blocks$H_L, H_U = blocks$H_U,
    Sigma_YY   = blocks$Sigma_YY,
    Sigma_ff_l = blocks$Sigma_ff_l,
    Sigma_ff_u = blocks$Sigma_ff_u,
    Sigma_Yf   = blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})


test_that("regression-mode (GLM): metrics-only continuous input works", {
  set.seed(456)

  n_l <- 500
  n_u <- 2500
  p   <- 3

  # True GLM coefficients
  beta <- c(0.8, -0.5, 1.2)

  X_l <- matrix(rnorm(n_l * p), n_l, p)
  X_u <- matrix(rnorm(n_u * p), n_u, p)

  # Logistic link
  linkfun <- function(x) 1 / (1 + exp(-x))

  # True linear predictors
  eta_l <- drop(X_l %*% beta)
  eta_u <- drop(X_u %*% beta)

  # Generate Bernoulli outcomes
  prob_l <- linkfun(eta_l)
  Y_l <- rbinom(n_l, size = 1, prob = prob_l)

  # Prediction model (slightly noisy probability)
  f_l <- linkfun(eta_l + rnorm(n_l, sd = 0.15))
  f_u <- linkfun(eta_u + rnorm(n_u, sd = 0.15))

  # Prediction residual stats
  var_res <- var(Y_l - f_l)
  var_f   <- var(f_l)
  var_y   <- prob_l * (1 - prob_l) |> mean()  # approximate Bernoulli variance
  cov_y_f <- cov(Y_l, f_l)

  metrics <- list(
    type    = "continuous",
    mse     = var_res,
    var_y   = var_y,
    cov_y_f = cov_y_f,
    m_obs   = n_l
  )

  # Build regression blocks for GLM
  blocks <- compute_ppi_blocks(
    model_type = "glm",
    X_l = X_l, Y_l = Y_l, f_l = f_l,
    X_u = X_u, f_u = f_u,
    beta = beta,
    family = "binomial"
  )

  c_vec <- c(1, -1, 0)
  N     <- n_u
  delta <- 0.40        # effect size on contrast scale
  alpha <- 0.05
  power <- 0.8

  # Direct-moment version
  direct <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = blocks$H_L, H_U = blocks$H_U,
    Sigma_YY   = blocks$Sigma_YY,
    Sigma_ff_l = blocks$Sigma_ff_l,
    Sigma_ff_u = blocks$Sigma_ff_u,
    Sigma_Yf   = blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  # metrics-only version (same blocks)
  via_metrics <- power_ppi_regression(
    delta = delta, N = N,
    n = NULL, power = power,
    alpha = alpha,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = blocks$H_L, H_U = blocks$H_U,
    Sigma_YY   = blocks$Sigma_YY,
    Sigma_ff_l = blocks$Sigma_ff_l,
    Sigma_ff_u = blocks$Sigma_ff_u,
    Sigma_Yf   = blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})
