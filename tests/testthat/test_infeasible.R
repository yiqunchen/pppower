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

test_that("n_required_ppi_pp matches the PPI++ quadratic solution", {
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
    N = N,
    alpha = alpha,
    power = power,
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = var_f,
    var_res = var_res,
    warn_smallN = FALSE
  )

  expect_identical(result, as.integer(expected_n))
})

test_that("n_required_ppi_pp works with metrics input", {
  delta <- 0.18
  N <- 2500
  alpha <- 0.05
  power <- 0.85
  sigma_y2 <- 0.9
  sigma_f2 <- 0.35
  cov_y_f  <- 0.12
  var_res  <- sigma_y2 - sigma_f2

  metrics <- list(
    type = "continuous",
    mse = var_res,
    var_y = sigma_y2,
    cov_y_f = cov_y_f,
    m_obs = 200
  )

  direct <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = sigma_f2,
    var_res = var_res,
    warn_smallN = FALSE
  )

  via_metrics <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    metrics = metrics,
    metric_type = "continuous",
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})

test_that("n_required_ppi_pp works with binary metrics input", {
  delta <- 0.22
  N <- 1800
  alpha <- 0.05
  power <- 0.9
  p_y <- 0.4
  p_hat <- 0.42
  accuracy <- 0.78
  m_obs <- 240
  adj <- m_obs / (m_obs - 1)
  bias <- p_y - p_hat
  var_res <- adj * ((1 - accuracy) - bias^2)
  var_f <- adj * p_hat * (1 - p_hat)
  sigma_y2 <- p_y * (1 - p_y)
  sigma_f2 <- var_f
  tp_rate <- (p_hat + p_y + accuracy - 1) / 2
  cov_y_f <- tp_rate - p_y * p_hat

  metrics <- list(
    type = "hard",
    accuracy = accuracy,
    p_y = p_y,
    p_hat = p_hat,
    var_y = sigma_y2,
    m_obs = m_obs
  )

  direct <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    sigma_y2 = sigma_y2,
    sigma_f2 = sigma_f2,
    cov_y_f = cov_y_f,
    var_f = var_f,
    var_res = var_res,
    warn_smallN = FALSE
  )

  via_metrics <- n_required_ppi_pp(
    delta = delta,
    N = N,
    alpha = alpha,
    power = power,
    metrics = metrics,
    metric_type = "hard",
    warn_smallN = FALSE
  )

  expect_identical(via_metrics, direct)
})

test_that("n_required_ppi_pp caps infeasible requests and records achieved power", {
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
      delta = delta,
      N = N,
      alpha = alpha,
      power = power,
      sigma_y2 = sigma_y2,
      sigma_f2 = sigma_f2,
      cov_y_f = cov_y_f,
      var_f = var_f,
      var_res = var_res,
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

test_that("n_required_ppi_pp validates required moments", {
  expect_error(
    n_required_ppi_pp(
      delta = 0.2,
      N = 500,
      alpha = 0.05,
      power = 0.8,
      sigma_y2 = 1.0,
      sigma_f2 = 0,
      cov_y_f = 0.1,
      var_f = 0,
      var_res = 1.0
    ),
    regexp = "sigma_f2 must be positive"
  )

  expect_error(
    n_required_ppi_pp(
      delta = 0.2,
      N = 500,
      alpha = 0.05,
      power = 0.8,
      sigma_y2 = 1.0,
      sigma_f2 = 0.3,
      cov_y_f = 5,
      var_f = 0.3,
      var_res = 0.7
    ),
    regexp = "Infeasible"
  )
})
