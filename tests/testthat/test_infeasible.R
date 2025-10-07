test_that("n_required_PP correctly detects infeasible configurations", {

  # Case where unlabeled term dominates: var_f / N >= delta^2 / c0^2, also for small N
  delta  <- 0.1
  var_f  <- 2.0
  var_res <- 1.0
  N <- 100
  alpha <- 0.05
  power <- 0.8

  expect_error(
    n_required_PP(
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
    n_required_PP(
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
    n_required_PP(
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
    n_required_PP(
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


test_that("n_required_PP returns expected sample size for the mean estimator", {
  delta <- 0.2
  N <- 1000
  alpha <- 0.05
  power <- 0.8
  var_f <- 0.3
  var_res <- 1.1

  c0 <- stats::qnorm(1 - alpha / 2) + stats::qnorm(power)
  denom <- (delta^2 / c0^2) - (var_f / N)
  expected_n <- as.integer(ceiling(var_res / denom))

  result <- n_required_PP(
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

test_that("n_required_PP matches analytical solution for PP-OLS contrasts", {
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

  result <- n_required_PP(
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

test_that("n_required_PP caps infeasible requests and records achieved power", {
  delta <- 0.05
  N <- 120
  alpha <- 0.05
  power <- 0.9
  var_unlabeled <- 0.01
  var_labeled <- 2

  expect_warning(
    capped <- n_required_PP(
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

test_that("n_required_PP validates required inputs by type", {
  expect_error(
    n_required_PP(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_res = 1
    ),
    regexp = "var_f"
  )

  expect_error(
    n_required_PP(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_f = -0.1,
      var_res = 1
    ),
    regexp = "var_f must be >= 0"
  )

  expect_error(
    n_required_PP(
      delta = 0.2,
      N = 600,
      type = "mean",
      var_f = 0.1,
      var_res = 0
    ),
    regexp = "var_res must be > 0"
  )

  expect_error(
    n_required_PP(
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
    n_required_PP(
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
    n_required_PP(
      delta = 0.2,
      N = 600,
      type = "custom",
      var_unlabeled = -1,
      var_labeled = 1
    ),
    regexp = "var_unlabeled must be non-negative"
  )

  expect_error(
    n_required_PP(
      delta = 0.2,
      N = 600,
      type = "custom",
      var_unlabeled = 1,
      var_labeled = 0
    ),
    regexp = "var_labeled must be positive"
  )
})
