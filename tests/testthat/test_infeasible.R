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
