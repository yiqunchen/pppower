set.seed(20260307)

make_ols_blocks <- function(n = 240, N = 1600) {
  p <- 3
  beta <- c(0.2, 0.5, -0.3)

  X_l <- cbind(1, matrix(stats::rnorm(n * (p - 1)), nrow = n))
  X_u <- cbind(1, matrix(stats::rnorm(N * (p - 1)), nrow = N))

  y_l <- as.numeric(X_l %*% beta + stats::rnorm(n, sd = 1.0))
  f_l <- as.numeric(X_l %*% beta + stats::rnorm(n, sd = 0.7))
  f_u <- as.numeric(X_u %*% beta + stats::rnorm(N, sd = 0.7))

  list(
    N = N,
    n = n,
    beta = beta,
    X_l = X_l,
    Y_l = y_l,
    f_l = f_l,
    X_u = X_u,
    f_u = f_u,
    blocks = compute_ppi_blocks(
      model_type = "ols",
      X_l = X_l,
      Y_l = y_l,
      f_l = f_l,
      X_u = X_u,
      f_u = f_u,
      beta = beta
    )
  )
}

test_that("power_ppi_regression returns valid probabilities and grows with effect size", {
  dat <- make_ols_blocks()
  c_vec <- c(0, 1, 0)

  p_small <- power_ppi_regression(
    delta = 0.10,
    N = dat$N,
    n = dat$n,
    alpha = 0.05,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = dat$blocks$H_L,
    H_U = dat$blocks$H_U,
    Sigma_YY = dat$blocks$Sigma_YY,
    Sigma_ff_l = dat$blocks$Sigma_ff_l,
    Sigma_ff_u = dat$blocks$Sigma_ff_u,
    Sigma_Yf = dat$blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  p_large <- power_ppi_regression(
    delta = 0.35,
    N = dat$N,
    n = dat$n,
    alpha = 0.05,
    lambda_mode = "vanilla",
    c = c_vec,
    H_L = dat$blocks$H_L,
    H_U = dat$blocks$H_U,
    Sigma_YY = dat$blocks$Sigma_YY,
    Sigma_ff_l = dat$blocks$Sigma_ff_l,
    Sigma_ff_u = dat$blocks$Sigma_ff_u,
    Sigma_Yf = dat$blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  expect_true(is.finite(p_small) && p_small > 0 && p_small < 1)
  expect_true(is.finite(p_large) && p_large > 0 && p_large <= 1)
  expect_gt(p_large, p_small)
})

test_that("required n from power_ppi_regression achieves requested power", {
  c_vec <- c(1, 0, -1)
  H <- diag(3)
  SYY <- diag(3)
  Sff <- 0.5 * diag(3)
  SYf <- 0.2 * diag(3)

  target_power <- 0.8

  n_req <- power_ppi_regression(
    delta = 0.20,
    N = 2400,
    n = NULL,
    power = target_power,
    alpha = 0.05,
    lambda_mode = "oracle",
    c = c_vec,
    H_L = H,
    H_U = H,
    Sigma_YY = SYY,
    Sigma_ff_l = Sff,
    Sigma_ff_u = Sff,
    Sigma_Yf = SYf,
    warn_smallN = FALSE
  )

  achieved <- power_ppi_regression(
    delta = 0.20,
    N = 2400,
    n = n_req,
    power = NULL,
    alpha = 0.05,
    lambda_mode = "oracle",
    c = c_vec,
    H_L = H,
    H_U = H,
    Sigma_YY = SYY,
    Sigma_ff_l = Sff,
    Sigma_ff_u = Sff,
    Sigma_Yf = SYf,
    warn_smallN = FALSE
  )

  expect_gte(n_req, 2)
  expect_true(is.finite(achieved) && achieved >= target_power)
})

test_that("compute_ppi_blocks integrates with binomial GLM regression power", {
  n <- 260
  N <- 1400
  p <- 3
  beta <- c(-0.3, 0.6, -0.4)

  X_l <- cbind(1, matrix(stats::rnorm(n * (p - 1)), nrow = n))
  X_u <- cbind(1, matrix(stats::rnorm(N * (p - 1)), nrow = N))

  prob_l <- stats::plogis(as.numeric(X_l %*% beta))
  prob_u <- stats::plogis(as.numeric(X_u %*% beta))

  Y_l <- stats::rbinom(n, size = 1, prob = prob_l)
  f_l <- pmin(pmax(prob_l + stats::rnorm(n, sd = 0.03), 1e-4), 1 - 1e-4)
  f_u <- pmin(pmax(prob_u + stats::rnorm(N, sd = 0.03), 1e-4), 1 - 1e-4)

  blocks <- compute_ppi_blocks(
    model_type = "glm",
    X_l = X_l,
    Y_l = Y_l,
    f_l = f_l,
    X_u = X_u,
    f_u = f_u,
    beta = beta,
    family = "binomial"
  )

  pwr <- power_ppi_regression(
    delta = 0.18,
    N = N,
    n = n,
    alpha = 0.05,
    lambda_mode = "oracle",
    c = c(0, 1, 0),
    H_L = blocks$H_L,
    H_U = blocks$H_U,
    Sigma_YY = blocks$Sigma_YY,
    Sigma_ff_l = blocks$Sigma_ff_l,
    Sigma_ff_u = blocks$Sigma_ff_u,
    Sigma_Yf = blocks$Sigma_Yf,
    warn_smallN = FALSE
  )

  expect_true(is.finite(pwr) && pwr > 0 && pwr < 1)
})
