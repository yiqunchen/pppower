test_that("power_ppi_ols outputs probability in (0,1)", {
  delta <- 0.2
  V_u <- matrix(c(1.0, 0.2, 0.2, 0.9), ncol = 2)
  V_l <- matrix(c(0.8, 0.1, 0.1, 0.6), ncol = 2)
  N <- 1500
  n <- 250
  c_vec <- c(1, -0.4)

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec)
  expect_true(out > 0 && out < 1)
})

test_that("power_ppi_ols equals alpha when effect is null", {
  delta <- 0
  V_u <- diag(c(0.9, 1.3))
  V_l <- matrix(c(0.7, 0.05, 0.05, 0.5), ncol = 2)
  N <- 1800
  n <- 200
  c_vec <- c(0.5, 0.8)
  alpha <- 0.05

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec, alpha)
  expect_equal(out, alpha, tolerance = 1e-3)
})

test_that("power_ppi_ols approaches one for large effects", {
  delta <- 6
  V_u <- diag(c(1.1, 0.8))
  V_l <- diag(c(0.6, 0.7))
  N <- 1200
  n <- 180
  c_vec <- c(1, 0.5)

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec)
  expect_gt(out, 0.995)
})

test_that("Monte Carlo power aligns with analytical power for PPI-OLS", {
  skip_on_cran()
  set.seed(20240501)

  n_total <- 4500
  x1 <- rnorm(n_total)
  x2 <- rnorm(n_total)
  eta <- 0.7 + 0.4 * x1 - 0.3 * x2
  y <- eta + rnorm(n_total, sd = 1.0)
  fhat_cf <- eta + rnorm(n_total, sd = 0.6)

  df <- data.frame(y = y, fhat_cf = fhat_cf, x1 = x1, x2 = x2)
  formula <- ~ x1 + x2
  X_full <- stats::model.matrix(formula, df)
  theta_full <- ols_fit(X_full, df$y)$coef

  c_vec <- c(0, 1, 0)  # coefficient on x1
  theta0 <- as.numeric(theta_full["x1"]) - 0.15

  N <- 2200
  n <- 260
  alpha <- 0.05
  R_sim <- 1200

  res <- simulate_power_ppi_ols(
    df = df,
    formula = formula,
    N = N,
    n = n,
    c = c_vec,
    theta0 = theta0,
    alpha = alpha,
    R = R_sim,
    seed = 321
  )

  empirical <- res$empirical_power
  analytical <- res$analytical_power
  diff <- abs(empirical - analytical)
  mc_se <- sqrt(analytical * (1 - analytical) / R_sim)
  tol <- max(3 * mc_se, 0.02)

  expect_true(
    diff < tol,
    info = sprintf(
      "Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      empirical, analytical, diff, tol
    )
  )
})

test_that("Monte Carlo agreement holds across OLS parameter grid", {
  skip_on_cran()
  set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))

  n_total <- 5000
  x1 <- rnorm(n_total)
  x2 <- rnorm(n_total)
  x3 <- rnorm(n_total)
  eta <- 0.5 + 0.35 * x1 - 0.25 * x2 + 0.15 * x3
  y <- eta + rnorm(n_total, sd = 1.1)
  fhat_cf <- eta + rnorm(n_total, sd = 0.7)

  df <- data.frame(y = y, fhat_cf = fhat_cf, x1 = x1, x2 = x2, x3 = x3)
  formula <- ~ x1 + x2 + x3
  X_full <- stats::model.matrix(formula, df)
  theta_full <- ols_fit(X_full, df$y)$coef

  c_vec <- c(0, 1, 0, 0)  # focus on x1 coefficient
  alpha <- 0.05
  R_sim <- 800

  grid <- expand.grid(
    delta = c(0.08, 0.15, 0.25),
    N = c(1600, 2300),
    n = c(180, 280),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    theta0 <- as.numeric(theta_full["x1"]) - g$delta

    res <- simulate_power_ppi_ols(
      df = df,
      formula = formula,
      N = g$N,
      n = g$n,
      c = c_vec,
      theta0 = theta0,
      alpha = alpha,
      R = R_sim,
      seed = 1000 + i
    )

    empirical <- res$empirical_power
    analytical <- res$analytical_power
    diff <- abs(empirical - analytical)
    mc_se <- sqrt(analytical * (1 - analytical) / R_sim)

    ok <- (diff < 4 * mc_se) || (diff < 0.03) ||
      (analytical > 0.98 && empirical > 0.98)

    msg <- sprintf(
      "Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, 4*MCSE=%.3f, delta=%.3f, N=%d, n=%d",
      empirical, analytical, diff, 4 * mc_se, g$delta, g$N, g$n
    )
    expect_true(ok, info = msg)
  }
})
