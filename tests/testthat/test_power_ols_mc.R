N <- 1500    # unlabeled sample size
n <- 200     # labeled sample size
alpha <- 0.05
R_sim <- 2000

test_that("power_ppi_ols outputs probability in (0,1)", {
  delta <- 0.2
  V_u <- matrix(c(1.0, 0.2, 0.2, 0.9), ncol = 2)
  V_l <- matrix(c(0.8, 0.1, 0.1, 0.6), ncol = 2)
  c_vec <- c(1, -0.4)

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec)
  expect_true(out > 0 && out < 1)
})

test_that("power_ppi_ols equals alpha when effect is null", {
  delta <- 0
  V_u <- diag(c(0.9, 1.3))
  V_l <- matrix(c(0.7, 0.05, 0.05, 0.5), ncol = 2)
  c_vec <- c(0.5, 0.8)
  alpha <- 0.05

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec, alpha)
  expect_equal(out, alpha, tolerance = 1e-3)
})

test_that("power_ppi_ols approaches one for large effects", {
  delta <- 6
  V_u <- diag(c(1.1, 0.8))
  V_l <- diag(c(0.6, 0.7))
  c_vec <- c(1, 0.5)

  out <- power_ppi_ols(delta, V_u, V_l, N, n, c_vec)
  expect_gt(out, 0.995)
})

test_that("Monte Carlo power aligns with analytical power for PPI and PPI++ OLS", {
  seed_base <- 20240618L
  set.seed(seed_base)

  beta_true <- c(0.7, 0.4, -0.3)  # intercept, x1, x2
  c_vec <- c(0, 1, 0)
  theta0 <- beta_true[2] - 0.15
  delta_pop <- beta_true[2] - theta0

  var_eta <- beta_true[2]^2 + beta_true[3]^2
  sigma_f2 <- 0.6^2
  sigma_y2 <- 1^2

  var_f_pop   <- var_eta + sigma_f2
  var_res_pop <- sigma_y2 + sigma_f2

  H_pop <- diag(3)
  Sigma_YY_pop <- sigma_y2 * H_pop
  Sigma_ff_pop <- sigma_f2 * H_pop
  Sigma_Yf_pop <- matrix(0, nrow = 3, ncol = 3)
  lambda_pop <- 0

  moments_ppi <- list(
    delta    = delta_pop,
    H        = H_pop,
    Sigma_YY = Sigma_YY_pop
  )

  res_ppi <- simulate_power_ppi_ols(
    R = R_sim,
    n = n,
    N = N,
    contrast = c_vec,
    alpha = alpha,
    family = stats::gaussian(),
    moments = moments_ppi,
    seed = seed_base + 1L
  )

  empirical_ppi <- res_ppi$empirical_power
  analytical_ppi <- res_ppi$analytical_power
  diff_ppi <- abs(empirical_ppi - analytical_ppi)
  mc_se_ppi <- sqrt(analytical_ppi * (1 - analytical_ppi) / R_sim)
  tol_ppi <- max(3 * mc_se_ppi, 0.01)

  expect_true(
    diff_ppi < tol_ppi,
    info = sprintf(
      "PPI-OLS: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      empirical_ppi, analytical_ppi, diff_ppi, tol_ppi
    )
  )

  moments_ppiplus <- list(
    delta = delta_pop,
    H_L = H_pop,
    H_U = H_pop,
    Sigma_YY = Sigma_YY_pop,
    Sigma_ff_l = Sigma_ff_pop,
    Sigma_ff_u = Sigma_ff_pop,
    Sigma_Yf = Sigma_Yf_pop,
    lambda = lambda_pop
  )

  res_ppiplus <- simulate_power_ppi_pp_ols(
    R = R_sim,
    n = n,
    N = N,
    contrast = c_vec,
    alpha = alpha,
    family = stats::gaussian(),
    lambda_type = "plugin",
    moments = moments_ppiplus,
    seed = seed_base + 2L
  )

  empirical_ppiplus <- res_ppiplus$empirical_power
  analytical_ppiplus <- res_ppiplus$analytical_power
  diff_ppiplus <- abs(empirical_ppiplus - analytical_ppiplus)
  mc_se_ppiplus <- sqrt(analytical_ppiplus * (1 - analytical_ppiplus) / R_sim)
  tol_ppiplus <- max(3 * mc_se_ppiplus, 0.01)

  expect_true(
    diff_ppiplus < tol_ppiplus,
    info = sprintf(
      "PPI++-OLS: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      empirical_ppiplus, analytical_ppiplus, diff_ppiplus, tol_ppiplus
    )
  )
})

test_that("Monte Carlo agreement holds across OLS parameter grid", {
  set.seed(as.integer(Sys.Date()))

  beta_grid <- c(0.5, 0.35, -0.25, 0.15)
  c_vec <- c(0, 1, 0, 0)

  var_eta_grid <- sum(beta_grid[-1]^2)
  sigma_f2_grid <- 0.7^2
  sigma_y2_grid <- 1.1^2

  var_f_pop <- var_eta_grid + sigma_f2_grid
  var_res_pop <- sigma_y2_grid + sigma_f2_grid

  H_pop <- diag(4)
  Sigma_YY_pop <- sigma_y2_grid * H_pop
  Sigma_ff_pop <- sigma_f2_grid * H_pop
  Sigma_Yf_pop <- matrix(0, nrow = 4, ncol = 4)

  grid <- expand.grid(
    delta = c(0.08, 0.15, 0.25),
    N     = c(1600, 2300),
    n     = c(180, 280),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
  moments_ppi <- list(
    delta    = g$delta,
    H        = H_pop,
    Sigma_YY = Sigma_YY_pop
  )

  moments_ppiplus <- list(
    delta = g$delta,
    H_L = H_pop,
    H_U = H_pop,
    Sigma_YY = Sigma_YY_pop,
    Sigma_ff_l = Sigma_ff_pop,
    Sigma_ff_u = Sigma_ff_pop,
    Sigma_Yf = Sigma_Yf_pop
  )

    res <- simulate_power_ppi_ols(
      R = R_sim,
      n = g$n,
      N = g$N,
      contrast = c_vec,
      alpha = alpha,
      family = stats::gaussian(),
      moments = moments_ppi,
      seed = 1000 + i
    )

    empirical <- res$empirical_power
    analytical <- res$analytical_power
    diff <- abs(empirical - analytical)
    mc_se <- sqrt(analytical * (1 - analytical) / R_sim)

    ok <- (diff < 3 * mc_se) || (diff < 0.03) ||
      (analytical > 0.98 && empirical > 0.98)

    msg <- sprintf(
      "Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, 3*MCSE=%.3f, delta=%.3f, N=%d, n=%d",
      empirical, analytical, diff, 3 * mc_se, g$delta, g$N, g$n
    )
    expect_true(ok, info = msg)

    res_ppiplus <- simulate_power_ppi_pp_ols(
      R = R_sim,
      n = g$n,
      N = g$N,
      contrast = c_vec,
      alpha = alpha,
      family = stats::gaussian(),
      lambda_type = "plugin",
      moments = moments_ppiplus,
      seed = 2000 + i
    )

    empirical_ppiplus <- res_ppiplus$empirical_power
    analytical_ppiplus <- res_ppiplus$analytical_power
    diff_ppiplus <- abs(empirical_ppiplus - analytical_ppiplus)
    mc_se_ppiplus <- sqrt(analytical_ppiplus * (1 - analytical_ppiplus) / R_sim)

    ok_ppiplus <- (diff_ppiplus < 3 * mc_se_ppiplus) || (diff_ppiplus < 0.03) ||
      (analytical_ppiplus > 0.98 && empirical_ppiplus > 0.98)

    msg_ppiplus <- sprintf(
      "PPI++ OLS: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, 3*MCSE=%.3f, delta=%.3f, N=%d, n=%d",
      empirical_ppiplus, analytical_ppiplus, diff_ppiplus, 3 * mc_se_ppiplus, g$delta, g$N, g$n
    )
    expect_true(ok_ppiplus, info = msg_ppiplus)
  }
})

test_that("Monte Carlo agreement holds across OLS parameter grid (Binomial)", {
  set.seed(as.integer(Sys.Date()))

  beta_grid <- c(-0.2, 0.8, -0.4, 0.2)    # logistic coefficients
  c_vec <- c(0, 1, 0, 0)

  # Pick a representative mean probability (e.g., marginal mean)
  p_mean <- 0.6
  var_y_pop <- p_mean * (1 - p_mean)       # Bernoulli variance

  # Prediction variance (variance of f(X) = P(Y=1 | X))
  sigma_f2_grid <- 0.20
  var_f_pop <- sigma_f2_grid

  # Residual variance for Bernoulli: Var(Y - f(X)) = p(1-p) - Var(f(X))
  var_res_pop <- var_y_pop - var_f_pop

  # Fisher information matrix for logistic regression ~ expected Hessian
  p <- length(beta_grid)
  H_pop <- var_y_pop * diag(p)           # approx. E[ X'WX ]
  Sigma_YY_pop  <- var_y_pop * diag(p)
  Sigma_ff_pop  <- sigma_f2_grid * diag(p)
  Sigma_Yf_pop  <- sigma_f2_grid * diag(p) # assume unbiased predictor

  grid <- expand.grid(
    delta = c(0.08, 0.15, 0.25),
    N     = c(1600, 2300),
    n     = c(180, 280),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]

    moments_ppi <- list(
      delta    = g$delta,
      H        = H_pop,
      Sigma_YY = Sigma_YY_pop
    )

    moments_ppiplus <- list(
      delta       = g$delta,
      H_L         = H_pop,
      H_U         = H_pop,
      Sigma_YY    = Sigma_YY_pop,
      Sigma_ff_l  = Sigma_ff_pop,
      Sigma_ff_u  = Sigma_ff_pop,
      Sigma_Yf    = Sigma_Yf_pop
    )

    res <- simulate_power_ppi_ols(
      R = R_sim,
      n = g$n,
      N = g$N,
      contrast = c_vec,
      alpha = alpha,
      family = stats::binomial(),
      moments = moments_ppi,
      seed = 1000 + i
    )

    empirical <- res$empirical_power
    analytical <- res$analytical_power
    diff <- abs(empirical - analytical)
    mc_se <- sqrt(analytical * (1 - analytical) / R_sim)

    ok <- (diff < 3 * mc_se) || (diff < 0.03) ||
      (analytical > 0.98 && empirical > 0.98)

    msg <- sprintf(
      "Binomial PPI OLS: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, 3*MCSE=%.3f, delta=%.3f, N=%d, n=%d",
      empirical, analytical, diff, 3 * mc_se, g$delta, g$N, g$n
    )
    expect_true(ok, info = msg)

    res_ppiplus <- simulate_power_ppi_pp_ols(
      R = R_sim,
      n = g$n,
      N = g$N,
      contrast = c_vec,
      alpha = alpha,
      family = stats::binomial(),
      lambda_type = "plugin",
      moments = moments_ppiplus,
      seed = 2000 + i
    )

    empirical_ppiplus <- res_ppiplus$empirical_power
    analytical_ppiplus <- res_ppiplus$analytical_power
    diff_ppiplus <- abs(empirical_ppiplus - analytical_ppiplus)
    mc_se_ppiplus <- sqrt(analytical_ppiplus * (1 - analytical_ppiplus) / R_sim)

    ok_ppiplus <- (diff_ppiplus < 3 * mc_se_ppiplus) || (diff_ppiplus < 0.03) ||
      (analytical_ppiplus > 0.98 && empirical_ppiplus > 0.98)

    msg_ppiplus <- sprintf(
      "Binomial PPI++ OLS: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, 3*MCSE=%.3f, delta=%.3f, N=%d, n=%d",
      empirical_ppiplus, analytical_ppiplus, diff_ppiplus, 3 * mc_se_ppiplus, g$delta, g$N, g$n
    )
    expect_true(ok_ppiplus, info = msg_ppiplus)
  }
})


test_that("ppi_plus_ols plugin lambda matches plug-in formula and user override", {
  set.seed(as.integer(Sys.time()))
  beta_true <- c(0.4, -0.9)
  x_l <- rnorm(n)
  x_u <- rnorm(N)

  X_l <- model.matrix(~ x_l)
  X_u <- model.matrix(~ x_u)
  linpred_l <- drop(X_l %*% beta_true)
  linpred_u <- drop(X_u %*% beta_true)

  Y_l <- linpred_l + rnorm(n, sd = 1.1)
  f_l <- linpred_l + rnorm(n, sd = 0.7)
  f_u <- linpred_u + rnorm(N, sd = 0.7)

  data_l <- data.frame(y = Y_l, f = f_l, x = x_l)
  data_u <- data.frame(f = f_u, x = x_u)
  contrast <- c(0, 1)  # slope on x

  fit_plus <- ppi_plus_ols_fit(
    y = y,
    f_col = "f",
    formula = ~ x,
    data_l = data_l,
    data_u = data_u,
    lambda_type = "plugin",
    contrast = contrast
  )

  expect_s3_class(fit_plus, c("ppi_plus_ols_fit", "ppi_ols_fit"))
  expect_true(fit_plus$lambda >= 0 && fit_plus$lambda <= 1)

  # manual plug-in lambda (mirrors ppi_plus_ols internal calculation)
  base <- pppower:::ppi_ols(X_l, Y_l, f_l, X_u, f_u)
  theta_init <- base$theta_hat
  H_total <- crossprod(rbind(X_l, X_u)) / (n + N)
  H_inv <- solve(H_total)

  res_Y <- Y_l - drop(X_l %*% theta_init)
  res_f_l <- f_l - drop(X_l %*% theta_init)
  res_f_u <- f_u - drop(X_u %*% theta_init)

  Sigma_YY <- pppower:::meat_matrix(X_l, res_Y)
  Sigma_ff_l <- pppower:::meat_matrix(X_l, res_f_l)
  Sigma_ff_u <- pppower:::meat_matrix(X_u, res_f_u)
  Sigma_ff <- (Sigma_ff_l + Sigma_ff_u) / 2

  WX_Y <- X_l * as.numeric(res_Y)
  WX_f <- X_l * as.numeric(res_f_l)
  Sigma_Yf <- crossprod(WX_Y, WX_f) / n

  numer <- as.numeric(t(contrast) %*% H_inv %*% Sigma_Yf %*% H_inv %*% contrast)
  denom <- (1 + n / N) * as.numeric(t(contrast) %*% H_inv %*% Sigma_ff %*% H_inv %*% contrast)
  lambda_star <- if (denom <= 0) 0 else numer / denom
  lambda_star <- max(min(lambda_star, 1), 0)

  expect_equal(fit_plus$lambda, lambda_star, tolerance = 1e-10)

  # user-specified lambda should override the plug-in result
  fit_user <- ppi_plus_ols_fit(
    y = y,
    f_col = "f",
    formula = ~ x,
    data_l = data_l,
    data_u = data_u,
    lambda = 0.35,
    lambda_type = "user"
  )
  expect_equal(fit_user$lambda, 0.35, tolerance = 1e-12)

  X_l_formula <- stats::model.matrix(~ x, data_l)
  n_l <- nrow(X_l_formula)
  fit_label <- pppower:::ols_fit(X_l_formula, Y_l)
  fit_zero <- ppi_plus_ols_fit(
    y = y,
    f_col = "f",
    formula = ~ x,
    data_l = data_l,
    data_u = data_u,
    lambda = 0,
    lambda_type = "user"
  )

  H_L <- crossprod(X_l_formula) / n_l
  bread_L <- solve(H_L)
  Sigma_YY <- pppower:::meat_matrix(X_l_formula, fit_label$resid)
  expected_vcov_zero <- bread_L %*% (Sigma_YY / n_l) %*% bread_L

  expect_equal(unname(fit_zero$coef), unname(fit_label$coef), tolerance = 1e-10)
  expect_equal(unname(as.matrix(fit_zero$vcov)), unname(expected_vcov_zero), tolerance = 1e-8)

  X_u_formula <- stats::model.matrix(~ x, data_u)         # unlabeled design
  H_U <- crossprod(X_u_formula) / nrow(X_u_formula)       # X_u^TX_u/N
  G_U <- crossprod(X_u_formula, f_u) / nrow(X_u_formula)

  X_l_formula <- stats::model.matrix(~ x, data_l)         # labeled design
  H_L <- crossprod(X_l_formula) / nrow(X_l_formula)
  G_L <- crossprod(X_l_formula, Y_l) / nrow(X_l_formula)
  W_L <- crossprod(X_l_formula, f_l) / nrow(X_l_formula)

  fit_one <- ppi_plus_ols_fit(
    y = y,
    f_col = "f",
    formula = ~ x,
    data_l = data_l,
    data_u = data_u,
    lambda = 1,
    lambda_type = "user"
  )

  expected_coef_one <- drop(solve(H_U, G_L + G_U - W_L))

  expect_equal(unname(fit_one$coef), unname(expected_coef_one), tolerance = 1e-10)
})
