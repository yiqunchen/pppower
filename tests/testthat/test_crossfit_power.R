test_that("Analytical power equals Monte Carlo power under cross-fitted data", {
  set.seed(as.integer(Sys.time()))
  df_bin <- simulate_crossfit_data(
    n = 3000, p = 5, family = binomial(), K = 5, seed = as.integer(Sys.time())
  )

  N <- 1500    # unlabeled sample size
  n <- 200     # labeled sample size
  alpha <- 0.05
  theta_true <- mean(df_bin$y)
  theta0 <- theta_true - 0.2
  delta <- theta_true - theta0

  var_f <- var(df_bin$fhat_cf)
  var_res <- var(df_bin$y - df_bin$fhat_cf)

  power_exact <- power_ppi_mean(delta, var_f, var_res, N, n, alpha)

  res_mc <- simulate_power_ppi_mean(
    df_bin, N, n, alpha = alpha, R = 1000,
    theta0 = theta0, seed = as.integer(Sys.time())
  )

  empirical <- res_mc$empirical_power
  analytical <- res_mc$analytical_power

  diff <- abs(empirical - analytical)
  se_mc <- sqrt(empirical * (1 - empirical) / 1000)  # MC standard error
  tol <- max(3 * se_mc, 0.01)
  expect_true(
  diff < tol,
  label = sprintf("Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < tol=%.3f",
                  empirical, analytical, diff, tol)
)
})

test_that("Analytical power approximates Monte Carlo power (cross-fitted Gaussian DGP)", {
  set.seed(as.integer(Sys.time()))

  # Moderate-signal Gaussian DGP to avoid power saturation
  df <- simulate_crossfit_data(
    n = 4000, p = 6, family = gaussian(),
    K = 5, seed = as.integer(Sys.time())
  )

  theta  <- mean(df$y)
  theta0 <- theta - 0.25
  delta  <- theta - theta0

  var_f   <- var(df$fhat_cf)
  var_res <- var(df$y - df$fhat_cf)

  N <- 2000
  n <- 150
  alpha <- 0.05

  analytical <- power_ppi_mean(delta, var_f, var_res, N, n, alpha)
  empirical  <- simulate_power(delta, var_f, var_res, N, n, alpha, R = 50000)[2]

  diff  <- abs(empirical - analytical)
  se_mc <- sqrt(empirical * (1 - empirical) / 50000)

  expect_lt(
    diff, max(3 * se_mc, 0.01),
    sprintf(
      "Gaussian DGP: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < 3*MCSE=%.3f",
      empirical, analytical, diff, 3 * se_mc
    )
  )
})

test_that("Analytical power ≈ Monte Carlo power (cross-fitted Binomial DGP)", {
  set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))

  # Weak logistic signal to get midrange power (avoid near 1)
  df <- simulate_crossfit_data(
    n = 4000, p = 6, family = binomial(),
    K = 5, seed = 20251005
  )

  theta  <- mean(df$y)
  theta0 <- theta - 0.12
  delta  <- theta - theta0

  var_f   <- var(df$fhat_cf)
  var_res <- var(df$y - df$fhat_cf)

  N <- 3000
  n <- 250
  alpha <- 0.05

  analytical <- power_ppi_mean(delta, var_f, var_res, N, n, alpha)
  empirical  <- simulate_power(delta, var_f, var_res, N, n, alpha, R = 40000)[2]

  diff  <- abs(empirical - analytical)
  se_mc <- sqrt(empirical * (1 - empirical) / 40000)

  expect_lt(
    diff, max(3 * se_mc, 0.01),
    sprintf(
      "Binary DGP: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < 3*MCSE=%.3f",
      empirical, analytical, diff, 3 * se_mc
    )
  )
})