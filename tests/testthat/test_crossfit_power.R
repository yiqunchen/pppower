super_pop_n <- 20000
N <- 1500    # unlabeled sample size
n <- 200     # labeled sample size
alpha <- 0.05
sim <- 20000

test_that("Analytical power equals Monte Carlo power under cross-fitted data", {
  set.seed(as.integer(Sys.time()))
  df_bin <- simulate_crossfit_data(
    n = super_pop_n, p = 5, family = binomial(), K = 5, seed = as.integer(Sys.time())
  )
  theta <- mean(df_bin$y)
  theta0 <- theta - 0.25
  delta <- theta - theta0

  var_f <- var(df_bin$fhat_cf)
  var_res <- var(df_bin$y - df_bin$fhat_cf)

  power_exact <- power_ppi_mean(
  delta = delta,
  var_f = var_f,
  var_res = var_res,
  N = N,
  n = n,
  alpha = alpha
  )
  
  res_mc <- simulate_power_ppi_mean(
    df_bin, N, n, alpha = alpha, R = sim,
    theta0 = theta0, seed = as.integer(Sys.time())
  )

  empirical <- res_mc$empirical_power
  analytical <- res_mc$analytical_power

  diff <- abs(empirical - analytical)
  se_mc <- sqrt(empirical * (1 - empirical) / sim)  # MC standard error
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
    n = super_pop_n, p = 6, family = gaussian(),
    K = 5, seed = as.integer(Sys.time())
  )

  theta  <- mean(df$y)
  theta0 <- theta - 0.25
  delta  <- theta - theta0

  var_f   <- var(df$fhat_cf)
  var_res <- var(df$y - df$fhat_cf)
  cov_yf  <- cov(df$y, df$fhat_cf)
  
  res <- simulate_power(
    delta,
    N = N,
    n = n,
    alpha = alpha,
    R = sim,
    var_f = var_f,
    var_res = var_res,
    cov_y_f = cov_yf
  )

  analytical_pp     <- res["Exact_PP"]
  analytical_ppplus <- res["Exact_PPplus"]
  empirical_pp      <- res["Empirical_PP"]
  empirical_ppplus  <- res["Empirical_PPplus"]

  diff_pp <- abs(empirical_pp - analytical_pp)
  diff_ppplus <- abs(empirical_ppplus - analytical_ppplus)

  se_mc_pp <- sqrt(empirical_pp * (1 - empirical_pp) / sim)
  tol_pp <- max(3 * se_mc_pp, 0.01)
  se_mc_ppplus <- sqrt(empirical_ppplus * (1 - empirical_ppplus) / sim)
  tol_ppplus <- max(3 * se_mc_ppplus, 0.01)

  expect_lt(
    abs(diff_pp), tol_pp,
    sprintf(
      "Gaussian DGP (PP): Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < tol=%.3f",
      empirical_pp, analytical_pp, abs(diff_pp), tol_pp
    )
  )

  expect_lt(
    abs(diff_ppplus), tol_ppplus,
    sprintf(
      "Gaussian DGP (PP++): Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < tol=%.3f",
      empirical_ppplus, analytical_ppplus, abs(diff_ppplus), tol_ppplus
    )
  )
})

test_that("Analytical power approximates Monte Carlo power (cross-fitted Binomial DGP)", {
  set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))

  # Weak logistic signal to get midrange power (avoid near 1)
  df <- simulate_crossfit_data(
    n = super_pop_n, p = 6, family = binomial(),
    K = 5, seed = as.integer(Sys.time())
  )

  theta  <- mean(df$y)
  theta0 <- theta - 0.25
  delta  <- theta - theta0

  var_f   <- var(df$fhat_cf)
  var_res <- var(df$y - df$fhat_cf)
  cov_yf  <- cov(df$y, df$fhat_cf)

  res <- simulate_power(
    delta,
    N = N,
    n = n,
    alpha = alpha,
    R = sim,
    var_f = var_f,
    var_res = var_res,
    cov_y_f = cov_yf
  )

  analytical_pp     <- res["Exact_PP"]
  analytical_ppplus <- res["Exact_PPplus"]
  empirical_pp      <- res["Empirical_PP"]
  empirical_ppplus  <- res["Empirical_PPplus"]

  diff_pp <- abs(empirical_pp - analytical_pp)
  diff_ppplus <- abs(empirical_ppplus - analytical_ppplus)

  se_mc_pp <- sqrt(empirical_pp * (1 - empirical_pp) / sim)
  tol_pp <- max(3 * se_mc_pp, 0.01)
  se_mc_ppplus <- sqrt(empirical_ppplus * (1 - empirical_ppplus) / sim)
  tol_ppplus <- max(3 * se_mc_ppplus, 0.01)

  expect_lt(
    abs(diff_pp), tol_pp,
    sprintf(
      "Binary DGP (PP): Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < tol=%.3f",
      empirical_pp, analytical_pp, abs(diff_pp), tol_pp
    )
  )

  expect_lt(
    abs(diff_ppplus), tol_ppplus,
    sprintf(
      "Binary DGP (PP++): Empirical=%.3f, Analytical=%.3f, |diff|=%.3f < tol=%.3f",
      empirical_ppplus, analytical_ppplus, abs(diff_ppplus), tol_ppplus
    )
  )
})