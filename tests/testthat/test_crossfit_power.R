# Population parameters
N <- 1500    # unlabeled sample size
n <- 200     # labeled sample size
alpha <- 0.05
R_sim <- 2000

test_that("Monte Carlo power aligns with analytical power for PPI and PPI++ mean estimators", {
  set.seed(20240620L)

  delta_pop    <- 0.2
  var_f_pop    <- 0.45
  var_res_pop  <- 0.8
  sigma_y2_pop <- var_f_pop + var_res_pop
  cov_yf_pop   <- var_f_pop
  sigma_f2_pop <- var_f_pop

  # PPI (vanilla)
  moments_ppi <- list(
    delta   = delta_pop,
    var_f   = var_f_pop,
    var_res = var_res_pop
  )

  res_ppi <- simulate_ppi_vanilla_mean(
    R = R_sim,
    n = n,
    N = N,
    alpha = alpha,
    family = stats::gaussian(),
    moments = moments_ppi,
    seed = 20240621L
  )

  empirical_ppi    <- res_ppi$empirical_power
  analytical_ppi   <- res_ppi$analytical_power
  diff_ppi         <- abs(empirical_ppi - analytical_ppi)
  mc_se_ppi        <- sqrt(analytical_ppi * (1 - analytical_ppi) / R_sim)
  tol_ppi          <- max(3 * mc_se_ppi, 0.01)

  expect_true(
    diff_ppi < tol_ppi,
    info = sprintf(
      "PPI-mean: Empirical=%.4f, Analytical=%.4f, |diff|=%.4f, tol=%.4f",
      empirical_ppi, analytical_ppi, diff_ppi, tol_ppi
    )
  )

  # PPI++
  moments_ppiplus <- list(
    delta     = delta_pop,
    sigma_y2  = sigma_y2_pop,
    sigma_f2  = sigma_f2_pop,
    cov_y_f   = cov_yf_pop,
    var_f     = var_f_pop,
    var_res   = var_res_pop
  )

  res_ppiplus <- simulate_ppi_mean(
    R = R_sim,
    n = n,
    N = N,
    alpha = alpha,
    family = stats::gaussian(),
    moments = moments_ppiplus,
    lambda_type = "plugin",
    seed = 20240622L
  )

  empirical_ppiplus    <- res_ppiplus$empirical_power
  analytical_ppiplus   <- res_ppiplus$analytical_power
  diff_ppiplus         <- abs(empirical_ppiplus - analytical_ppiplus)
  mc_se_ppiplus        <- sqrt(analytical_ppiplus * (1 - analytical_ppiplus) / R_sim)
  tol_ppiplus          <- max(3 * mc_se_ppiplus, 0.01)

  expect_true(
    diff_ppiplus < tol_ppiplus,
    info = sprintf(
      "PPI++-mean: Empirical=%.4f, Analytical=%.4f, |diff|=%.4f, tol=%.4f",
      empirical_ppiplus, analytical_ppiplus, diff_ppiplus, tol_ppiplus
    )
  )
})
