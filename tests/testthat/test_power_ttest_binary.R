test_that("binary two-sample PPI++ power matches Monte Carlo", {
  skip_on_cran()

  cfg <- list(
    p_A = 0.34,
    p_B = 0.26,
    sens = 0.85,
    spec = 0.85,
    n = 200,
    N = 5000,
    alpha = 0.05,
    R = 1500
  )

  sim <- simulate_ppi_ttest_binary(
    R = cfg$R,
    p_A = cfg$p_A,
    p_B = cfg$p_B,
    sens_A = cfg$sens,
    spec_A = cfg$spec,
    sens_B = cfg$sens,
    spec_B = cfg$spec,
    n_A = cfg$n,
    n_B = cfg$n,
    N_A = cfg$N,
    N_B = cfg$N,
    alpha = cfg$alpha,
    seed = 20240701
  )

  diff <- abs(sim$empirical_power - sim$analytical_power)
  tol <- max(3 * sim$mc_se, 0.01)

  expect_true(
    diff < tol,
    info = sprintf(
      "Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      sim$empirical_power, sim$analytical_power, diff, tol
    )
  )
})
