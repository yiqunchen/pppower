test_that("EIF binary power matches Monte Carlo", {
  skip_on_cran()

  cfg <- list(
    p = 0.55,
    sens = 0.82,
    spec = 0.84,
    N = 4000,
    m_cal = 400,
    delta = 0.05,
    alpha = 0.05,
    R = 1200
  )

  analytical <- power_eif_binary(
    delta = cfg$delta,
    N = cfg$N,
    m_cal = cfg$m_cal,
    p = cfg$p,
    sens = cfg$sens,
    spec = cfg$spec,
    alpha = cfg$alpha
  )

  sim <- simulate_eif_binary(
    R = cfg$R,
    delta = cfg$delta,
    N = cfg$N,
    m_cal = cfg$m_cal,
    p = cfg$p,
    sens = cfg$sens,
    spec = cfg$spec,
    alpha = cfg$alpha,
    seed = 20240702
  )

  diff <- abs(sim$empirical_power - analytical)
  tol <- max(3 * sim$mc_se, 0.02)

  expect_true(
    diff < tol,
    info = sprintf(
      "Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      sim$empirical_power, analytical, diff, tol
    )
  )
})

test_that("power_eif_binary with m_cal=NULL inverts the power equation", {
  skip_on_cran()

  cfg <- list(
    p = 0.5,
    sens = 0.8,
    spec = 0.8,
    N = 5000,
    delta = 0.04,
    alpha = 0.05,
    target_power = 0.8
  )

  n_req <- power_eif_binary(
    delta = cfg$delta,
    N = cfg$N,
    m_cal = NULL,
    power = cfg$target_power,
    p = cfg$p,
    sens = cfg$sens,
    spec = cfg$spec,
    alpha = cfg$alpha
  )

  achieved <- power_eif_binary(
    delta = cfg$delta,
    N = cfg$N,
    m_cal = n_req,
    p = cfg$p,
    sens = cfg$sens,
    spec = cfg$spec,
    alpha = cfg$alpha
  )

  expect_gt(achieved, cfg$target_power - 0.02)
})
