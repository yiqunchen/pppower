test_that("power_ppi_mean formula yields reasonable power range", {
  out <- power_ppi_mean(delta = 0.2, var_f = 0.3, var_res = 1.0, N = 1000, n = 200)
  expect_true(out > 0 && out < 1)
})

test_that("power_ppi_mean equals alpha when delta = 0 (two-sided test)", {
  out <- power_ppi_mean(delta = 0, var_f = 0.3, var_res = 1.0, N = 1000, n = 200, alpha = 0.05)
  expect_equal(out, 0.05, tolerance = 1e-3)
})

test_that("power_ppi_mean approaches 1 when delta is large", {
  out <- power_ppi_mean(delta = 5, var_f = 0.3, var_res = 1.0, N = 1000, n = 200)
  expect_gt(out, 0.999)
})

test_that("Monte Carlo estimate matches analytical power", {
  sim <- simulate_power(delta = 0.2, var_f = 0.4, var_res = 1.0, N = 1000, n = 200, R = 50000)
  expect_equal(unname(sim["Empirical"]), unname(sim["Exact"]), tolerance = 0.01)
})

test_that("Monte Carlo agreement holds across a range of parameters", {
  set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))

  grid <- expand.grid(
    delta   = c(0.1, 0.2, 0.5),
    var_f   = c(0.2, 0.4),
    var_res = c(0.2, 0.5, 0.8, 1.0),
    N       = c(500, 1000, 3000),
    n       = c(50, 100, 200),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    sim <- simulate_power(
      delta = g$delta, var_f = g$var_f, var_res = g$var_res,
      N = g$N, n = g$n, alpha = 0.05, R = 50000
    )
    empirical <- unname(sim["Empirical"])
    expected  <- unname(sim["Exact"])

    # Monte Carlo SE for a Bernoulli proportion with R draws
    R <- 5e4
    se_mc <- sqrt(expected * (1 - expected) / R)

    eps <- .Machine$double.eps^0.5
    ok <- (abs(empirical - expected) < 4 * se_mc) ||
        (abs(empirical - expected) < eps) ||
        (expected > 0.99 & empirical > 0.99)

    msg <- sprintf(
      "Empirical=%.4f, Exact=%.4f, |diff|=%.4f, 4*MCSE=%.4f, cfg: delta=%.2f var_f=%.2f var_res=%.2f N=%d n=%d",
      empirical, expected, abs(empirical - expected), 4 * se_mc,
      g$delta, g$var_f, g$var_res, g$N, g$n
    )
    expect_true(ok, info = msg)
  }
})