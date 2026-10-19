test_that("plugin_correction lowers power and adds about one label for Gaussian data", {
  args <- list(delta = 0.2, N = 1000, sigma_y2 = 1, sigma_f2 = 1, cov_y_f = 0.7)
  p0 <- do.call(power_ppi_mean, c(args, n = 110))
  p1 <- do.call(power_ppi_mean, c(args, n = 110, plugin_correction = TRUE))
  expect_lt(p1, p0)
  n0 <- suppressWarnings(do.call(power_ppi_mean, c(args, list(n = NULL, power = 0.8))))
  n1 <- suppressWarnings(do.call(power_ppi_mean, c(args, list(n = NULL, power = 0.8, plugin_correction = TRUE))))
  expect_gte(n1, n0)
  expect_lte(n1 - n0, 2)
  expect_gte(do.call(power_ppi_mean, c(args, n = n1, plugin_correction = TRUE)), 0.8)
})

test_that("plugin variance term matches c * s_eps2 * kappa / (n (n - 3))", {
  n <- 30; N <- 500; rho <- 0.6
  v <- pppower:::plugin_variance_term(n, N, 1, 1, rho, kappa = 2)
  expect_equal(v, N / (n + N) * (1 - rho^2) * 2 / (n * (n - 3)))
  expect_equal(pppower:::plugin_variance_term(n, N, 1, 1, rho, kappa = 0.5),
               pppower:::plugin_variance_term(n, N, 1, 1, rho, kappa = 1))
})

test_that("corrected standard error is larger and uses a t reference", {
  set.seed(1)
  n <- 20; N <- 1000
  y <- rnorm(n); f <- 0.6 * y + rnorm(n, sd = 0.8); fu <- rnorm(N, sd = sqrt(0.36 + 0.64))
  a <- ppi_mean_test(y, f, fu)
  b <- ppi_mean_test(y, f, fu, se_type = "corrected")
  expect_equal(a$estimate, b$estimate)
  expect_gt(b$se, a$se)
  expect_equal(a$df, Inf)
  expect_equal(b$df, n - 2)
  expect_equal(b$p_value, 2 * pt(-abs(b$z_stat), n - 2))
  # a fixed lambda needs no correction
  c1 <- ppi_mean_test(y, f, fu, lambda = 0.5, se_type = "corrected")
  c0 <- ppi_mean_test(y, f, fu, lambda = 0.5)
  expect_equal(c1$se, c0$se)
})

test_that("plug-in lambda is 0 when all labeled predictions are equal", {
  set.seed(3)
  out <- ppi_mean_test(Y_L = rbinom(15, 1, 0.3), f_L = rep(1, 15), f_U = rbinom(500, 1, 0.3))
  expect_equal(out$lambda, 0)
  expect_true(is.finite(out$estimate))
})
