test_that("rho_lower_bound matches the Fisher z formula and keeps the sign", {
  expect_equal(rho_lower_bound(0.7, 50), tanh(atanh(0.7) - qnorm(0.9) / sqrt(47)))
  expect_equal(rho_lower_bound(-0.7, 50), -rho_lower_bound(0.7, 50))
  expect_equal(rho_lower_bound(0.1, 10), 0)
  expect_lt(rho_lower_bound(0.7, 50, confidence = 0.95), rho_lower_bound(0.7, 50))
  expect_error(rho_lower_bound(0.7, 3), "pilot_n")
  expect_error(rho_lower_bound(0.7, 50, confidence = 1), "confidence")
})

test_that("pilot_n in power_ppi_mean plans at the lower bound", {
  rho_l <- rho_lower_bound(0.7, 50)
  n_pilot <- power_ppi_mean(delta = 0.2, N = 1000, power = 0.8,
                            sigma_y2 = 1, sigma_f2 = 1, cov_y_f = 0.7, pilot_n = 50)
  n_bound <- power_ppi_mean(delta = 0.2, N = 1000, power = 0.8,
                            sigma_y2 = 1, sigma_f2 = 1, cov_y_f = rho_l)
  n_point <- power_ppi_mean(delta = 0.2, N = 1000, power = 0.8,
                            sigma_y2 = 1, sigma_f2 = 1, cov_y_f = 0.7)
  expect_equal(n_pilot, n_bound)
  expect_gt(n_pilot, n_point)
  p_pilot <- power_ppi_mean(delta = 0.2, N = 1000, n = 120,
                            sigma_y2 = 2, sigma_f2 = 0.5, cov_y_f = 0.7, pilot_n = 50)
  p_bound <- power_ppi_mean(delta = 0.2, N = 1000, n = 120, sigma_y2 = 2, sigma_f2 = 0.5,
                            cov_y_f = rho_lower_bound(0.7, 50) * sqrt(2 * 0.5))
  expect_equal(p_pilot, p_bound)
  # combines with the plug-in correction
  n_both <- power_ppi_mean(delta = 0.2, N = 1000, power = 0.8, sigma_y2 = 1, sigma_f2 = 1,
                           cov_y_f = 0.7, pilot_n = 50, plugin_correction = TRUE)
  expect_gte(n_both, n_pilot)
})

test_that("pilot_n is passed through by the paired and two-sample functions", {
  rho_l <- rho_lower_bound(0.6, 40)
  expect_equal(
    power_ppi_paired(delta = 0.3, N = 500, power = 0.8, sigma_D2 = 1, rho_D = 0.6, pilot_n = 40),
    power_ppi_paired(delta = 0.3, N = 500, power = 0.8, sigma_D2 = 1, rho_D = rho_l))
  expect_equal(
    power_ppi_ttest(delta = 0.3, n_A = 80, n_B = 80, N_A = 500, N_B = 500,
                    sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = 0.6,
                    sigma_y2_B = 1, sigma_f2_B = 1, cov_yf_B = 0.6,
                    pilot_n_A = 40, pilot_n_B = 40),
    power_ppi_ttest(delta = 0.3, n_A = 80, n_B = 80, N_A = 500, N_B = 500,
                    sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = rho_l,
                    sigma_y2_B = 1, sigma_f2_B = 1, cov_yf_B = rho_l))
})

test_that("ppi_design accepts pilot_n", {
  a <- ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49, pilot_n = 50)
  b <- ppi_design("mean", delta = 0.2, N = 1000, r2 = rho_lower_bound(0.7, 50)^2)
  expect_equal(a$labeled_n_required, b$labeled_n_required)
  expect_equal(ppi_design("inputs", pilot_r = 0.7, pilot_n = 50)$r2_lower_bound,
               rho_lower_bound(0.7, 50)^2)
})
