test_that("agent tools are defined and match the design functions", {
  skip_if_not_installed("ellmer")
  skip_if_not_installed("jsonlite")
  tools <- pppower_mcp_tools()
  nms <- vapply(tools, function(t) t@name, character(1))
  expect_setequal(nms, c("ppi_mean_design", "ppi_two_sample_design",
                         "ppi_two_proportion_design", "ppi_paired_design",
                         "ppi_2x2_design", "ppi_planning_inputs"))
  call_tool <- function(.tool, ...) {
    jsonlite::fromJSON(tools[[which(nms == .tool)]](...))
  }

  # one-sample mean: matches power_ppi_mean
  res <- call_tool("ppi_mean_design", delta = 0.2, N = 1000, r2 = 0.49)
  n_direct <- ceiling(suppressWarnings(power_ppi_mean(
    delta = 0.2, N = 1000, n = NULL, power = 0.8,
    sigma_y2 = 1, sigma_f2 = 1, cov_y_f = 0.7)))
  expect_equal(res$labeled_n_required, n_direct)
  expect_equal(res$classical_n_required, 197)

  # binary inputs
  res_b <- call_tool("ppi_mean_design", delta = 0.05, N = 1000,
                     prevalence = 0.3, sensitivity = 0.85, specificity = 0.85)
  expect_equal(res_b$labeled_n_required, 458)

  # power mode returns a probability
  res_p <- call_tool("ppi_mean_design", delta = 0.2, N = 1000, r2 = 0.49, n = 110)
  expect_gte(res_p$power_ppi, 0.8)
  expect_lt(res_p$power_classical, res_p$power_ppi)

  # 2x2 relative risk matches power_ppi_2x2
  res_rr <- call_tool("ppi_2x2_design", p_exposed = 0.4, p_control = 0.2, N = 1000,
                      sensitivity = 0.8, specificity = 0.8, effect_measure = "RR")
  expect_equal(res_rr$labeled_n_required_total,
               suppressWarnings(power_ppi_2x2(p_exp = 0.4, p_ctrl = 0.2, N = 1000,
                                              n = NULL, power = 0.8, sens = 0.8,
                                              spec = 0.8, effect_measure = "RR")))

  # planning inputs: conservative R^2 is below the point estimate
  inp <- call_tool("ppi_planning_inputs", pilot_r = 0.7, pilot_n = 60)
  expect_lt(inp$r2_lower_bound, inp$r2_point)

  # missing inputs give an informative error
  expect_error(tools[[which(nms == "ppi_mean_design")]](delta = 0.2, N = 1000), "r2")
})

test_that("ppi_design gives the same answers as the design functions", {
  res <- ppi_design("mean", delta = 0.2, N = 1000, r2 = 0.49)
  expect_equal(res$labeled_n_required, 110L)
  expect_equal(res$classical_n_required, 197L)
  expect_equal(ppi_design("paired", delta = 0.3, N = 500, sigma_d2 = 1, r2_d = 0.25)$labeled_pairs_required,
               as.integer(ceiling(suppressWarnings(power_ppi_paired(delta = 0.3, N = 500, n = NULL, power = 0.8,
                                                                    sigma_D2 = 1, rho_D = 0.5)))))
  expect_error(ppi_design("mean", delta = 0.2, N = 1000), "r2")
  expect_error(ppi_design("unknown"))
})
