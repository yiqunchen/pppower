test_that("power_ppi_vanilla_mean formula yields reasonable power range", {
  out <- power_ppi_vanilla_mean(delta = 0.2, var_f = 0.3, var_res = 1.0, N = 1000, n = 200)
  expect_true(out > 0 && out < 1)
})

test_that("power_ppi_vanilla_mean equals alpha when delta = 0 (two-sided test)", {
  out <- power_ppi_vanilla_mean(delta = 0, var_f = 0.3, var_res = 1.0, N = 1000, n = 200, alpha = 0.05)
  expect_equal(out, 0.05, tolerance = 1e-3)
})

test_that("power_ppi_vanilla_mean approaches 1 when delta is large", {
  out <- power_ppi_vanilla_mean(delta = 5, var_f = 0.3, var_res = 1.0, N = 1000, n = 200)
  expect_gt(out, 0.999)
})

test_that("Monte Carlo power aligns with analytical power for PPI and PPI++ mean", {
  R_sim <- 2000
  N <- 1500
  n <- 200
  alpha <- 0.05

  delta_pop <- 0.2
  var_f_pop <- 0.45
  var_res_pop <- 0.8

  z_alpha <- stats::qnorm(1 - alpha / 2)

  se_ppi <- sqrt(var_f_pop / N + var_res_pop / n)
  mu_ppi <- abs(delta_pop) / se_ppi
  analytical_ppi <- 1 - stats::pnorm(z_alpha - mu_ppi) + stats::pnorm(-z_alpha - mu_ppi)

  u <- (seq_len(R_sim) - 0.5) / R_sim
  theta_ppi <- delta_pop + se_ppi * stats::qnorm(u)
  z_ppi <- theta_ppi / se_ppi
  empirical_ppi <- mean(abs(z_ppi) > z_alpha)
  diff_ppi <- abs(empirical_ppi - analytical_ppi)
  mc_se_ppi <- sqrt(analytical_ppi * (1 - analytical_ppi) / R_sim)
  tol_ppi <- max(3 * mc_se_ppi, 0.01)

  expect_true(
    diff_ppi < tol_ppi,
    info = sprintf(
      "PPI-mean: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      empirical_ppi, analytical_ppi, diff_ppi, tol_ppi
    )
  )

  sigma_y2_pop <- var_f_pop + var_res_pop
  sigma_f2_pop <- var_f_pop
  cov_yf_pop <- var_f_pop
  lambda_pop <- ppi_pp_lambda_star(n = n, N = N, cov_y_f = cov_yf_pop, sigma_f2 = sigma_f2_pop)

  var_ppiplus <- ppi_pp_variance(
    n = n,
    N = N,
    sigma_y2 = sigma_y2_pop,
    sigma_f2 = sigma_f2_pop,
    cov_y_f = cov_yf_pop,
    lambda = lambda_pop,
    lambda_type = "user"
  )
  se_ppiplus <- sqrt(var_ppiplus)
  mu_ppiplus <- abs(delta_pop) / se_ppiplus
  analytical_ppiplus <- 1 - stats::pnorm(z_alpha - mu_ppiplus) + stats::pnorm(-z_alpha - mu_ppiplus)

  theta_ppiplus <- delta_pop + se_ppiplus * stats::qnorm(u)
  z_ppiplus <- theta_ppiplus / se_ppiplus
  empirical_ppiplus <- mean(abs(z_ppiplus) > z_alpha)
  diff_ppiplus <- abs(empirical_ppiplus - analytical_ppiplus)
  mc_se_ppiplus <- sqrt(analytical_ppiplus * (1 - analytical_ppiplus) / R_sim)
  tol_ppiplus <- max(3 * mc_se_ppiplus, 0.01)

  expect_true(
    diff_ppiplus < tol_ppiplus,
    info = sprintf(
      "PPI++-mean: Empirical=%.3f, Analytical=%.3f, |diff|=%.3f, tol=%.3f",
      empirical_ppiplus, analytical_ppiplus, diff_ppiplus, tol_ppiplus
    )
  )
})

test_that("power resolvers handle metrics input across metric types", {
  delta <- 0.28
  N <- 1200
  n <- 150
  alpha <- 0.05
  precision <- 0.78
  recall <- 0.60
  p_y <- 0.30

  scenarios <- list(
    list(
      metric_type = "continuous",
      metrics = list(
        type = "continuous",
        mse = 0.55,
        var_y = 1.0,
        cov_y_f = 0.30,
        bias = 0,
        m_obs = n
      )
    ),
    list(
      metric_type = "classification",
      metrics = list(
          precision = precision,
          recall = recall,
          p_y = p_y,
          m_obs = n
        )
    ),
    list(
      metric_type = "classification",
      metrics = list(
        tp = 48,
        fp = 12,
        fn = 10,
        m_obs = n,
        p_y = p_y
      )
    )
  )

  for (sc in scenarios) {
    metric_type <- sc$metric_type
    metrics <- sc$metrics
    info_msg <- paste("metric_type =", metric_type)

    resolved_pp <- resolve_ppi_variances(
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = metrics$m_obs,
      correction = TRUE
    )

    moments_ppplus <- resolve_ppi_pp_moments(
      var_f = resolved_pp$var_f,
      var_res = resolved_pp$var_res,
      metrics = metrics,
      metric_type = metric_type,
      m_labeled = metrics$m_obs,
      correction = TRUE
    )

    # Vanilla PPI via internal function
    expect_equal(
      power_ppi_vanilla_mean(
        delta = delta,
        N = N,
        n = n,
        alpha = alpha,
        metrics = metrics,
        metric_type = metric_type
      ),
      power_ppi_vanilla_mean(
        delta = delta,
        var_f = resolved_pp$var_f,
        var_res = resolved_pp$var_res,
        N = N,
        n = n,
        alpha = alpha
      ),
      tolerance = 1e-10,
      info = info_msg
    )

    # PPI++ via new unified power_ppi_mean
    expect_equal(
      power_ppi_mean(
        delta = delta,
        N = N,
        n = n,
        sigma_y2 = moments_ppplus$sigma_y2,
        sigma_f2 = moments_ppplus$sigma_f2,
        cov_y_f = moments_ppplus$cov_y_f,
        var_f = moments_ppplus$var_f,
        var_res = moments_ppplus$var_res
      ),
      power_ppi_mean(
        delta = delta,
        N = N,
        n = n,
        sigma_y2 = moments_ppplus$sigma_y2,
        sigma_f2 = moments_ppplus$sigma_f2,
        cov_y_f = moments_ppplus$cov_y_f,
        var_f = moments_ppplus$var_f,
        var_res = moments_ppplus$var_res
      ),
      tolerance = 1e-10,
      info = info_msg
    )
  }
})


test_that("Monte Carlo agreement holds across a range of parameters", {
  set.seed(as.integer(format(Sys.Date(), "%Y%m%d")))
  n_sim <- 50000
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
    var_y  <- g$var_f + g$var_res
    cov_yf <- g$var_f

    res <- simulate_power(
      delta   = g$delta,
      var_f   = g$var_f,
      var_res = g$var_res,
      cov_y_f = cov_yf,
      metrics = list(
        type  = "continuous",
        mse   = g$var_res,
        var_y = var_y,
        m_obs = g$n
      ),
      metric_type = "continuous",
      N       = g$N,
      n       = g$n,
      alpha   = 0.05,
      R       = n_sim
    )

    empirical_pp <- unname(res["Empirical_PP"])
    exact_pp  <- unname(res["Exact_PP"])

    se_mc <- sqrt(exact_pp * (1 - exact_pp) / n_sim)
    eps <- .Machine$double.eps^0.5
    ok <- (abs(empirical_pp - exact_pp) < 4 * se_mc) ||
      (abs(empirical_pp - exact_pp) < eps) ||
      (exact_pp > 0.99 & empirical_pp > 0.99)

    msg <- sprintf(
      "Empirical_PP=%.4f, Exact_PP=%.4f, |diff|=%.4f, 4*MCSE=%.4f, cfg: delta=%.2f var_f=%.2f var_res=%.2f N=%d n=%d",
      empirical_pp, exact_pp, abs(empirical_pp - exact_pp), 4 * se_mc,
      g$delta, g$var_f, g$var_res, g$N, g$n
    )
    expect_true(ok, info = msg)

    if ("Exact_PPplus" %in% names(res)) {
      empirical_ppplus <- unname(res["Empirical_PPplus"])
      exact_ppplus <- unname(res["Exact_PPplus"])
      rel_err <- abs(empirical_ppplus - exact_ppplus) / max(abs(exact_ppplus), eps)
      expect_lt(
        rel_err,
        0.05,
        sprintf(
          "Empirical_PPplus=%.4f, Exact_PPplus=%.4f, rel.err=%.4f, cfg: delta=%.2f var_f=%.2f var_res=%.2f N=%d n=%d",
          empirical_ppplus, exact_ppplus, rel_err,
          g$delta, g$var_f, g$var_res, g$N, g$n
        )
      )
    }
  }
})
