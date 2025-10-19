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

    expect_equal(
      power_ppi_mean(
        delta = delta,
        N = N,
        n = n,
        alpha = alpha,
        metrics = metrics,
        metric_type = metric_type
      ),
      power_ppi_mean(
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

    expect_equal(
      power_ppi_pp_mean(
        delta = delta,
        N = N,
        n = n,
        metrics = metrics,
        metric_type = metric_type
      ),
      power_ppi_pp_mean(
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
