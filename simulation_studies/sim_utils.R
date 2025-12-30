#' Simulation Utilities for PPPower Validation
#'
#' Simple, focused functions for power analysis validation.
#' No external dependencies beyond base R.

#' Generate bivariate normal (Y, f) data
#'
#' @param n Sample size
#' @param mu_y True mean of Y (theta_star)
#' @param sigma_y SD of Y
#' @param sigma_f SD of f (predictor)
#' @param rho Correlation between Y and f
#' @return data.frame with columns y and f
generate_bivariate_normal <- function(n, mu_y = 0, sigma_y = 1, sigma_f = 1, rho = 0.5) {
  # Generate correlated normals using Cholesky
  z1 <- rnorm(n)
  z2 <- rnorm(n)

  y <- mu_y + sigma_y * z1
  f <- sigma_f * (rho * z1 + sqrt(1 - rho^2) * z2)

  data.frame(y = y, f = f)
}

#' Single Monte Carlo replicate for PPI mean estimation
#'
#' @param n Labeled sample size
#' @param N Unlabeled sample size
#' @param mu_y True mean (theta_star)
#' @param sigma_y SD of Y
#' @param sigma_f SD of f
#' @param rho Correlation between Y and f
#' @param alpha Significance level
#' @param use_ppi_plus Use PPI++ (TRUE) or vanilla PPI (FALSE)
#' @return List with theta_hat, se, z_stat, reject
ppi_mean_one_rep <- function(n, N, mu_y = 0, sigma_y = 1, sigma_f = 1, rho = 0.5,
                              alpha = 0.05, use_ppi_plus = FALSE) {
  # Generate labeled data
  dat_L <- generate_bivariate_normal(n, mu_y, sigma_y, sigma_f, rho)

  # Generate unlabeled data (only f observed)
  dat_U <- generate_bivariate_normal(N, mu_y, sigma_y, sigma_f, rho)

  y_L <- dat_L$y
  f_L <- dat_L$f
  f_U <- dat_U$f

  if (use_ppi_plus) {
    # PPI++ estimator with plug-in lambda
    cov_yf <- cov(y_L, f_L)
    var_f <- var(c(f_L, f_U))  # pooled variance
    r <- n / N

    # Optimal lambda (clipped to [0, 1])
    lambda_raw <- cov_yf / ((1 + r) * var_f)
    lambda <- pmin(1, pmax(0, lambda_raw))

    # PPI++ estimator
    theta_hat <- mean(y_L) + lambda * (mean(f_U) - mean(f_L))

    # Variance components
    var_res <- var(y_L - lambda * f_L)
    se <- sqrt(var_res / n + (lambda^2) * var(f_U) / N)

  } else {
    # Vanilla PPI estimator (lambda = 1)
    theta_hat <- mean(f_U) + mean(y_L - f_L)

    # Variance components
    var_f_U <- var(f_U)
    var_res <- var(y_L - f_L)
    se <- sqrt(var_f_U / N + var_res / n)
    lambda <- 1
  }

  # Test H0: theta = 0 vs H1: theta != 0
  z_stat <- theta_hat / se
  z_crit <- qnorm(1 - alpha / 2)
  reject <- abs(z_stat) > z_crit

  list(
    theta_hat = theta_hat,
    se = se,
    z_stat = z_stat,
    reject = reject,
    lambda = lambda
  )
}

#' Compute empirical power via Monte Carlo
#'
#' @param R Number of MC replicates
#' @param n Labeled sample size
#' @param N Unlabeled sample size
#' @param delta Effect size (mu_y under H1)
#' @param sigma_y SD of Y
#' @param sigma_f SD of f
#' @param rho Correlation between Y and f
#' @param alpha Significance level
#' @param use_ppi_plus Use PPI++ (TRUE) or vanilla PPI (FALSE)
#' @param seed Random seed for reproducibility
#' @return List with empirical_power, mc_se, analytical_power
compute_empirical_power <- function(R = 500, n, N, delta, sigma_y = 1, sigma_f = 1,
                                     rho = 0.5, alpha = 0.05, use_ppi_plus = FALSE,
                                     seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # Run MC replicates
  rejects <- replicate(R, {
    res <- ppi_mean_one_rep(n, N, mu_y = delta, sigma_y, sigma_f, rho, alpha, use_ppi_plus)
    res$reject
  })

  emp_power <- mean(rejects)
  mc_se <- sqrt(emp_power * (1 - emp_power) / R)

  # Analytical power
  cov_yf <- rho * sigma_y * sigma_f
  var_eps <- sigma_y^2 - 2 * cov_yf + sigma_f^2  # Var(Y - f)
  r <- n / N

  if (use_ppi_plus) {
    # PPI++ optimal lambda
    lambda_star <- cov_yf / ((1 + r) * sigma_f^2)
    # PPI++ variance
    var_theta <- sigma_y^2 / n - (cov_yf^2 / sigma_f^2) * N / (n * (n + N))
    se_analytical <- sqrt(var_theta)
  } else {
    # Vanilla PPI variance
    se_analytical <- sqrt(sigma_f^2 / N + var_eps / n)
  }

  z_alpha <- qnorm(1 - alpha / 2)
  analytical_power <- 1 - pnorm(z_alpha - delta / se_analytical) +
                      pnorm(-z_alpha - delta / se_analytical)

  list(
    empirical_power = emp_power,
    mc_se = mc_se,
    analytical_power = analytical_power,
    se_analytical = se_analytical,
    n = n,
    N = N,
    delta = delta,
    rho = rho,
    method = if (use_ppi_plus) "PPI++" else "PPI"
  )
}

#' Simple power curve comparison
#'
#' @param n_values Vector of labeled sample sizes
#' @param delta_values Vector of effect sizes
#' @param N Unlabeled sample size
#' @param sigma_y SD of Y
#' @param sigma_f SD of f
#' @param rho Correlation between Y and f
#' @param R MC replicates per configuration
#' @param use_ppi_plus Use PPI++ or vanilla PPI
#' @return data.frame with results
power_curve_grid <- function(n_values, delta_values, N, sigma_y = 1, sigma_f = 1,
                              rho = 0.5, R = 200, use_ppi_plus = FALSE, seed = 42) {
  results <- list()
  k <- 1

  for (n in n_values) {
    for (delta in delta_values) {
      res <- compute_empirical_power(
        R = R, n = n, N = N, delta = delta,
        sigma_y = sigma_y, sigma_f = sigma_f, rho = rho,
        use_ppi_plus = use_ppi_plus, seed = seed + k
      )
      results[[k]] <- as.data.frame(res)
      k <- k + 1
    }
  }

  do.call(rbind, results)
}
