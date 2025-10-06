simulate_power <- function(delta, var_f, var_res, N, n, alpha = 0.05, R = 100000) {
  se <- sqrt(var_f / N + var_res / n)
  mu <- abs(delta) / se
  Z <- rnorm(R, mean = mu, sd = 1)
  c_empirical = mean(abs(Z) > qnorm(1 - alpha / 2))
  c_exact     = power_ppi_mean(delta, var_f, var_res, N, n, alpha)
  c(Exact = c_exact, Empirical = c_empirical)
}

simulate_crossfit_data <- function(n = 2000, p = 5,
                                   family = stats::binomial(),
                                   K = 5, seed = 1) {
  set.seed(seed)
  X <- matrix(rnorm(n * p), n, p)
  colnames(X) <- paste0("x", 1:p)

  beta <- seq(1, length.out = p, by = 0.5)

  if (identical(family$family, "binomial")) {
    # Logistic DGP
    eta <- drop(X %*% beta)
    prob <- plogis(eta)
    y <- rbinom(n, size = 1, prob = prob)
  } else {
    # Gaussian DGP
    mu <- drop(X %*% beta)
    y <- mu + rnorm(n, sd = 1.0)
  }

  # Cross-fitted predictions (out-of-fold)
  fhat_cf <- crossfit_glm(X, y, K = K, family = family, seed = seed)

  # Return tidy data.frame
  df <- data.frame(y = y, X, fhat_cf = fhat_cf)
  attr(df, "family") <- family$family
  attr(df, "K") <- K
  df
}
