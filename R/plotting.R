
#' Power curve for the PP mean estimator
#'
#' @description
#' Computes Monte Carlo power values for a grid of labeled sample sizes using
#' `simulate_power()` given pre-computed variance components.
#'
#' @param n_grid Integer vector of candidate labeled sample sizes.
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param N Unlabeled sample size.
#' @param var_f Variance of \eqn{f(X)}.
#' @param var_res Variance of residuals \eqn{Y - f(X)}.
#' @param alpha Two-sided significance level.
#' @param R Number of Monte Carlo replicates passed to `simulate_power()`.
#' @param seed Optional RNG seed for reproducibility.
#'
#' @return Data frame with columns `n`, `power_empirical`, `power_exact`, `delta`,
#'   `N`, `alpha`, `var_f`, and `var_res` ready for plotting power curves.
#'
#' @examples
#' power_curve_mean(
#'   n_grid = seq(50, 200, by = 25),
#'   delta = 0.2,
#'   N = 2000,
#'   var_f = 0.4,
#'   var_res = 1.0,
#'   R = 5000,
#'   seed = 123
#' )
#'
#' @export
power_curve_mean <- function(n_grid,
                             delta,
                             N,
                             var_f,
                             var_res,
                             alpha = 0.05,
                             R = 10000,
                             seed = NULL) {
  if (!is.numeric(n_grid) || length(n_grid) == 0L) {
    stop("n_grid must be a numeric vector of candidate sample sizes.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta)) {
    stop("delta must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || N <= 0) {
    stop("N must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(var_f) || length(var_f) != 1L || var_f < 0) {
    stop("var_f must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(var_res) || length(var_res) != 1L || var_res <= 0) {
    stop("var_res must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || R <= 0) {
    stop("R must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be a finite numeric scalar when provided.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  n_grid <- as.integer(sort(unique(n_grid)))
  if (any(is.na(n_grid)) || any(n_grid <= 0)) {
    stop("n_grid must contain positive integers only.", call. = FALSE)
  }

  sims <- lapply(
    n_grid,
    function(n) {
      simulate_power(
        delta = delta,
        var_f = var_f,
        var_res = var_res,
        N = N,
        n = n,
        alpha = alpha,
        R = as.integer(R)
      )
    }
  )
  powers_empirical <- vapply(sims, function(x) unname(x["Empirical"]), numeric(1L))
  powers_exact <- vapply(sims, function(x) unname(x["Exact"]), numeric(1L))

  data.frame(
    n = n_grid,
    power_empirical = powers_empirical,
    power_exact = powers_exact,
    delta = rep(delta, length(n_grid)),
    N = rep(N, length(n_grid)),
    alpha = rep(alpha, length(n_grid)),
    var_f = rep(var_f, length(n_grid)),
    var_res = rep(var_res, length(n_grid)),
    row.names = NULL
  )
}

#' Power curve using Gaussian/Binomial DGPs
#'
#' @description
#' Generates a synthetic superpopulation via `simulate_crossfit_data()` and
#' computes Monte Carlo power curves for the PP mean estimator over a grid of
#' labeled sample sizes using `simulate_power_ppi_mean()`.
#'
#' @inheritParams power_curve_mean
#' @param theta0 Null value for the mean estimand.
#' @param family GLM family passed to the DGP; defaults to `stats::gaussian()`
#'   (with identity link) or `stats::binomial()` (logistic).
#' @param superpop_n Size of the synthetic superpopulation used to estimate
#'   variance components.
#' @param p Number of covariates in the DGP.
#' @param K Number of folds for cross-fitting.
#' @param seed Random seed forwarded to `simulate_crossfit_data()`.
#' @param R Number of Monte Carlo replicates passed to `simulate_power_ppi_mean()`.
#'
#' @return Data frame with the columns produced by `power_curve_mean()` plus
#'   `family`, `theta`, and `theta0` for downstream plotting.
#'
#' @examples
#' power_curve_mean_dgp(
#'   n_grid = seq(50, 200, by = 25),
#'   N = 5000,
#'   theta0 = 0,
#'   family = stats::gaussian(),
#'   R = 1000
#' )
#'
#' @export
power_curve_mean_dgp <- function(n_grid,
                                 N,
                                 theta0 = 0,
                                 family = stats::gaussian(),
                                 superpop_n = 5000,
                                 p = 5,
                                 K = 5,
                                 seed = 1,
                                 alpha = 0.05,
                                 R = 2000) {
  if (!inherits(family, "family")) {
    stop("family must be a valid stats::family object.", call. = FALSE)
  }
  if (!is.numeric(theta0) || length(theta0) != 1L || !is.finite(theta0)) {
    stop("theta0 must be a finite numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(superpop_n) || length(superpop_n) != 1L || superpop_n <= 0) {
    stop("superpop_n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || p <= 0) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || K <= 1) {
    stop("K must be an integer greater than 1.", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || R <= 0) {
    stop("R must be a positive numeric scalar.", call. = FALSE)
  }

  n_grid <- as.integer(sort(unique(n_grid)))
  if (any(is.na(n_grid)) || any(n_grid <= 0)) {
    stop("n_grid must contain positive integers only.", call. = FALSE)
  }

  df <- simulate_crossfit_data(
    n = as.integer(superpop_n),
    p = as.integer(p),
    family = family,
    K = as.integer(K),
    seed = seed
  )

  theta <- mean(df$y)
  delta <- theta - theta0
  var_f <- stats::var(df$fhat_cf)
  var_res <- stats::var(df$y - df$fhat_cf)

  seeds <- seed + seq_along(n_grid)
  sims <- lapply(
    seq_along(n_grid),
    function(i) {
      simulate_power_ppi_mean(
        df = df,
        N = N,
        n = n_grid[i],
        alpha = alpha,
        R = as.integer(R),
        theta0 = theta0,
        seed = seeds[i]
      )
    }
  )

  power_empirical <- vapply(sims, function(x) x$empirical_power, numeric(1L))
  power_exact <- vapply(sims, function(x) x$analytical_power, numeric(1L))

  data.frame(
    n = n_grid,
    power_empirical = power_empirical,
    power_exact = power_exact,
    delta = rep(delta, length(n_grid)),
    N = rep(N, length(n_grid)),
    alpha = rep(alpha, length(n_grid)),
    var_f = rep(var_f, length(n_grid)),
    var_res = rep(var_res, length(n_grid)),
    family = rep(family$family, length(n_grid)),
    theta = rep(theta, length(n_grid)),
    theta0 = rep(theta0, length(n_grid)),
    row.names = NULL
  )
}

#' Convenience wrapper for Gaussian DGP power curves
#'
#' @inherit power_curve_mean_dgp params return
#'
#' @export
power_curve_gaussian <- function(n_grid,
                                 N,
                                 theta0 = 0,
                                 superpop_n = 5000,
                                 p = 5,
                                 K = 5,
                                 seed = 1,
                                 alpha = 0.05,
                                 R = 2000) {
  power_curve_mean_dgp(
    n_grid = n_grid,
    N = N,
    theta0 = theta0,
    family = stats::gaussian(),
    superpop_n = superpop_n,
    p = p,
    K = K,
    seed = seed,
    alpha = alpha,
    R = R
  )
}

#' Convenience wrapper for Binomial DGP power curves
#'
#' @inherit power_curve_mean_dgp params return
#'
#' @export
power_curve_binomial <- function(n_grid,
                                 N,
                                 theta0 = 0,
                                 superpop_n = 5000,
                                 p = 5,
                                 K = 5,
                                 seed = 1,
                                 alpha = 0.05,
                                 R = 2000) {
  power_curve_mean_dgp(
    n_grid = n_grid,
    N = N,
    theta0 = theta0,
    family = stats::binomial(),
    superpop_n = superpop_n,
    p = p,
    K = K,
    seed = seed,
    alpha = alpha,
    R = R
  )
}


#' Type I error curve for the PP mean estimator
#'
#' @description
#' Computes empirical and analytical Type I error estimates across a grid of
#' effect sizes using Monte Carlo via `simulate_power()`. When the null is true
#' (`effect_size = 0`), the curve reports the Type I error; for other effect
#' sizes the values coincide with the rejection probability (i.e., power).
#'
#' @param effect_grid Numeric vector of effect sizes \eqn{\theta - \theta_0} to
#'   evaluate.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param var_f Variance of \eqn{f(X)}.
#' @param var_res Variance of residuals \eqn{Y - f(X)}.
#' @param alpha Two-sided significance level.
#' @param R Number of Monte Carlo replicates passed to `simulate_power()`.
#' @param seed Optional RNG seed for reproducibility.
#'
#' @return Data frame with columns `effect_size`, `type1_empirical`,
#'   `type1_exact`, `N`, `n`, `alpha`, `var_f`, and `var_res`.
#'
#' @examples
#' type1_error_curve_mean(
#'   effect_grid = seq(-0.4, 0.4, by = 0.05),
#'   N = 4000,
#'   n = 200,
#'   var_f = 0.4,
#'   var_res = 1.1,
#'   R = 2000
#' )
#'
#' @export
type1_error_curve_mean <- function(effect_grid,
                                   N,
                                   n,
                                   var_f,
                                   var_res,
                                   alpha = 0.05,
                                   R = 10000,
                                   seed = NULL) {
  if (!is.numeric(effect_grid) || length(effect_grid) == 0L) {
    stop("effect_grid must be a numeric vector of effect sizes.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || N <= 0) {
    stop("N must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("n must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(var_f) || length(var_f) != 1L || var_f < 0) {
    stop("var_f must be a non-negative numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(var_res) || length(var_res) != 1L || var_res <= 0) {
    stop("var_res must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || R <= 0) {
    stop("R must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be a finite numeric scalar when provided.", call. = FALSE)
    }
    set.seed(as.integer(seed))
  }

  effect_grid <- sort(effect_grid)
  sims <- lapply(
    effect_grid,
    function(delta) {
      simulate_power(
        delta = delta,
        var_f = var_f,
        var_res = var_res,
        N = N,
        n = n,
        alpha = alpha,
        R = as.integer(R)
      )
    }
  )

  type1_empirical <- vapply(sims, function(x) unname(x["Empirical"]), numeric(1L))
  type1_exact <- vapply(sims, function(x) unname(x["Exact"]), numeric(1L))

  data.frame(
    effect_size = effect_grid,
    type1_empirical = type1_empirical,
    type1_exact = type1_exact,
    N = rep(N, length(effect_grid)),
    n = rep(as.integer(n), length(effect_grid)),
    alpha = rep(alpha, length(effect_grid)),
    var_f = rep(var_f, length(effect_grid)),
    var_res = rep(var_res, length(effect_grid)),
    row.names = NULL
  )
}

#' Type I error curve using Gaussian/Binomial DGPs
#'
#' @description
#' Simulates superpopulation data via `simulate_crossfit_data()` and evaluates
#' empirical and analytical Type I error across an effect-size grid by repeatedly
#' calling `simulate_power_ppi_mean()`.
#'
#' @inheritParams power_curve_mean_dgp
#' @param effect_grid Numeric vector of effect sizes \eqn{\theta - \theta_0} to
#'   evaluate.
#' @param n Labeled sample size used inside `simulate_power_ppi_mean()`.
#'
#' @return Data frame with columns `effect_size`, `type1_empirical`,
#'   `type1_exact`, `theta`, `theta0`, and additional metadata mirroring the
#'   inputs.
#'
#' @examples
#' type1_error_curve_mean_dgp(
#'   effect_grid = seq(-0.4, 0.4, by = 0.05),
#'   N = 4000,
#'   n = 200,
#'   family = stats::gaussian(),
#'   R = 2000
#' )
#'
#' @export
type1_error_curve_mean_dgp <- function(effect_grid,
                                       N,
                                       n,
                                       family = stats::gaussian(),
                                       superpop_n = 5000,
                                       p = 5,
                                       K = 5,
                                       seed = 1,
                                       alpha = 0.05,
                                       R = 2000) {
  if (!is.numeric(effect_grid) || length(effect_grid) == 0L) {
    stop("effect_grid must be a numeric vector of effect sizes.", call. = FALSE)
  }
  if (!inherits(family, "family")) {
    stop("family must be a valid stats::family object.", call. = FALSE)
  }
  if (!is.numeric(N) || length(N) != 1L || N <= 0) {
    stop("N must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    stop("n must be a positive numeric scalar.", call. = FALSE)
  }
  if (!is.numeric(superpop_n) || length(superpop_n) != 1L || superpop_n <= 0) {
    stop("superpop_n must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(p) || length(p) != 1L || p <= 0) {
    stop("p must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(K) || length(K) != 1L || K <= 1) {
    stop("K must be an integer greater than 1.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    stop("alpha must lie in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(R) || length(R) != 1L || R <= 0) {
    stop("R must be a positive numeric scalar.", call. = FALSE)
  }

  effect_grid <- sort(effect_grid)
  df <- simulate_crossfit_data(
    n = as.integer(superpop_n),
    p = as.integer(p),
    family = family,
    K = as.integer(K),
    seed = seed
  )

  theta <- mean(df$y)
  seeds <- seed + seq_along(effect_grid)
  sims <- lapply(
    seq_along(effect_grid),
    function(i) {
      theta0_i <- theta - effect_grid[i]
      simulate_power_ppi_mean(
        df = df,
        N = N,
        n = n,
        alpha = alpha,
        R = as.integer(R),
        theta0 = theta0_i,
        seed = seeds[i]
      )
    }
  )

  type1_empirical <- vapply(sims, function(x) x$empirical_power, numeric(1L))
  type1_exact <- vapply(sims, function(x) x$analytical_power, numeric(1L))
  theta0_vals <- vapply(effect_grid, function(delta) theta - delta, numeric(1L))
  var_f <- stats::var(df$fhat_cf)
  var_res <- stats::var(df$y - df$fhat_cf)

  data.frame(
    effect_size = effect_grid,
    type1_empirical = type1_empirical,
    type1_exact = type1_exact,
    theta = rep(theta, length(effect_grid)),
    theta0 = theta0_vals,
    N = rep(N, length(effect_grid)),
    n = rep(as.integer(n), length(effect_grid)),
    alpha = rep(alpha, length(effect_grid)),
    var_f = rep(var_f, length(effect_grid)),
    var_res = rep(var_res, length(effect_grid)),
    family = rep(family$family, length(effect_grid)),
    row.names = NULL
  )
}

#' Plot Type I error curves
#'
#' @description
#' Creates a simple line plot of empirical and/or analytical Type I error (or
#' rejection probabilities) against effect size.
#'
#' @param curve_df Data frame returned by `type1_error_curve_mean()` or
#'   `type1_error_curve_mean_dgp()`.
#' @param empirical Logical; include the Monte Carlo estimate (`type1_empirical`).
#' @param exact Logical; include the analytical estimate (`type1_exact`).
#' @param add_reference Logical; add a horizontal line at the nominal level
#'   `alpha` when available.
#' @param empirical_col,exact_col Column names to use for empirical and exact
#'   curves. Override only if you have renamed the defaults.
#' @param legend_pos Character or numeric legend position passed to
#'   `graphics::legend()`, ignored when fewer than two curves are drawn.
#' @param ... Additional arguments forwarded to the initial `graphics::plot()`
#'   call (e.g., `main`, `xlab`, `ylab`).
#'
#' @return The input data frame, invisibly.
#'
#' @export
plot_type1_error_curve <- function(curve_df,
                                   empirical = TRUE,
                                   exact = TRUE,
                                   add_reference = TRUE,
                                   empirical_col = "type1_empirical",
                                   exact_col = "type1_exact",
                                   legend_pos = "topright",
                                   ...) {
  if (!empirical && !exact) {
    stop("At least one of empirical or exact must be TRUE.", call. = FALSE)
  }
  if (!is.data.frame(curve_df)) {
    stop("curve_df must be a data frame.", call. = FALSE)
  }
  if (!("effect_size" %in% names(curve_df))) {
    stop("curve_df must contain an 'effect_size' column.", call. = FALSE)
  }
  x <- curve_df$effect_size
  if (!is.numeric(x)) {
    stop("'effect_size' column must be numeric.", call. = FALSE)
  }

  y_emp <- if (empirical) curve_df[[empirical_col]] else NULL
  y_exact <- if (exact) curve_df[[exact_col]] else NULL
  xlab_default <- "Effect size (theta - theta0)"
  ylab_default <- "Type I error / rejection probability"

  if (empirical && !is.null(y_emp)) {
    graphics::plot(
      x,
      y_emp,
      type = "l",
      col = "steelblue",
      lwd = 2,
      xlab = xlab_default,
      ylab = ylab_default,
      ...
    )
  } else if (exact && !is.null(y_exact)) {
    graphics::plot(
      x,
      y_exact,
      type = "l",
      col = "firebrick",
      lwd = 2,
      lty = 2,
      xlab = xlab_default,
      ylab = ylab_default,
      ...
    )
  }

  if (empirical && exact && !is.null(y_emp) && !is.null(y_exact)) {
    graphics::lines(x, y_emp, col = "steelblue", lwd = 2)
    graphics::lines(x, y_exact, col = "firebrick", lwd = 2, lty = 2)
    graphics::legend(
      legend_pos,
      legend = c("Empirical", "Analytical"),
      col = c("steelblue", "firebrick"),
      lty = c(1, 2),
      lwd = 2,
      bty = "n"
    )
  } else {
    if (empirical && !is.null(y_emp)) {
      graphics::lines(x, y_emp, col = "steelblue", lwd = 2)
    }
    if (exact && !is.null(y_exact)) {
      graphics::lines(x, y_exact, col = "firebrick", lwd = 2, lty = 2)
    }
  }

  if (add_reference && "alpha" %in% names(curve_df)) {
    graphics::abline(h = unique(curve_df$alpha), col = "gray60", lty = 3)
  }

  invisible(curve_df)
}

