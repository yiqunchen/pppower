
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
  powers_empirical <- vapply(sims, function(x) unname(x["Empirical_PP"]), numeric(1L))
  powers_exact <- vapply(sims, function(x) unname(x["Exact_PP"]), numeric(1L))

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

#' Power curve for PPI mean with supplied population moments
#'
#' @description
#' Computes **analytical** and **Monte Carlo** power for the prediction-powered
#' mean estimator across a grid of labeled sample sizes `n_grid`, given
#' population moments and an effect size. This function does **not** simulate a
#' superpopulation; instead, you provide the required moments:
#' \eqn{\sigma_f^2 = Var(f(X))}, \eqn{\sigma_{\mathrm{res}}^2 = Var(Y - f(X))},
#' and the effect size \eqn{\delta = \theta - \theta_0}. For each `n` in
#' `n_grid`, analytical power is computed from the normal theory formula using
#' \eqn{\sqrt{\sigma_f^2/N + \sigma_{\mathrm{res}}^2/n}}, and empirical power is
#' estimated via Monte Carlo using `simulate_ppi_vanilla_mean()`.
#'
#' @param n_grid Numeric vector of labeled sample sizes to evaluate.
#' @param N Unlabeled sample size.
#' @param theta0 Null value for the mean estimand \eqn{\theta_0}.
#' @param var_f Population variance of predictions \eqn{Var(f(X))}.
#' @param var_res Population residual variance \eqn{Var(Y - f(X))}.
#' @param delta Effect size \eqn{\theta - \theta_0}.
#' @param family GLM family (default `stats::gaussian()`). Included for API
#'   consistency; the analytical calculation depends only on `var_f`, `var_res`,
#'   `N`, `n`, and `delta`.
#' @param R Number of Monte Carlo replicates for empirical power (default 1000).
#' @param alpha Two-sided test level (default 0.05).
#' @param seed RNG seed used inside the per-`n` simulations.
#'
#' @return A data frame with one row per `n` in `n_grid`, containing:
#'   \describe{
#'     \item{n}{Labeled sample size.}
#'     \item{analytical}{Analytical power from the normal theory formula.}
#'     \item{empirical}{Monte Carlo power (rejection rate).}
#'     \item{abs_diff}{Absolute difference |empirical - analytical|.}
#'   }
#'
#' @examples
#' power_curve_mean_dgp(
#'   n_grid = seq(50, 200, by = 25),
#'   N = 3000,
#'   theta0 = 0,
#'   var_f = 0.45,
#'   var_res = 0.80,
#'   delta = 0.2,
#'   R = 1000
#' )
#' @export
power_curve_mean_dgp <- function(
  n_grid = seq(50, 200, by = 25),
  N,
  theta0 = 0,
  var_f = 0.45,
  var_res = 0.80,
  delta = 0.2,
  family = stats::gaussian(),
  R = 1000,
  alpha = 0.05,
  seed = 1
) {
  if (!is.numeric(n_grid) || length(n_grid) == 0L) {
    stop("n_grid must be a numeric vector of sample sizes.", call. = FALSE)
  }

  results <- lapply(
    n_grid,
    function(n_i) {
      moments_ppi <- list(
        delta   = delta,
        var_f   = var_f,
        var_res = var_res
      )

      res <- simulate_ppi_vanilla_mean(
        R       = R,
        n       = n_i,
        N       = N,
        alpha   = alpha,
        family  = family,
        moments = moments_ppi,
        theta0  = theta0,
        seed    = seed + n_i
      )

      data.frame(
        n             = n_i,
        analytical    = res$analytical_power,
        empirical     = res$empirical_power,
        abs_diff      = abs(res$empirical_power - res$analytical_power),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(rbind, results)
}

# Thin wrapper functions removed - use power_curve_mean_dgp() directly with family argument


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

  type1_empirical <- vapply(sims, function(x) unname(x["Empirical_PP"]), numeric(1L))
  type1_exact <- vapply(sims, function(x) unname(x["Exact_PP"]), numeric(1L))

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

#' Type I error and power curve for the PPI mean estimator
#'
#' @description
#' Computes empirical and analytical Type I error (or power) across a grid of
#' effect sizes for the prediction-powered mean estimator, **without simulating
#' a superpopulation**. Instead of generating data, you supply the key population
#' quantities directly:
#' \eqn{\sigma_f^2 = Var(f(X))} and \eqn{\sigma_{\mathrm{res}}^2 = Var(Y - f(X))}.
#' Analytical power is computed from the normal theory formula
#' \eqn{\sqrt{\sigma_f^2/N + \sigma_{\mathrm{res}}^2/n}}, while empirical power is
#' estimated via repeated Monte Carlo simulations using
#' `simulate_ppi_vanilla_mean()`.
#'
#' @param effect_grid Numeric vector of effect sizes \eqn{\delta = \theta - \theta_0}
#'   over which to evaluate Type I error / power.
#' @param N Unlabeled sample size.
#' @param n Labeled sample size.
#' @param var_f Variance of predictions \eqn{Var(f(X))}.
#' @param var_res Residual variance \eqn{Var(Y - f(X))}.
#' @param theta True mean of \eqn{Y} under the alternative. Defaults to 0.
#'   Only affects \eqn{\theta_0 = \theta - \delta}; it does not change power.
#' @param alpha Two-sided significance level (default 0.05).
#' @param R Number of Monte Carlo replicates (default 2000).
#' @param seed RNG seed for reproducibility (default 1).
#'
#' @return A data frame with one row per value in `effect_grid`, containing:
#' \describe{
#'   \item{effect_size}{The effect size \eqn{\delta}.}
#'   \item{type1_empirical}{Empirical power (or Type I error when \eqn{\delta = 0}).}
#'   \item{type1_exact}{Analytical power from the normal theory formula.}
#'   \item{theta}{True mean of \eqn{Y}.}
#'   \item{theta0}{Null value \eqn{\theta_0}.}
#'   \item{N}{Unlabeled sample size.}
#'   \item{n}{Labeled sample size.}
#'   \item{alpha}{Significance level.}
#'   \item{var_f}{Prediction variance.}
#'   \item{var_res}{Residual variance.}
#'   \item{family}{Distribution family used (`"gaussian"`).}
#' }
#'
#' @examples
#' type1_error_curve_mean_dgp(
#'   effect_grid = seq(-0.4, 0.4, by = 0.05),
#'   N           = 4000,
#'   n           = 200,
#'   var_f       = 0.45,
#'   var_res     = 0.80,
#'   theta       = 0,
#'   R           = 2000
#' )
#'
#' @export
type1_error_curve_mean_dgp <- function(effect_grid,
                                       N,
                                       n,
                                       var_f,
                                       var_res,
                                       theta = 0,    # can just be 0, since only delta matters
                                       alpha = 0.05,
                                       R = 2000,
                                       seed = 1) {
  if (!is.numeric(effect_grid) || length(effect_grid) == 0L) {
    stop("effect_grid must be a numeric vector of effect sizes.", call. = FALSE)
  }
  if (missing(var_f) || missing(var_res)) {
    stop("You must supply var_f and var_res when not using superpopulation.", call. = FALSE)
  }

  effect_grid <- sort(effect_grid)
  seeds <- seed + seq_along(effect_grid)

  sims <- lapply(
    seq_along(effect_grid),
    function(i) {
      delta_i   <- effect_grid[i]
      theta0_i  <- theta - delta_i

      moments_ppi <- list(
        delta   = delta_i,
        var_f   = var_f,
        var_res = var_res
      )

      simulate_ppi_vanilla_mean(
        R        = as.integer(R),
        n        = n,
        N        = N,
        alpha    = alpha,
        family   = stats::gaussian(),
        moments  = moments_ppi,
        theta0   = theta0_i,
        seed     = seeds[i]
      )
    }
  )

  # Extract empirical and analytical power
  type1_empirical <- vapply(sims, function(x) x$empirical_power, numeric(1))
  type1_exact     <- vapply(sims, function(x) x$analytical_power, numeric(1))
  theta0_vals     <- theta - effect_grid

  # Build output data frame
  data.frame(
    effect_size     = effect_grid,
    type1_empirical = type1_empirical,
    type1_exact     = type1_exact,
    theta           = rep(theta, length(effect_grid)),
    theta0          = theta0_vals,
    N               = rep(N, length(effect_grid)),
    n               = rep(n, length(effect_grid)),
    alpha           = rep(alpha, length(effect_grid)),
    var_f           = rep(var_f, length(effect_grid)),
    var_res         = rep(var_res, length(effect_grid)),
    family          = rep("gaussian", length(effect_grid)),
    row.names       = NULL
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
      col = pppower_colors$empirical,
      lwd = 2.5,
      xlab = xlab_default,
      ylab = ylab_default,
      ...
    )
  } else if (exact && !is.null(y_exact)) {
    graphics::plot(
      x,
      y_exact,
      type = "l",
      col = pppower_colors$analytical,
      lwd = 2.5,
      lty = 2,
      xlab = xlab_default,
      ylab = ylab_default,
      ...
    )
  }

  if (empirical && exact && !is.null(y_emp) && !is.null(y_exact)) {
    graphics::lines(x, y_emp, col = pppower_colors$empirical, lwd = 2.5)
    graphics::lines(x, y_exact, col = pppower_colors$analytical, lwd = 2.5, lty = 2)
    graphics::legend(
      legend_pos,
      legend = c("Empirical", "Analytical"),
      col = c(pppower_colors$empirical, pppower_colors$analytical),
      lty = c(1, 2),
      lwd = 2.5,
      bty = "n"
    )
  } else {
    if (empirical && !is.null(y_emp)) {
      graphics::lines(x, y_emp, col = pppower_colors$empirical, lwd = 2.5)
    }
    if (exact && !is.null(y_exact)) {
      graphics::lines(x, y_exact, col = pppower_colors$analytical, lwd = 2.5, lty = 2)
    }
  }

  if (add_reference && "alpha" %in% names(curve_df)) {
    graphics::abline(h = unique(curve_df$alpha), col = pppower_colors$reference, lty = 3)
  }

  invisible(curve_df)
}
