
#' Power for prediction-powered linear regression (OLS) estimator
#' @description
#' Computes the analytical (normal-theory) power for a two-sided Wald test on a
#' linear contrast \eqn{c^\top \beta} in the \emph{prediction-powered linear regression (PPI-OLS)} framework.
#' @param delta Numeric scalar. True effect size for the linear contrast
#' \eqn{c^\top (\beta - \beta_0)}.
#' @param V_u Numeric \eqn{p \times p} matrix.
#' Sandwich variance for the unlabeled (imputed) OLS term.
#' @param V_l Numeric \eqn{p \times p} matrix.
#' Sandwich variance for the labeled (rectifier) OLS term.
#' @param N Integer. Unlabeled sample size.
#' @param n Integer. Labeled sample size.
#' @param c Numeric vector of length \eqn{p}.
#' Linear contrast to test; e.g. to test coefficient \eqn{j}, set \code{c = e_j}.
#' @param alpha Numeric scalar in (0, 1). Two-sided significance level.
#'
#' @return
#' Numeric scalar: the analytical (normal-theory) power \eqn{1-\beta}
#' of the two-sided test.
#' 
#' @examples
#' # Example: 2D regression with modest effect
#' set.seed(1)
#' p <- 2
#' V_u <- matrix(c(1.0, 0.2, 0.2, 0.8), ncol = p)
#' V_l <- matrix(c(0.9, 0.1, 0.1, 0.7), ncol = p)
#' N <- 2000
#' n <- 200
#' c <- c(1, 0)     # test first coefficient
#' delta <- 0.15
#'
#' power_ppi_ols(delta, V_u, V_l, N, n, c, alpha = 0.05)
#'
#' # Test second coefficient (different direction)
#' c2 <- c(0, 1)
#' power_ppi_ols(delta = 0.10, V_u, V_l, N, n, c2)
#' 
#' @export

power_ppi_ols <- function(delta, V_u, V_l, N, n, c, alpha = 0.05) {
  # variance of c' theta_hat: v = c'(V_u/N + V_l/n)c
  v <- as.numeric(drop(c %*% (V_u/N + V_l/n) %*% c))
  se <- sqrt(v)
  z_alpha <- stats::qnorm(1 - alpha/2)
  mu <- abs(delta) / se
  # exact two-sided normal power
  1 - stats::pnorm(z_alpha - mu) + stats::pnorm(-z_alpha - mu)
}


#' @keywords internal
#' Monte Carlo power for PPI-OLS (iid sample from the DGP)
simulate_power_ppi_ols <- function(R,
                                   n,
                                   N,
                                   contrast,
                                   alpha = 0.05,
                                   family = stats::gaussian(),
                                   moments = NULL,
                                   seed = 1) {
  if (!is.null(seed)) set.seed(seed)

  # CASE 1: population moments provided 
  if (!is.null(moments)) {
    required <- c("delta", "H", "Sigma_YY")
    missing <- setdiff(required, names(moments))
    if (length(missing) > 0) {
      stop("moments list must include: ", paste(required, collapse = ", "), call. = FALSE)
    }

    delta     <- moments$delta
    H         <- moments$H
    Sigma_YY  <- moments$Sigma_YY

    # Variance of contrast estimator under vanilla PPI
    V_contrast <- drop(t(contrast) %*% solve(H) %*% (Sigma_YY / n) %*% solve(H) %*% contrast)
    se_contrast <- sqrt(V_contrast)

    zcrit <- stats::qnorm(1 - alpha / 2)
    mu <- abs(delta) / se_contrast

    # Analytical power BEFORE MC
    analytical_power <- 1 - stats::pnorm(zcrit - mu) + stats::pnorm(-zcrit - mu)

    # Monte Carlo simulation directly from the implied Wald statistic distribution
    u <- (seq_len(R) - 0.5) / R  # deterministic inverse CDF draws
    theta_diff <- delta + se_contrast * stats::qnorm(u)
    z_stat <- theta_diff / se_contrast
    reject <- abs(z_stat) > zcrit

    details <- lapply(
      seq_len(R),
      function(i) list(
        theta_diff = theta_diff[i],
        z_stat = z_stat[i],
        se = se_contrast,
        reject = reject[i]
      )
    )

    return(list(
      empirical_power   = mean(reject),
      analytical_power  = analytical_power,
      se_contrast       = se_contrast,
      details           = details
    ))
  }

  # population moments not provided 
  results <- replicate(
    R,
    simulate_ppi_ols_rep(
      n_l = n,
      n_u = N,
      contrast = contrast,
      family = family,
      method = "ppi",
      alpha = alpha,
      seed = sample.int(.Machine$integer.max, 1)
    ),
    simplify = FALSE
  )

  list(
    empirical_power   = mean(vapply(results, `[[`, logical(1), "reject")),
    analytical_power  = NULL,
    details           = results
  )
}