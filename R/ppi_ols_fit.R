#' Fit Prediction-Powered OLS (user-facing)
#'
#' @param y        Name of outcome in labeled data (string or symbol)
#' @param f_col    Name of column(s) containing predictions (in labeled/unlabeled)
#' @param formula  RHS-only formula for covariates, e.g. ~ x1 + x2
#' @param data_l   Labeled data.frame (must contain y, f_col, covariates)
#' @param data_u   Unlabeled data.frame (must contain f_col, covariates)
#' @return An object of class \code{"ppi_ols_fit"} with coef, vcov, se, n, N, call.
#' @export
ppi_ols_fit <- function(y, f_col, formula, data_l, data_u) {
  y <- rlang::as_name(rlang::ensym(y))

  X_l <- stats::model.matrix(formula, data_l)
  X_u <- stats::model.matrix(formula, data_u)

  Y_l <- data_l[[y]]
  f_l <- data_l[[f_col]]
  f_u <- data_u[[f_col]]

  fit <- ppi_ols(X_l, Y_l, f_l, X_u, f_u)

  structure(
    list(
      coef   = fit$theta_hat,
      vcov   = fit$V_theta,
      se     = fit$se,
      pieces = fit$pieces,  # keeps V_u, V_l, n, N for downstream power calcs
      n      = nrow(X_l),
      N      = nrow(X_u),
      call   = match.call(),
      terms  = stats::terms(formula)
    ),
    class = "ppi_ols_fit"
  )
}

#' @export
coef.ppi_ols_fit <- function(object, ...) object$coef

#' @export
vcov.ppi_ols_fit <- function(object, ...) object$vcov

#' @export
summary.ppi_ols_fit <- function(object, ...) {
  z  <- object$coef / object$se
  p  <- 2 * stats::pnorm(-abs(z))
  out <- data.frame(Estimate = object$coef, Std.Error = object$se, z, `Pr(>|z|)` = p)
  attr(out, "n") <- object$n; attr(out, "N") <- object$N
  out
}

#' @export
confint.ppi_ols_fit <- function(object, parm = seq_along(object$coef), level = 0.95, ...) {
  a <- (1 - level)/2; z <- stats::qnorm(1 - a)
  b <- cbind(object$coef[parm] - z*object$se[parm],
             object$coef[parm] + z*object$se[parm])
  colnames(b) <- c("lower","upper"); b
}

#' Wald test for a linear contrast of OLS
#'
#' @param fit A `ppi_ols_fit` object returned by [ppi_ols_fit()].
#' @param A Numeric matrix whose rows encode the linear contrast(s).
#' @param b Numeric value or vector specifying the null hypothesis (defaults to 0).
#'
#' @return A list with elements `estimate`, `se`, `z`, and `p.value`.
#' @export
ppi_ols_wald <- function(fit, A, b = 0) {
  th <- drop(A %*% fit$coef)
  V  <- drop(A %*% fit$vcov %*% t(A))
  se <- sqrt(V); z <- (th - b)/se; p <- 2*stats::pnorm(-abs(z))
  list(estimate = th, se = se, z = z, p.value = p)
}