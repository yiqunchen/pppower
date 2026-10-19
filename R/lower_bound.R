#' Lower confidence bound for the prediction-outcome correlation
#'
#' Planning with the point estimate of \eqn{\rho_{Yf}} from a small pilot gives
#' an underpowered study about half of the time, because the estimate is as
#' likely to be too high as too low. This function returns the one-sided Fisher
#' \eqn{z} lower bound
#' \deqn{\tanh\{\operatorname{atanh}(|\hat\rho|) - z_{c}/\sqrt{m - 3}\},}
#' truncated at 0 and carrying the sign of \eqn{\hat\rho}, where \eqn{m} is the
#' number of labeled pairs in the pilot and \eqn{c} is the confidence level.
#' Square it to get a lower bound for \eqn{R^2}.
#'
#' The power functions [power_ppi_mean()], [power_ppi_paired()], and
#' [power_ppi_ttest()] apply this bound when a pilot size is given.
#'
#' @param rho_hat Pilot estimate of the correlation between the outcome and the
#'   predictions.
#' @param pilot_n Number of labeled pairs in the pilot (more than 3).
#' @param confidence One-sided confidence level (default 0.9).
#'
#' @return The lower bound, with the sign of `rho_hat`.
#' @export
#'
#' @examples
#' rho_lower_bound(0.7, pilot_n = 50)
#' rho_lower_bound(0.7, pilot_n = 50)^2   # lower bound for R^2
rho_lower_bound <- function(rho_hat, pilot_n, confidence = 0.9) {
  if (length(pilot_n) != 1 || !is.finite(pilot_n) || pilot_n <= 3) {
    stop("pilot_n must be a single number greater than 3.", call. = FALSE)
  }
  if (!is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must be between 0 and 1.", call. = FALSE)
  }
  if (any(abs(rho_hat) > 1)) stop("rho_hat must lie in [-1, 1].", call. = FALSE)
  r <- pmin(abs(rho_hat), 1 - 1e-12)
  lower <- pmax(tanh(atanh(r) - stats::qnorm(confidence) / sqrt(pilot_n - 3)), 0)
  sign(rho_hat) * lower
}

# Replace Cov(Y, f) by its value at the lower confidence bound for rho.
cov_at_lower_bound <- function(sigma_y2, sigma_f2, cov_y_f, pilot_n, confidence) {
  scale <- sqrt(sigma_y2 * sigma_f2)
  rho_lower_bound(cov_y_f / scale, pilot_n, confidence) * scale
}
