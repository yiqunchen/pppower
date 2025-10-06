#' Solve for required n (labeled sample size) for a desired power
#'
#' @param delta Effect size (θ - θ_0)
#' @param var_f Variance of f(X), the predictor function
#' @param var_res Variance of (Y - f(X)), the residual
#' @param N Unlabeled sample size
#' @param alpha Significance level
#' @param power Target power (1 - β)
#' @return Required labeled sample size n
#' @export
#' @examples
#' n_required_PP(delta = 0.1, var_f = 0.2, var_res = 0.5, N = 5000, alpha = 0.05, power = 0.8)
n_required_PP <- function(delta,
                          var_f,
                          var_res,
                          N,
                          alpha = 0.05,
                          power = 0.80,
                          warn_smallN = TRUE,
                          smallN_threshold = 500,
                          mode = c("error", "cap")) {
  mode <- match.arg(mode)

  # Optional advisory for small unlabeled pool
  if (warn_smallN && is.finite(N) && N < smallN_threshold) {
    warning("Unlabeled N is too small. The finite-N term Var(f)/N may be non-negligible.",
            call. = FALSE)
  }

  c0 <- qnorm(1 - alpha/2) + qnorm(power)
  denom <- (delta^2 / c0^2) - (var_f / N)

  if (denom <= 0) {
    stop("Infeasible: unlabeled term dominates; increase N or intended effect size (delta), ",
         "or lower target power / increase alpha.", call. = FALSE)
  }

  n_star <- ceiling(var_res / denom)

  if (n_star <= N) {
    return(n_star)
  }

  # n_star > N: handle per mode
  if (mode == "error") {
    stop(sprintf("Infeasible: required n = %d exceeds N = %d. ",
                 n_star, N),
         "Increase N or delta, or relax power/alpha. ",
         "Alternatively, set mode='cap' to return n = N with achieved power.",
         call. = FALSE)
  }

  # mode == "cap": cap at N and report achieved power as an attribute
  n_capped <- as.integer(N)
  ach_power <- power_ppi_mean(delta, var_f, var_res, N, n_capped, alpha)
  warning(sprintf("Required n = %d exceeds N = %d. Capping to n = N.\nAchieved power = %.4f (target: %.4f).",
                  n_star, N, ach_power, power),
          call. = FALSE)
  # return n with an attribute for user to inspect
  attr(n_capped, "achieved_power") <- ach_power
  n_capped
}
