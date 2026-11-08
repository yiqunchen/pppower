# Cross-check plugin PPI++ against the two-phase (double sampling) regression
# estimator in the survey package (Lumley). The GREG/calibration estimator of
# a two-phase design with f as the phase-1 auxiliary variable and simple random
# phase-2 sampling is algebraically identical to PPI++ with the plugin lambda.
# Usage: Rscript simulation_studies/external_checks/check_survey_twophase.R
suppressMessages(library(survey))
set.seed(11)

one <- function(n, N, rho, binary = FALSE) {
  M <- n + N
  if (binary) {
    Y <- rbinom(M, 1, 0.3); f <- ifelse(Y == 1, rbinom(M, 1, 0.9), rbinom(M, 1, 0.1))
  } else {
    z1 <- rnorm(M); Y <- z1; f <- rho * z1 + sqrt(1 - rho^2) * rnorm(M)
  }
  lab <- seq_len(M) <= n
  d <- data.frame(id = 1:M, Y = ifelse(lab, Y, NA), f = f, lab = lab)
  des <- suppressWarnings(twophase(id = list(~id, ~id), subset = ~lab, data = d, method = "simple"))
  est <- suppressWarnings(svymean(~Y, calibrate(des, phase = 2, formula = ~f, calfun = "linear")))
  b <- cov(Y[lab], f[lab]) / var(f[lab])
  ppi <- mean(Y[lab]) + b * N / (n + N) * (mean(f[!lab]) - mean(f[lab]))
  c(survey = coef(est)[[1]], ppi = ppi, se = SE(est)[[1]])
}

cat("Point estimates (survey GREG vs PPI++ plugin):\n")
for (s in list(c(30, 500, 0.7, 0), c(80, 2000, 0.5, 0), c(60, 1000, 0, 1))) {
  r <- one(s[1], s[2], s[3], s[4] == 1)
  cat(sprintf("  n=%d N=%d binary=%d: diff = %.1e\n", s[1], s[2], s[4], abs(r[1] - r[2])))
}

cat("Variance (Gaussian, rho = 0.7, N = 1000):\n")
rho <- 0.7; N <- 1000
for (n in c(20, 50, 100)) {
  R <- 400
  rr <- t(replicate(R, one(n, N, rho)))
  ok <- is.finite(rr[, 3])
  Vopt <- (1 - rho^2) / n + rho^2 / (n + N)
  Vpl <- Vopt + (N / (n + N)) * (1 - rho^2) / (n * (n - 3))
  cat(sprintf("  n=%3d: mean survey SE^2 = %.5f (%d/%d finite), V_opt = %.5f, V_plugin = %.5f\n",
              n, mean(rr[ok, 3]^2), sum(ok), R, Vopt, Vpl))
}
