# Compare pppower's design functions with classical double sampling formulas.
#   Two-sample and paired continuous: Cochran's double sampling variance,
#     sigma^2 (1 - rho^2) / n + rho^2 sigma^2 / (n + N), per group or for the
#     paired differences.
#   Two-sample binary and 2x2 (RR, OR): Tenenbein's variance per exposure group,
#     p q {(1 - K) / n + K / (n + N)}, with reliability K, combined by the
#     delta method for log RR and log OR.
# Usage: Rscript simulation_studies/external_checks/check_designs_classical.R
suppressPackageStartupMessages(library(pppower))
z <- qnorm(0.975)
pw <- function(delta, V) pnorm(-z + abs(delta) / sqrt(V)) + pnorm(-z - abs(delta) / sqrt(V))
cochran_var <- function(s2, rho, n, Ntot) s2 * (1 - rho^2) / n + rho^2 * s2 / Ntot
tenenbein_K <- function(p, sens, spec) {
  pf <- sens * p + (1 - spec) * (1 - p)
  p * (sens^2 / pf + (1 - sens)^2 / (1 - pf) - 1) / (1 - p)
}
tenenbein_var <- function(p, K, n, Ntot) p * (1 - p) * ((1 - K) / n + K / Ntot)
solve_n <- function(Vfun, delta, power = 0.8) {
  S2 <- (delta / (z + qnorm(power)))^2
  ceiling(uniroot(function(n) Vfun(n) - S2, c(2, 1e7))$root)
}
out <- list()

# 1. Two-sample continuous (power at fixed n, unequal groups)
for (s in list(c(0.5, 0.7, 40, 60, 300, 500), c(0.9, 0.3, 80, 80, 1000, 200))) {
  rA <- s[1]; rB <- s[2]; nA <- s[3]; nB <- s[4]; NA_ <- s[5]; NB <- s[6]
  p_pp <- power_ppi_ttest(delta = 0.4, n_A = nA, n_B = nB, N_A = NA_, N_B = NB,
                          sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = rA,
                          sigma_y2_B = 1.5, sigma_f2_B = 1, cov_yf_B = rB * sqrt(1.5))
  p_cl <- pw(0.4, cochran_var(1, rA, nA, nA + NA_) + cochran_var(1.5, rB, nB, nB + NB))
  out[[length(out) + 1]] <- data.frame(design = "two-sample mean", setting = sprintf("rho=(%.1f,%.1f) n=(%d,%d) N=(%d,%d)", rA, rB, nA, nB, NA_, NB),
                                       quantity = "power", pppower = p_pp, classical = p_cl)
}

# 2. Paired continuous (required n)
for (s in list(c(0.5, 500), c(0.8, 2000))) {
  n_pp <- ceiling(suppressWarnings(power_ppi_paired(delta = 0.3, N = s[2], n = NULL, power = 0.8,
                                                    sigma_D2 = 1, rho_D = s[1])))
  n_cl <- solve_n(function(n) cochran_var(1, s[1], n, n + s[2]), 0.3)
  out[[length(out) + 1]] <- data.frame(design = "paired mean", setting = sprintf("rho_D=%.1f N=%d", s[1], s[2]),
                                       quantity = "n", pppower = n_pp, classical = n_cl)
}

# 3. Two-sample binary (power at fixed n)
for (s in list(c(0.35, 0.25, 0.85, 100, 1000), c(0.10, 0.05, 0.95, 300, 3000))) {
  pA <- s[1]; pB <- s[2]; a <- s[3]; n <- s[4]; N <- s[5]
  p_pp <- power_ppi_ttest_binary(p_A = pA, p_B = pB, n_A = n, n_B = n, N_A = N, N_B = N,
                                 sens_A = a, spec_A = a, delta = pA - pB)
  V <- tenenbein_var(pA, tenenbein_K(pA, a, a), n, n + N) + tenenbein_var(pB, tenenbein_K(pB, a, a), n, n + N)
  out[[length(out) + 1]] <- data.frame(design = "two-sample binary", setting = sprintf("p=(%.2f,%.2f) acc=%.2f n=%d N=%d", pA, pB, a, n, N),
                                       quantity = "power", pppower = p_pp, classical = pw(pA - pB, V))
}

# 4. 2x2 table, relative risk and odds ratio (required total n)
for (s in list(c(0.40, 0.20, 0.80, 1000, 0.5), c(0.30, 0.20, 0.90, 2000, 0.3))) {
  p1 <- s[1]; p0 <- s[2]; a <- s[3]; N <- s[4]; pi <- s[5]
  K1 <- tenenbein_K(p1, a, a); K0 <- tenenbein_K(p0, a, a)
  Vg <- function(n, p, K, frac) tenenbein_var(p, K, frac * n, frac * (n + N))
  for (m in c("RR", "OR")) {
    n_pp <- suppressWarnings(power_ppi_2x2(p_exp = p1, p_ctrl = p0, N = N, n = NULL, power = 0.8,
                                           prev_exp = pi, sens = a, spec = a, effect_measure = m))
    if (m == "RR") {
      n_cl <- solve_n(function(n) Vg(n, p1, K1, pi) / p1^2 + Vg(n, p0, K0, 1 - pi) / p0^2, log(p1 / p0))
    } else {
      n_cl <- solve_n(function(n) Vg(n, p1, K1, pi) / (p1 * (1 - p1))^2 +
                        Vg(n, p0, K0, 1 - pi) / (p0 * (1 - p0))^2,
                      log(p1 / (1 - p1)) - log(p0 / (1 - p0)))
    }
    out[[length(out) + 1]] <- data.frame(design = paste("2x2", m), setting = sprintf("p=(%.2f,%.2f) acc=%.2f N=%d exposed=%.1f", p1, p0, a, N, pi),
                                         quantity = "n", pppower = n_pp, classical = n_cl)
  }
}

# 5. PROCOVA (Schuler et al. 2022; EMA handbook): randomized trial, 1:1, every
#    participant has the outcome and a prognostic score, nu^2 = (2 sigma)^2 (1 - R^2) / n.
#    This is the two-sample PPI++ design with N -> infinity.
for (s in list(c(0.5, 100), c(0.8, 300))) {
  R <- s[1]; n <- s[2]; nu <- sqrt(4 * (1 - R^2) / n)
  p_cl <- pnorm(qnorm(0.025) + 0.3 / nu) + pnorm(qnorm(0.025) - 0.3 / nu)
  p_pp <- power_ppi_ttest(delta = 0.3, n_A = n / 2, n_B = n / 2, N_A = 1e9, N_B = 1e9,
                          sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = R,
                          sigma_y2_B = 1, sigma_f2_B = 1, cov_yf_B = R)
  out[[length(out) + 1]] <- data.frame(design = "PROCOVA (N -> inf)", setting = sprintf("R=%.1f n=%d", R, n),
                                       quantity = "power", pppower = p_pp, classical = p_cl)
}

res <- do.call(rbind, out)
res$difference <- res$pppower - res$classical
print(res, digits = 4, row.names = FALSE)
