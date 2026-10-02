# Check the browser calculator (vignettes/sample-size-calculator.Rmd) against
# pppower for the difference-in-means, paired, and paired-proportion panels.
# Needs node. Usage: Rscript simulation_studies/external_checks/check_calculator.R
suppressPackageStartupMessages(library(pppower))
harness <- file.path("simulation_studies", "external_checks", "calculator_harness.js")
run_js <- function(over) {
  out <- system2("node", c(harness, shQuote(jsonlite::toJSON(over, auto_unbox = TRUE))), stdout = TRUE)
  jsonlite::fromJSON(out)
}
num <- function(x) as.numeric(gsub(",", "", x))
solve_n <- function(pf, target = 0.8) pppower:::.mcp_solve_n(pf, target)
z <- function(a, pw) (qnorm(1 - a / 2) + qnorm(pw))
rows <- list()

# difference in means
for (s in list(list(delta = 0.3, sd_A = 1, sd_B = 1, r2_A = 0.5, r2_B = 0.5, N_A = 1000, N_B = 1000),
               list(delta = 2, sd_A = 5, sd_B = 7, r2_A = 0.2, r2_B = 0.4, N_A = 300, N_B = 2000),
               list(delta = 0.5, sd_A = 1, sd_B = 1, r2_A = 0.9, r2_B = 0.1, N_A = 5000, N_B = 500))) {
  js <- run_js(setNames(s, paste0("dm_", names(s))))$dm
  pw <- function(n, rA, rB) power_ppi_ttest(delta = s$delta, n_A = n, n_B = n, N_A = s$N_A, N_B = s$N_B,
    sigma_y2_A = s$sd_A^2, sigma_f2_A = s$sd_A^2, cov_yf_A = sqrt(rA) * s$sd_A^2,
    sigma_y2_B = s$sd_B^2, sigma_f2_B = s$sd_B^2, cov_yf_B = sqrt(rB) * s$sd_B^2)
  rows[[length(rows) + 1]] <- data.frame(panel = "difference in means",
    setting = paste(names(s), unlist(s), sep = "=", collapse = " "),
    js_ppi = num(js$n_ppi), r_ppi = solve_n(function(n) pw(n, s$r2_A, s$r2_B)),
    js_classical = num(js$n_classic), r_classical = solve_n(function(n) pw(n, 0, 0)))
}

# paired mean difference
for (s in list(list(delta = 0.3, sd_d = 1, r2_d = 0.25, N = 500),
               list(delta = 1, sd_d = 4, r2_d = 0.6, N = 3000),
               list(delta = 0.2, sd_d = 1, r2_d = 0.05, N = 200))) {
  js <- run_js(setNames(s, paste0("pm_", names(s))))$pm
  r_ppi <- ceiling(suppressWarnings(power_ppi_paired(delta = s$delta, N = s$N, n = NULL, power = 0.8,
                                                     sigma_D2 = s$sd_d^2, rho_D = sqrt(s$r2_d))))
  rows[[length(rows) + 1]] <- data.frame(panel = "paired mean difference",
    setting = paste(names(s), unlist(s), sep = "=", collapse = " "),
    js_ppi = num(js$n_ppi), r_ppi = r_ppi,
    js_classical = num(js$n_classic), r_classical = ceiling(s$sd_d^2 * z(0.05, 0.8)^2 / s$delta^2))
}

# paired proportions
for (s in list(list(pA = 40, pB = 30, rho = 0.3, sens = 0.85, spec = 0.85, N = 1000),
               list(pA = 20, pB = 12, rho = 0.5, sens = 0.95, spec = 0.9, N = 5000),
               list(pA = 55, pB = 50, rho = 0.1, sens = 0.7, spec = 0.7, N = 800))) {
  js <- run_js(setNames(s, paste0("pb_", names(s))))$pb
  pA <- s$pA / 100; pB <- s$pB / 100
  pf <- function(n) suppressWarnings(power_ppi_paired_binary(delta = pA - pB, N = s$N, n = n, p_A = pA, p_B = pB,
                                                             rho_within = s$rho, sens = s$sens, spec = s$spec))
  p11 <- pA * pB + s$rho * sqrt(pA * (1 - pA) * pB * (1 - pB))
  sD2 <- (pA - p11) + (pB - p11) - (pA - pB)^2
  rows[[length(rows) + 1]] <- data.frame(panel = "paired proportions",
    setting = paste(names(s), unlist(s), sep = "=", collapse = " "),
    js_ppi = num(js$n_ppi), r_ppi = solve_n(pf),
    js_classical = num(js$n_classic), r_classical = ceiling(sD2 * z(0.05, 0.8)^2 / (pA - pB)^2))
}

res <- do.call(rbind, rows)
res$diff_ppi <- res$js_ppi - res$r_ppi
res$diff_classical <- res$js_classical - res$r_classical
print(res, row.names = FALSE)
cat("max |difference| =", max(abs(c(res$diff_ppi, res$diff_classical))), "\n")
