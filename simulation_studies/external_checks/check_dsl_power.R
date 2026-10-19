# Compare dsl::power_dsl (Egami et al.) predicted SEs with pppower. Needs dsl:
# remotes::install_github("naoki-egami/dsl"). Usage: Rscript simulation_studies/external_checks/check_dsl_power.R
suppressMessages({library(dsl); library(pppower)})
set.seed(3)
n_tot <- 20000; n_lab <- 2000; rho <- 0.7
z <- rnorm(n_tot); Y <- z; f <- rho * z + sqrt(1 - rho^2) * rnorm(n_tot)
lab <- sample(n_tot, n_lab)
d <- data.frame(Y = NA_real_, f = f, x = rnorm(n_tot))
d$Y[lab] <- Y[lab]
out <- suppressMessages(dsl(model = "lm", formula = Y ~ 1, predicted_var = "Y", prediction = "f",
                            data = d, sl_method = "lm", feature = "f", sample_split = 5, seed = 1))
pw <- suppressMessages(power_dsl(labeled_size = c(3000, 5000, 8000), dsl_out = out))
se_dsl <- pw$predicted_se[, 1]
sizes <- pw$labeled_size
V <- function(n) var_n <- (1 - rho^2) / n + rho^2 / n_tot      # true moments, N = n_tot - n
se_ppp <- sapply(sizes, function(n) sqrt(1 / n - rho^2 * (n_tot - n) / (n * n_tot)))
# pppower power at each size for delta = 0.03, versus DSL-implied power
pw_ppp <- sapply(sizes, function(n) power_ppi_mean(delta = 0.03, N = n_tot - n, n = n, sigma_y2 = 1, sigma_f2 = 1, cov_y_f = rho))
pw_dsl <- pnorm(-qnorm(0.975) + 0.03 / se_dsl) + pnorm(-qnorm(0.975) - 0.03 / se_dsl)
print(data.frame(labeled = sizes, se_dsl = se_dsl, se_pppower = se_ppp, ratio = se_dsl / se_ppp,
                 power_dsl = pw_dsl, power_pppower = pw_ppp), digits = 4)
