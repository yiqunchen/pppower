#!/usr/bin/env Rscript
# =============================================================================
# Planning with a lower confidence bound for the prediction-outcome correlation
# =============================================================================
# Usage: Rscript simulation_studies/lower_bound_r2.R
# For prototypical one-sample designs, draw a pilot of m labeled pairs, plan n*
# with the pilot estimate rho_hat and with the one-sided 90% Fisher z lower
# bound tanh(atanh(rho_hat) - z_0.9 / sqrt(m - 3)), and compute the power that
# each plan achieves at the true rho (with the closed-form power of Section 3.1,
# which Section 5 validates). Writes
#   lower_bound_r2/lower_bound_r2.csv
#   table_lower_bound_r2.tex

script_dir <- if (dir.exists("simulation_studies")) "simulation_studies" else "."
out_dir <- file.path(script_dir, "lower_bound_r2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

R_reps <- 10000
alpha <- 0.05; target <- 0.80; delta <- 0.2; N <- 1000; sigma_y2 <- 1
confidence <- 0.90
z <- qnorm(1 - alpha / 2) + qnorm(target)
S2 <- (delta / z)^2
set.seed(20261019)

# closed-form n* (Section 3.1) and the power at a given n under the true rho
n_star <- function(rho) {
  b <- sigma_y2 - S2 * N
  ceiling((b + sqrt(b^2 + 4 * S2 * N * sigma_y2 * (1 - rho^2))) / (2 * S2))
}
power_at <- function(n, rho) {
  v <- sigma_y2 / n - rho^2 * sigma_y2 * N / (n * (n + N))
  pnorm(-qnorm(1 - alpha / 2) + delta / sqrt(v)) + pnorm(-qnorm(1 - alpha / 2) - delta / sqrt(v))
}
# sample correlation of m bivariate normal pairs, R_reps times
pilot_rho <- function(rho, m) {
  y <- matrix(rnorm(m * R_reps), m)
  f <- rho * y + sqrt(1 - rho^2) * matrix(rnorm(m * R_reps), m)
  yc <- sweep(y, 2, colMeans(y)); fc <- sweep(f, 2, colMeans(f))
  colSums(yc * fc) / sqrt(colSums(yc^2) * colSums(fc^2))
}

n_cl <- ceiling(sigma_y2 / S2)
rows <- list()
for (rho in c(0.5, 0.7, 0.9)) {
  for (m in c(30, 50, 100, 200)) {
    r_hat <- pmin(pmax(pilot_rho(rho, m), 0), 0.999)
    r_low <- pmax(tanh(atanh(r_hat) - qnorm(confidence) / sqrt(m - 3)), 0)
    for (plan in c("point", "lower")) {
      r_plan <- if (plan == "point") r_hat else r_low
      n <- n_star(r_plan)
      pw <- power_at(n, rho)
      rows[[length(rows) + 1]] <- data.frame(
        rho = rho, r2 = rho^2, m = m, plan = plan, n_true = n_star(rho), n_cl = n_cl,
        mean_n = mean(n), mean_power = mean(pw), p10_power = unname(quantile(pw, 0.1)),
        prob_target = mean(pw >= target - 0.005), R = R_reps)
    }
  }
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out_dir, "lower_bound_r2.csv"), row.names = FALSE)

# table: one row per (R^2, m), point estimate versus lower bound
w <- reshape(res[, c("r2", "m", "plan", "n_true", "n_cl", "mean_n", "mean_power", "p10_power", "prob_target")],
             idvar = c("r2", "m", "n_true", "n_cl"), timevar = "plan", direction = "wide")
w <- w[order(w$r2, w$m), ]
tex <- c("\\begin{tabular}{rrrrrrrrr}", "\\toprule",
         " & & & \\multicolumn{3}{c}{Plan with $\\hat\\rho_{Yf}$} & \\multicolumn{3}{c}{Plan with the lower bound} \\\\",
         "\\cmidrule(lr){4-6}\\cmidrule(lr){7-9}",
         "$R^2$ & $m$ & $n^\\star$ & mean $n$ & power & $P(\\geq 0.80)$ & mean $n$ & power & $P(\\geq 0.80)$ \\\\",
         "\\midrule")
for (i in seq_len(nrow(w))) {
  x <- w[i, ]
  tex <- c(tex, sprintf("%.2f & %d & %d & %.0f & %.2f & %.2f & %.0f & %.2f & %.2f \\\\",
                        x$r2, x$m, x$n_true, x$mean_n.point, x$mean_power.point, x$prob_target.point,
                        x$mean_n.lower, x$mean_power.lower, x$prob_target.lower))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_lower_bound_r2.tex"))
print(res)
cat(sprintf("n_cl = %d\n", n_cl))
