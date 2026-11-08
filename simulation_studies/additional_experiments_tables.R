#!/usr/bin/env Rscript
# =============================================================================
# Tables for the three additional experiments (revision R1)
# =============================================================================
# Usage: Rscript simulation_studies/additional_experiments_tables.R
# Reads additional_experiments/*.csv (from run_additional_experiments.R) and writes
#   table_add_low_r2.tex       low prediction quality: variance ratio and Type I error
#   table_add_pilot_shift.tex  pilot-to-study shift: power and Type I error at n*
#   table_add_skew.tex         Type I error versus n and skewness

script_dir <- if (dir.exists("simulation_studies")) "simulation_studies" else "."
in_dir <- file.path(script_dir, "additional_experiments")
f2 <- function(z) sprintf("%.2f", z)
f3 <- function(z) sprintf("%.3f", z)
row_tex <- function(...) paste0(paste(..., sep = " & "), " \\\\")

# ---------------------------------------------------------------------------
# Low R^2: Var(plugin) / Var(classical) under H0, and the largest Type I error
# over the R^2 grid
# ---------------------------------------------------------------------------
x <- read.csv(file.path(in_dir, "low_r2_crossover.csv"))
x0 <- subset(x, hyp == "H0")
r2_grid <- sort(unique(x0$r2))
cells <- unique(x0[, c("design", "n")])
cells <- cells[order(cells$design, cells$n), ]
tex <- c("\\begin{tabular}{llrrrrrr}", "\\toprule",
         " & & \\multicolumn{4}{c}{Variance of plug-in / classical estimator} & \\multicolumn{2}{c}{Max.\\ Type I error} \\\\",
         "\\cmidrule(lr){3-6}\\cmidrule(lr){7-8}",
         paste0("Design & $n$ & $R^2=", r2_grid[1], "$ & ",
                paste0("$", r2_grid[-1], "$", collapse = " & "), " & plug-in & corrected \\\\"),
         "\\midrule")
for (i in seq_len(nrow(cells))) {
  s <- subset(x0, design == cells$design[i] & n == cells$n[i])
  ratio <- sapply(r2_grid, function(r) {
    s$mc_var[s$r2 == r & s$method == "plugin"] / s$mc_var[s$r2 == r & s$method == "classical"]
  })
  tex <- c(tex, row_tex(cells$design[i], cells$n[i], paste(f2(ratio), collapse = " & "),
                        f3(max(s$reject[s$method == "plugin"])),
                        f3(max(s$reject[s$method == "plugin_adj"]))))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_add_low_r2.tex"))

# ---------------------------------------------------------------------------
# Pilot-to-study shift: n* planned with rho_pilot, data from rho_study
# ---------------------------------------------------------------------------
x <- read.csv(file.path(in_dir, "pilot_study_shift.csv"))
cells <- unique(x[, c("design", "rho_pilot", "rho_study", "n")])
cells <- cells[order(cells$design, cells$rho_pilot, cells$rho_study), ]
pick <- function(s, h, m) s$reject[s$hyp == h & s$method == m]
tex <- c("\\begin{tabular}{llrrrrrr}", "\\toprule",
         " & & & & \\multicolumn{2}{c}{Power} & \\multicolumn{2}{c}{Type I error} \\\\",
         "\\cmidrule(lr){5-6}\\cmidrule(lr){7-8}",
         "Design & $\\rho_{\\mathrm{pilot}}$ & $\\rho_{\\mathrm{study}}$ & $n^\\star$ & corrected & pilot $\\lambda$ & corrected & pilot $\\lambda$ \\\\",
         "\\midrule")
for (i in seq_len(nrow(cells))) {
  s <- merge(x, cells[i, ])
  tex <- c(tex, row_tex(cells$design[i], cells$rho_pilot[i], cells$rho_study[i], cells$n[i],
                        f3(pick(s, "H1", "plugin_adj")), f3(pick(s, "H1", "pilot")),
                        f3(pick(s, "H0", "plugin_adj")), f3(pick(s, "H0", "pilot"))))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_add_pilot_shift.tex"))

# ---------------------------------------------------------------------------
# Skewness: Type I error at four of the six sample sizes
# ---------------------------------------------------------------------------
x <- read.csv(file.path(in_dir, "skew_calibration.csv"))
x <- subset(x, n %in% c(15, 30, 100, 200))
outcome <- function(s) ifelse(s == 0, "Gaussian", sprintf("Log-normal (sd %g)", s))
cells <- unique(x[, c("sigma_log", "design", "n")])
cells <- cells[order(cells$sigma_log, cells$design, cells$n), ]
methods <- c("classical", "oracle", "plugin", "plugin_adj", "pilot")
tex <- c("\\begin{tabular}{llrrrrrr}", "\\toprule",
         "Outcome & Design & $n$ & classical & oracle $\\lambda^\\star$ & plug-in & corrected & pilot $\\lambda$ \\\\",
         "\\midrule")
for (i in seq_len(nrow(cells))) {
  s <- merge(x, cells[i, ])
  tex <- c(tex, row_tex(outcome(cells$sigma_log[i]), cells$design[i], cells$n[i],
                        paste(f3(s$reject[match(methods, s$method)]), collapse = " & ")))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_add_skew.tex"))
cat("Additional experiment tables written\n")
