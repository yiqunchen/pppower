#!/usr/bin/env Rscript
# =============================================================================
# SETTING 23: ESTIMATED TUNING PARAMETER, SKEWED OUTCOMES, AND SMALL n
# =============================================================================
# Revision analyses for Reviewer 1, comments 1(a) and 3.
#   (a) Power and Type I error when lambda is estimated, for Gaussian,
#       log-normal (two skew levels), t5, binary, and rare binary outcomes,
#       in the one-sample and two-sample designs.
#   (b) Finite-sample variance: Monte Carlo variance of each estimator versus
#       the asymptotic variance and the plug-in corrected variance, the mean
#       estimated variance, and interval coverage.
#   (c) Achieved power at the planned n, planning with the oracle formula and
#       with the plug-in corrected formula.
# Methods: oracle, plugin, plugin_adj, crossfit, pilot, classical
# (see revision_utils.R).
#
# Usage: Rscript simulation_studies/power_validation_setting23_plugin_lambda.R
# =============================================================================

suppressPackageStartupMessages({ library(parallel) })
script_dir <- if (dir.exists("simulation_studies")) "simulation_studies" else "."
source(file.path(script_dir, "revision_utils.R"))

R_reps  <- as.integer(Sys.getenv("SETTING23_R", "2000"))
n_cores <- as.integer(Sys.getenv("SETTING23_CORES", max(1, detectCores() - 2)))
alpha   <- 0.05
methods <- c("oracle", "plugin", "plugin_adj", "crossfit", "pilot", "classical")

cat("Setting 23: Estimated tuning parameter\n")
cat(sprintf("R = %d replicates, %d cores\n", R_reps, n_cores))
t_start <- Sys.time()

# ---------------------------------------------------------------------------
# Distribution families. Each family has a function that returns the spec
# under the null and under the alternative.
# ---------------------------------------------------------------------------
rho_grid <- c(0.3, 0.5, 0.7, 0.9)
acc_grid <- c(0.70, 0.85, 0.95)

families <- list()
for (rho in rho_grid) {
  families[[length(families) + 1]] <- list(family = "Gaussian", type = "gaussian", rho = rho)
  families[[length(families) + 1]] <- list(family = "Log-normal (sd 0.5)", type = "lognormal",
                                           rho = rho, sigma_log = 0.5)
  families[[length(families) + 1]] <- list(family = "Log-normal (sd 1)", type = "lognormal",
                                           rho = rho, sigma_log = 1.0)
  families[[length(families) + 1]] <- list(family = "t5", type = "t5", rho = rho,
                                           rho_gauss = t5_rho_gauss(rho))
}
for (a in acc_grid) {
  families[[length(families) + 1]] <- list(family = "Binary (p = 0.35)", type = "binary",
                                           p0 = 0.35, p1 = 0.50, sens = a, spec = a)
  families[[length(families) + 1]] <- list(family = "Binary (p = 0.05)", type = "binary",
                                           p0 = 0.05, p1 = 0.10, sens = a, spec = a)
}

is_rare <- function(fam) fam$type == "binary" && fam$p0 < 0.1
n_grid_for <- function(fam) if (is_rare(fam)) c(50, 100, 200, 400) else c(15, 20, 30, 50, 100)
N_grid_for <- function(fam) if (is_rare(fam)) c(1000, 5000) else c(200, 1000)
m_pilot_for <- function(fam) if (is_rare(fam)) 1000 else 200

# spec and shift for a hypothesis ("H0" or "H1")
spec_for <- function(fam, hyp) {
  if (fam$type == "binary") {
    d <- fam; d$p <- if (hyp == "H1") fam$p1 else fam$p0
    return(list(dist = d, shift = 0))
  }
  list(dist = fam, shift = 0)
}

# Population moments, computed once per family and hypothesis
cat("Computing population moments...\n")
moments <- lapply(families, function(fam) {
  m0 <- pop_moments(spec_for(fam, "H0")$dist)
  m1 <- pop_moments(spec_for(fam, "H1")$dist)
  delta <- if (fam$type == "binary") fam$p1 - fam$p0 else 0.3 * sqrt(m0$sY2)
  list(H0 = m0, H1 = m1, delta = delta, theta0 = m0$theta)
})

# data for one group: returns per-method estimates
one_group <- function(fam, hyp, n, N, R, mom, extra_shift = 0) {
  sp <- spec_for(fam, hyp)
  shift <- sp$shift + extra_shift
  d <- gen_pairs(sp$dist, n, R, shift)
  fU <- gen_fU_bar(sp$dist, N, R, shift)
  lam_o <- mom$cov / ((1 + n / N) * mom$sf2)
  lam_pi <- pilot_lambda(sp$dist, m_pilot_for(fam), R, n, N, shift)
  ppi_methods(d$Y, d$f, fU, N, lam_o, lam_pi)
}

# ---------------------------------------------------------------------------
# Part (a) and (b): power, Type I error, and variance on the (n, N) grid
# ---------------------------------------------------------------------------
grid <- do.call(rbind, lapply(seq_along(families), function(i) {
  fam <- families[[i]]
  expand.grid(fam_id = i, n = n_grid_for(fam), N = N_grid_for(fam),
              hyp = c("H0", "H1"), design = c("one-sample", "two-sample"),
              stringsAsFactors = FALSE)
}))
cat(sprintf("Grid: %d configurations\n", nrow(grid)))

run_grid_row <- function(k) {
  g <- grid[k, ]; fam <- families[[g$fam_id]]; mo <- moments[[g$fam_id]]
  n <- g$n; N <- g$N; delta <- mo$delta
  if (g$design == "one-sample") {
    mom <- mo[[g$hyp]]
    shift <- if (fam$type != "binary" && g$hyp == "H1") delta else 0
    est <- one_group(fam, g$hyp, n, N, R_reps, mom, extra_shift = shift)
    theta_true <- mom$theta + shift
    theta0 <- mo$theta0
    V_or <- var_opt(mo$H1, n, N); V_pl <- var_plugin(mo$H1, n, N)
    V_or_true <- var_opt(mom, n, N); V_pl_true <- var_plugin(mom, n, N)
  } else {
    # group A follows the alternative under H1 and the null under H0;
    # group B always follows the null
    hypA <- g$hyp
    momA <- mo[[hypA]]; momB <- mo$H0
    shiftA <- if (fam$type != "binary" && hypA == "H1") delta else 0
    eA <- one_group(fam, hypA, n, N, R_reps, momA, extra_shift = shiftA)
    eB <- one_group(fam, "H0", n, N, R_reps, momB)
    est <- lapply(methods, function(mt) {
      a <- eA[[mt]]; b <- eB[[mt]]
      df <- if (mt == "classical") (a$v + b$v)^2 / (a$v^2 / (n - 1) + b$v^2 / (n - 1))
            else if (mt == "plugin_adj") 2 * (n - 2) else Inf
      list(th = a$th - b$th, v = a$v + b$v, df = df)
    })
    names(est) <- methods
    theta_true <- (momA$theta + shiftA) - momB$theta
    theta0 <- 0
    V_or <- var_opt(mo$H1, n, N) + var_opt(mo$H0, n, N)
    V_pl <- var_plugin(mo$H1, n, N) + var_plugin(mo$H0, n, N)
    V_or_true <- var_opt(momA, n, N) + var_opt(momB, n, N)
    V_pl_true <- var_plugin(momA, n, N) + var_plugin(momB, n, N)
  }
  res <- do.call(rbind, lapply(methods, function(mt) {
    s <- summarize_method(est[[mt]], theta0, theta_true, alpha)
    cbind(method = mt, s)
  }))
  cl_V <- if (g$design == "one-sample") mo$H1$sY2 / n else (mo$H1$sY2 + mo$H0$sY2) / n
  cbind(family = fam$family, type = fam$type, rho_target = if (is.null(fam$rho)) NA else fam$rho,
        sens = if (is.null(fam$sens)) NA else fam$sens,
        rho_H1 = mo$H1$rho, kappa = mo$H1$kappa, design = g$design, hyp = g$hyp,
        n = n, N = N, delta = delta, R = R_reps,
        power_theo_oracle = wald_power(delta, V_or, alpha),
        power_theo_plugin = wald_power(delta, V_pl, alpha),
        power_theo_classical = wald_power(delta, cl_V, alpha),
        V_opt = V_or_true, V_plugin = V_pl_true, res,
        stringsAsFactors = FALSE)
}

RNGkind("L'Ecuyer-CMRG"); set.seed(20260930)
res_list <- mclapply(seq_len(nrow(grid)), run_grid_row, mc.cores = n_cores, mc.set.seed = TRUE)
bad <- vapply(res_list, inherits, logical(1), "try-error")
if (any(bad)) { print(res_list[bad][1]); stop("Some configurations failed") }
results_setting23 <- do.call(rbind, res_list)
results_setting23$mc_se <- sqrt(results_setting23$reject * (1 - results_setting23$reject) / R_reps)
cat(sprintf("Grid done (%.1f min)\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))

# ---------------------------------------------------------------------------
# Part (c): achieved power at the planned sample size
# ---------------------------------------------------------------------------
plan_grid <- do.call(rbind, lapply(seq_along(families), function(i) {
  expand.grid(fam_id = i, N = N_grid_for(families[[i]]),
              design = c("one-sample", "two-sample"), stringsAsFactors = FALSE)
}))

run_plan_row <- function(k) {
  g <- plan_grid[k, ]; fam <- families[[g$fam_id]]; mo <- moments[[g$fam_id]]
  N <- g$N; delta <- mo$delta
  if (g$design == "one-sample") {
    Vo <- function(n) var_opt(mo$H1, n, N); Vp <- function(n) var_plugin(mo$H1, n, N)
  } else {
    Vo <- function(n) var_opt(mo$H1, n, N) + var_opt(mo$H0, n, N)
    Vp <- function(n) var_plugin(mo$H1, n, N) + var_plugin(mo$H0, n, N)
  }
  n_or <- max(plan_n(Vo, delta, 0.8, alpha, n_min = 6), 6)
  n_pl <- max(plan_n(Vp, delta, 0.8, alpha, n_min = 6), 6)
  out <- list()
  for (plan in c("oracle", "plugin")) {
    n <- if (plan == "oracle") n_or else n_pl
    if (g$design == "one-sample") {
      shift <- if (fam$type != "binary") delta else 0
      est <- one_group(fam, "H1", n, N, R_reps, mo$H1, extra_shift = shift)
      theta0 <- mo$theta0
    } else {
      shiftA <- if (fam$type != "binary") delta else 0
      eA <- one_group(fam, "H1", n, N, R_reps, mo$H1, extra_shift = shiftA)
      eB <- one_group(fam, "H0", n, N, R_reps, mo$H0)
      est <- lapply(methods, function(mt) {
        a <- eA[[mt]]; b <- eB[[mt]]
        df <- if (mt == "classical") (a$v + b$v)^2 / (a$v^2 / (n - 1) + b$v^2 / (n - 1))
              else if (mt == "plugin_adj") 2 * (n - 2) else Inf
        list(th = a$th - b$th, v = a$v + b$v, df = df)
      })
      names(est) <- methods
      theta0 <- 0
    }
    rej <- vapply(methods, function(mt) {
      e <- est[[mt]]
      df <- e$df; df[is.nan(df)] <- Inf
      crit <- qt(1 - alpha / 2, df)
      z <- abs(e$th - theta0) / sqrt(e$v); z[!is.finite(z)] <- 0
      mean(z > crit)
    }, numeric(1))
    out[[plan]] <- data.frame(family = fam$family, type = fam$type,
                              rho_target = if (is.null(fam$rho)) NA else fam$rho,
                              sens = if (is.null(fam$sens)) NA else fam$sens,
                              rho_H1 = mo$H1$rho, kappa = mo$H1$kappa,
                              design = g$design, N = N, plan = plan,
                              n_oracle_plan = n_or, n_plugin_plan = n_pl, n_used = n,
                              method = methods, achieved_power = rej,
                              stringsAsFactors = FALSE)
  }
  do.call(rbind, out)
}

plan_list <- mclapply(seq_len(nrow(plan_grid)), run_plan_row, mc.cores = n_cores, mc.set.seed = TRUE)
bad <- vapply(plan_list, inherits, logical(1), "try-error")
if (any(bad)) { print(plan_list[bad][1]); stop("Some planning configurations failed") }
plan_setting23 <- do.call(rbind, plan_list)
rownames(plan_setting23) <- NULL

settings_setting23 <- list(R = R_reps, alpha = alpha, rho_grid = rho_grid, acc_grid = acc_grid,
                           delta_rule = "0.3 * sd(Y) for continuous; p1 - p0 for binary",
                           pilot_m = "200 (1000 for rare binary)", methods = methods)

rds_dir <- file.path(script_dir, "rds_cache")
dir.create(rds_dir, showWarnings = FALSE)
saveRDS(list(results_setting23 = results_setting23, plan_setting23 = plan_setting23,
             settings_setting23 = settings_setting23),
        file.path(rds_dir, "power_validation_setting23_plugin_lambda.rds"))
cat(sprintf("Setting 23 done (%.1f min)\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
