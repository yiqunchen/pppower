#!/usr/bin/env Rscript
# Three follow-up checks for the R1 experiment ledger.
# Run from the repository root: Rscript simulation_studies/run_additional_experiments.R
# Outputs: simulation_studies/additional_experiments/*.csv

script_dir <- if (dir.exists("simulation_studies")) "simulation_studies" else "."
source(file.path(script_dir, "revision_utils.R"))
out_dir <- file.path(script_dir, "additional_experiments")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

R_reps <- as.integer(Sys.getenv("ADDITIONAL_R", "5000"))
if (!is.finite(R_reps) || R_reps < 100) stop("ADDITIONAL_R must be at least 100")
alpha <- 0.05
methods <- c("oracle", "plugin", "plugin_adj", "pilot", "classical")
set.seed(20261012)

group_est <- function(dist_study, dist_pilot, n, N, R, shift = 0, m_pilot = 200,
                      oracle_rho = dist_study$rho) {
  d <- gen_pairs(dist_study, n, R, shift)
  fU <- gen_fU_bar(dist_study, N, R, shift)
  lam_or <- oracle_rho * N / (n + N)
  lam_pilot <- pilot_lambda(dist_pilot, m_pilot, R, n, N)
  ppi_methods(d$Y, d$f, fU, N, lam_or, lam_pilot)
}

combine_est <- function(a, b, n) {
  out <- lapply(methods, function(mt) {
    ea <- a[[mt]]; eb <- b[[mt]]
    df <- if (mt == "classical") {
      (ea$v + eb$v)^2 / (ea$v^2 / (n - 1) + eb$v^2 / (n - 1))
    } else if (mt == "plugin_adj") 2 * (n - 2) else Inf
    list(th = ea$th - eb$th, v = ea$v + eb$v, df = df)
  })
  names(out) <- methods
  out
}

summarize_est <- function(est, theta0, theta_true, experiment, design, hyp, ...) {
  do.call(rbind, lapply(methods, function(mt) {
    z <- summarize_method(est[[mt]], theta0, theta_true, alpha)
    data.frame(experiment = experiment, design = design, hyp = hyp,
               method = mt, z, R = R_reps, ..., check.names = FALSE)
  }))
}

cat(sprintf("Additional experiments: %d replicates per cell\n", R_reps))

# A. Test the low-R^2 crossover at fixed labeled sample sizes. Report
# rejection under both hypotheses so power is interpreted alongside size.
cross <- list(); k <- 0L
for (design in c("one-sample", "two-sample")) {
  for (n in c(15L, 30L, 60L)) {
    for (rho in c(0.05, 0.10, 0.20, 0.30)) {
      dist <- list(type = "gaussian", rho = rho)
      N <- 1000L; delta <- 0.3
      for (hyp in c("H0", "H1")) {
        shift <- if (hyp == "H1") delta else 0
        a <- group_est(dist, dist, n, N, R_reps, shift)
        est <- if (design == "one-sample") a else {
          b <- group_est(dist, dist, n, N, R_reps)
          combine_est(a, b, n)
        }
        k <- k + 1L
        cross[[k]] <- summarize_est(est, 0, shift, "low_r2_crossover",
                                   design, hyp, n = n, N = N, rho = rho,
                                   r2 = rho^2, threshold_r2 = 1 / (n + 1),
                                   delta = delta)
      }
    }
  }
}
cross <- do.call(rbind, cross)
write.csv(cross, file.path(out_dir, "low_r2_crossover.csv"), row.names = FALSE)
cat("Low-R2 crossover complete\n")

# B. Pilot-to-study shift: plan n using the pilot population correlation,
# then draw an independent pilot for lambda and analyze study-population
# observations. H0 and H1 share the same n for each grid cell.
shift_results <- list(); k <- 0L
for (design in c("one-sample", "two-sample")) {
  for (rho_pilot in c(0.3, 0.6, 0.9)) {
    for (rho_study in c(0.3, 0.6, 0.9)) {
      N <- 500L; delta <- 0.2
      V_plan <- function(n) {
        v <- 1 / n - rho_pilot^2 * N / (n * (n + N))
        if (design == "two-sample") 2 * v else v
      }
      n <- plan_n(V_plan, delta, 0.8, alpha, n_min = 8)
      dist_pilot <- list(type = "gaussian", rho = rho_pilot)
      dist_study <- list(type = "gaussian", rho = rho_study)
      for (hyp in c("H0", "H1")) {
        effect <- if (hyp == "H1") delta else 0
        a <- group_est(dist_study, dist_pilot, n, N, R_reps, effect)
        est <- if (design == "one-sample") a else {
          b <- group_est(dist_study, dist_pilot, n, N, R_reps)
          combine_est(a, b, n)
        }
        k <- k + 1L
        shift_results[[k]] <- summarize_est(est, 0, effect, "pilot_study_shift",
                                            design, hyp, n = n, N = N,
                                            rho_pilot = rho_pilot,
                                            rho_study = rho_study,
                                            m_pilot = 200L, delta = delta)
      }
    }
  }
}
shift_results <- do.call(rbind, shift_results)
write.csv(shift_results, file.path(out_dir, "pilot_study_shift.csv"), row.names = FALSE)
cat("Pilot-to-study shift complete\n")

# C. Type I error versus n and skewness. The strongly skewed log-normal
# cells are compared with a milder log-normal and a Gaussian benchmark.
skew_results <- list(); k <- 0L
for (design in c("one-sample", "two-sample")) {
  for (sigma_log in c(0, 0.5, 1.0)) {
    dist <- if (sigma_log == 0) list(type = "gaussian", rho = 0.5) else
      list(type = "lognormal", rho = 0.5, sigma_log = sigma_log)
    mom <- pop_moments(dist, big = 500000)
    for (n in c(15L, 20L, 30L, 50L, 100L, 200L)) {
      N <- 1000L
      d <- gen_pairs(dist, n, R_reps)
      fU <- gen_fU_bar(dist, N, R_reps)
      lam_or <- mom$cov / ((1 + n / N) * mom$sf2)
      lam_pilot <- pilot_lambda(dist, 200L, R_reps, n, N)
      a <- ppi_methods(d$Y, d$f, fU, N, lam_or, lam_pilot)
      est <- if (design == "one-sample") a else {
        db <- gen_pairs(dist, n, R_reps)
        fUb <- gen_fU_bar(dist, N, R_reps)
        lam_pb <- pilot_lambda(dist, 200L, R_reps, n, N)
        b <- ppi_methods(db$Y, db$f, fUb, N, lam_or, lam_pb)
        combine_est(a, b, n)
      }
      theta0 <- if (design == "one-sample") mom$theta else 0
      k <- k + 1L
      skew_results[[k]] <- summarize_est(est, theta0, theta0,
                                         "skew_calibration", design, "H0",
                                         n = n, N = N, rho = 0.5,
                                         sigma_log = sigma_log, m_pilot = 200L)
    }
  }
}
skew_results <- do.call(rbind, skew_results)
write.csv(skew_results, file.path(out_dir, "skew_calibration.csv"), row.names = FALSE)
cat("Strong-skew calibration complete\n")

saveRDS(list(low_r2_crossover = cross, pilot_study_shift = shift_results,
             skew_calibration = skew_results,
             settings = list(R = R_reps, alpha = alpha, seed = 20261012)),
        file.path(out_dir, "additional_experiments.rds"))
cat("All additional experiments complete\n")
