#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive Power Validation Simulations (driver)
# =============================================================================
# Usage:
#   Rscript power_validation_comprehensive.R             # run sims + figures
#   Rscript power_validation_comprehensive.R --figures   # figures only (from cache)

script_path <- ""
if (sys.nframe() >= 1) {
  script_path <- suppressWarnings(tryCatch(normalizePath(sys.frame(1)$ofile), error = function(e) ""))
}
if (!nzchar(script_path)) {
  script_dir <- suppressWarnings(normalizePath("simulation_studies", mustWork = FALSE))
} else {
  script_dir <- dirname(script_path)
}

source(file.path(script_dir, "sim_utils.R"))

# Directory for intermediate RDS results (survives crashes)
rds_dir <- file.path(script_dir, "rds_cache")
dir.create(rds_dir, showWarnings = FALSE, recursive = TRUE)

# Check for --figures flag
args <- commandArgs(trailingOnly = TRUE)
figures_only <- any(grepl("^--figures", args))

cat("=============================================================================\n")
cat("Comprehensive Power Validation Simulations\n")
cat("=============================================================================\n")
cat(sprintf("Started: %s\n", Sys.time()))
cat(sprintf("Results cache: %s\n", rds_dir))
if (figures_only) cat("MODE: figures-only (loading from RDS cache)\n")
cat("\n")

start_time <- Sys.time()

if (figures_only) {
  # ---- Load all cached results into global env ----
  rds_files <- list.files(rds_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No cached results found in ", rds_dir, ". Run sims first.")
  for (f in rds_files) {
    objs <- readRDS(f)
    for (nm in names(objs)) {
      if (!is.null(objs[[nm]])) assign(nm, objs[[nm]], envir = .GlobalEnv)
    }
    cat(sprintf("  Loaded %s (%d objects)\n", basename(f), sum(!sapply(objs, is.null))))
  }
  cat("\n")
} else {
  # ---- Run simulations and cache results ----
  set.seed(2024)

  # Helper: run a setting script and save its results to RDS
  run_and_cache <- function(setting_file, result_names) {
    rds_path <- file.path(rds_dir, paste0(tools::file_path_sans_ext(basename(setting_file)), ".rds"))
    source(setting_file)
    results <- mget(result_names, envir = .GlobalEnv, ifnotfound = list(NULL))
    saveRDS(results, rds_path)
    cat(sprintf("  -> Cached %d objects to %s\n", sum(!sapply(results, is.null)), basename(rds_path)))
  }

  run_and_cache(file.path(script_dir, "power_validation_setting1_mean_continuous.R"),
                c("results_cont_mean", "settings_cont_mean"))
  run_and_cache(file.path(script_dir, "power_validation_setting2_mean_binary.R"),
                c("results_bin_mean"))
  run_and_cache(file.path(script_dir, "power_validation_setting3_ttest_continuous.R"),
                c("results_cont_ttest"))
  run_and_cache(file.path(script_dir, "power_validation_setting4_ttest_binary.R"),
                c("results_bin_ttest"))
  run_and_cache(file.path(script_dir, "power_validation_setting5_paired_continuous.R"),
                c("results_paired"))
  run_and_cache(file.path(script_dir, "power_validation_setting6_paired_binary.R"),
                c("results_paired_bin"))
  run_and_cache(file.path(script_dir, "power_validation_setting7_lognormal.R"),
                c("results_lognorm"))
  run_and_cache(file.path(script_dir, "power_validation_setting8_tdist.R"),
                c("results_tdist"))
  run_and_cache(file.path(script_dir, "power_validation_setting10_n_inversion.R"),
                c("results_n_inversion"))
  run_and_cache(file.path(script_dir, "power_validation_setting11_rule_of_thumb.R"),
                c("results_rule_of_thumb"))
  run_and_cache(file.path(script_dir, "power_validation_setting12_type1_error.R"),
                c("results_type1_error"))
  run_and_cache(file.path(script_dir, "power_validation_setting13_lambda_convergence.R"),
                c("results_lambda_convergence"))
  run_and_cache(file.path(script_dir, "power_validation_setting14_plugin_lambda.R"),
                c("results_plugin_lambda"))
  run_and_cache(file.path(script_dir, "power_validation_setting15_power_vs_delta.R"),
                c("results_power_vs_delta"))
  run_and_cache(file.path(script_dir, "power_validation_setting16_small_n.R"),
                c("results_small_n"))
  run_and_cache(file.path(script_dir, "power_validation_setting17_Nn_ratio.R"),
                c("results_Nn_ratio"))
  run_and_cache(file.path(script_dir, "power_validation_setting18_unequal_groups.R"),
                c("results_unequal"))
  run_and_cache(file.path(script_dir, "power_validation_setting19_misspecified_rho.R"),
                c("results_setting19"))
  run_and_cache(file.path(script_dir, "power_validation_setting20_ols_regression.R"),
                c("results_ols_regression"))
  run_and_cache(file.path(script_dir, "power_validation_setting21_glm_logistic.R"),
                c("results_glm_logistic"))
  run_and_cache(file.path(script_dir, "power_validation_setting22_2x2_table.R"),
                c("results_2x2_or", "results_2x2_rr"))
}

source(file.path(script_dir, "power_validation_plots.R"))
source(file.path(script_dir, "power_validation_summary.R"))

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat("=============================================================================\n")
cat(sprintf("Total runtime: %.1f minutes\n", total_time))
cat(sprintf("Completed: %s\n", Sys.time()))
cat("=============================================================================\n")
