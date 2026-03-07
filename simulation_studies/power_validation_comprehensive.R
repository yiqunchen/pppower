#!/usr/bin/env Rscript
# =============================================================================
# Comprehensive Power Validation Simulations (driver)
# =============================================================================

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

set.seed(2024)

cat("=============================================================================\n")
cat("Comprehensive Power Validation Simulations\n")
cat("=============================================================================\n")
cat(sprintf("Started: %s\n\n", Sys.time()))

start_time <- Sys.time()

source(file.path(script_dir, "power_validation_setting1_mean_continuous.R"))
source(file.path(script_dir, "power_validation_setting2_mean_binary.R"))
source(file.path(script_dir, "power_validation_setting3_ttest_continuous.R"))
source(file.path(script_dir, "power_validation_setting4_ttest_binary.R"))
source(file.path(script_dir, "power_validation_setting5_paired_continuous.R"))
source(file.path(script_dir, "power_validation_setting6_paired_binary.R"))
source(file.path(script_dir, "power_validation_setting7_lognormal.R"))
source(file.path(script_dir, "power_validation_setting8_tdist.R"))
source(file.path(script_dir, "power_validation_setting9_eif_binary.R"))
source(file.path(script_dir, "power_validation_setting10_n_inversion.R"))
source(file.path(script_dir, "power_validation_setting11_rule_of_thumb.R"))
source(file.path(script_dir, "power_validation_setting12_type1_error.R"))
source(file.path(script_dir, "power_validation_setting13_lambda_convergence.R"))

source(file.path(script_dir, "power_validation_plots.R"))
source(file.path(script_dir, "power_validation_summary.R"))

total_time <- difftime(Sys.time(), start_time, units = "mins")
cat("=============================================================================\n")
cat(sprintf("Total runtime: %.1f minutes\n", total_time))
cat(sprintf("Completed: %s\n", Sys.time()))
cat("=============================================================================\n")
