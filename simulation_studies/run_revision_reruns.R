#!/usr/bin/env Rscript
# =============================================================================
# Rerun selected validation settings for the R1 revision
# =============================================================================
# Usage: Rscript simulation_studies/run_revision_reruns.R 1 3 5 7
# Each setting is run with its own seed and cached to rds_cache/ in the same
# format as power_validation_comprehensive.R, so the figure script
# (power_validation_comprehensive.R --figures) picks up the new results.

script_dir <- "simulation_studies"
source(file.path(script_dir, "sim_utils.R"))
rds_dir <- file.path(script_dir, "rds_cache")

settings <- list(
  "1"  = list(file = "power_validation_setting1_mean_continuous.R",  objs = c("results_cont_mean", "settings_cont_mean")),
  "3"  = list(file = "power_validation_setting3_ttest_continuous.R", objs = c("results_cont_ttest")),
  "4"  = list(file = "power_validation_setting4_ttest_binary.R",     objs = c("results_bin_ttest")),
  "5"  = list(file = "power_validation_setting5_paired_continuous.R", objs = c("results_paired")),
  "6"  = list(file = "power_validation_setting6_paired_binary.R",    objs = c("results_paired_bin")),
  "7"  = list(file = "power_validation_setting7_lognormal.R",        objs = c("results_lognorm")),
  "8"  = list(file = "power_validation_setting8_tdist.R",            objs = c("results_tdist")),
  "14" = list(file = "power_validation_setting14_plugin_lambda.R",   objs = c("results_plugin_lambda", "settings_plugin_lambda")),
  "20" = list(file = "power_validation_setting20_ols_regression.R",  objs = c("results_ols_regression")),
  "21" = list(file = "power_validation_setting21_glm_logistic.R",    objs = c("results_glm_logistic"))
)

ids <- commandArgs(trailingOnly = TRUE)
for (id in ids) {
  s <- settings[[id]]
  if (is.null(s)) stop("Unknown setting ", id)
  set.seed(2024 + as.integer(id))
  t0 <- Sys.time()
  source(file.path(script_dir, s$file))
  res <- mget(s$objs, envir = .GlobalEnv, ifnotfound = list(NULL))
  saveRDS(res, file.path(rds_dir, paste0(tools::file_path_sans_ext(s$file), ".rds")))
  cat(sprintf("Setting %s cached (%.1f min)\n", id, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
