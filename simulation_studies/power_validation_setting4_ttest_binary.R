# =============================================================================
# SETTING 4: TWO-SAMPLE PROPORTION TEST - BINARY OUTCOME
# =============================================================================
cat("Setting 4: Two-Sample Proportion Test - Binary Outcome\n")
cat("-------------------------------------------------------\n")
setting4_start <- Sys.time()

# Locate utils if not already sourced
if (!exists("simulate_twosample_binary_power")) {
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
}

# Parameters
n_values_prop <- c(20, 40, 60, 80, 100)
N_prop_values <- c(200, 500)
p0_prop <- 0.3  # Base prevalence
delta_prop <- 0.08  # Difference in proportions
classifier_settings <- list(
  c(sens = 0.70, spec = 0.70),
  c(sens = 0.85, spec = 0.85),
  c(sens = 0.95, spec = 0.95)
)
alpha <- 0.05
R <- 1000  # Monte Carlo replications

results_bin_ttest <- do.call(
  rbind,
  lapply(classifier_settings, function(classifier) {
    sens <- classifier["sens"]
    spec <- classifier["spec"]
    simulate_twosample_binary_power(
      p0 = p0_prop,
      delta = delta_prop,
      sens = sens,
      spec = spec,
      n_values = n_values_prop,
      N_values = N_prop_values,
      R = R,
      alpha = alpha,
      seed = NULL,
      progress = TRUE
    )
  })
)

cat(sprintf("Done. (%.1f seconds)\n\n", difftime(Sys.time(), setting4_start, units = "secs")))
