# =============================================================================
# SUMMARY TABLES
# =============================================================================
cat("\n=============================================================================\n")
cat("VALIDATION SUMMARY\n")
cat("=============================================================================\n\n")

cat("Setting 1: Mean (Continuous) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_1 <- max(abs(results_cont_mean$power_theo_ppi - results_cont_mean$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_1))

cat("Setting 2: Mean (Binary) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_2 <- max(abs(results_bin_mean$power_theo_ppi - results_bin_mean$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_2))

cat("Setting 3: t-Test (Continuous) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_3 <- max(abs(results_cont_ttest$power_theo_ppi - results_cont_ttest$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_3))

cat("Setting 4: Proportion Test (Binary) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_4 <- max(abs(results_bin_ttest$power_theo_ppi - results_bin_ttest$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_4))

cat("Setting 5: Paired t-Test (Continuous) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_5 <- max(abs(results_paired$power_theo_ppi - results_paired$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_5))

cat("Setting 6: Paired Proportion Test (Binary) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_6 <- max(abs(results_paired_bin$power_theo_ppi - results_paired_bin$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_6))

cat("Setting 7: Mean (Log-Normal, Skewed) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_7 <- max(abs(results_lognorm$power_theo_ppi - results_lognorm$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_7))

cat("Setting 8: Mean (t-dist df=5, Heavy-Tailed) - Max |Theoretical - Empirical| for PPI++:\n")
max_diff_8 <- max(abs(results_tdist$power_theo_ppi - results_tdist$power_emp_ppi))
cat(sprintf("  %.4f\n\n", max_diff_8))

if (exists("results_n_inversion") && is.data.frame(results_n_inversion)) {
  cat("Setting 10: Sample Size Inversion - Achieved vs. Target Power:\n")
  valid_10 <- results_n_inversion[!is.na(results_n_inversion$analytical_achieved), ]
  cat(sprintf("  All analytical >= target? %s\n",
              ifelse(all(valid_10$analytical_achieved >= valid_10$target_power - 0.001), "YES", "NO")))
  cat(sprintf("  Max shortfall: %.4f\n\n",
              max(valid_10$target_power - valid_10$analytical_achieved)))
}

if (exists("results_rule_of_thumb") && is.data.frame(results_rule_of_thumb)) {
  cat("Setting 11: Rule of Thumb - Max |ratio - (1-rho^2)|:\n")
  max_diff_11 <- max(results_rule_of_thumb$abs_diff, na.rm = TRUE)
  cat(sprintf("  %.4f\n\n", max_diff_11))
}

if (exists("results_type1_error") && is.list(results_type1_error)) {
  cat("Setting 12: Type I Error Calibration:\n")
  if (!is.null(results_type1_error$continuous)) {
    cat(sprintf("  Continuous PPI++ range: [%.4f, %.4f]\n",
                min(results_type1_error$continuous$rejection_ppi),
                max(results_type1_error$continuous$rejection_ppi)))
  }
  if (!is.null(results_type1_error$binary)) {
    cat(sprintf("  Binary PPI++ range:     [%.4f, %.4f]\n\n",
                min(results_type1_error$binary$rejection_ppi),
                max(results_type1_error$binary$rejection_ppi)))
  }
}

if (exists("results_lambda_convergence") && is.data.frame(results_lambda_convergence)) {
  cat("Setting 13: Plugin Lambda Convergence (at n=500):\n")
  subset_2k <- results_lambda_convergence[results_lambda_convergence$n == 500, ]
  cat(sprintf("  Max RMSE: %.4f\n", max(subset_2k$rmse)))
  cat(sprintf("  Max |bias|: %.4f\n\n", max(abs(subset_2k$bias))))
}

if (exists("results_plugin_lambda") && is.data.frame(results_plugin_lambda)) {
  cat("Setting 14: Plugin Lambda Power - Max |theo - emp|:\n")
  cat(sprintf("  Oracle:  %.4f\n",
              max(abs(results_plugin_lambda$power_theo - results_plugin_lambda$power_emp_oracle))))
  cat(sprintf("  Plugin:  %.4f\n",
              max(abs(results_plugin_lambda$power_theo - results_plugin_lambda$power_emp_plugin))))
  cat(sprintf("  |Oracle - Plugin|: %.4f\n\n",
              max(abs(results_plugin_lambda$power_emp_oracle - results_plugin_lambda$power_emp_plugin))))
}

if (exists("results_power_vs_delta") && is.data.frame(results_power_vs_delta)) {
  cat("Setting 15: Power vs Delta - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_power_vs_delta$power_theo_ppi - results_power_vs_delta$power_emp_ppi))))
}

if (exists("results_small_n") && is.data.frame(results_small_n)) {
  cat("Setting 16: Small n Regime - Max |theo - emp| PPI++:\n")
  small <- results_small_n[results_small_n$n <= 25, ]
  large <- results_small_n[results_small_n$n >= 50, ]
  cat(sprintf("  n <= 25: %.4f\n", max(abs(small$power_theo_ppi - small$power_emp_ppi))))
  cat(sprintf("  n >= 50: %.4f\n\n", max(abs(large$power_theo_ppi - large$power_emp_ppi))))
}

if (exists("results_Nn_ratio") && is.data.frame(results_Nn_ratio)) {
  cat("Setting 17: N/n Ratio Sensitivity - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n", max(abs(results_Nn_ratio$power_theo_ppi - results_Nn_ratio$power_emp_ppi))))
  r1 <- results_Nn_ratio[results_Nn_ratio$Nn_ratio == 1, ]
  cat(sprintf("  At N/n=1, PPI++ ~ Classical (max gap: %.4f)\n\n",
              max(abs(r1$power_emp_ppi - r1$power_emp_classical))))
}

if (exists("results_unequal") && is.data.frame(results_unequal)) {
  cat("Setting 18: Unequal Groups - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_unequal$power_theo_ppi - results_unequal$power_emp_ppi))))
}

if (exists("results_setting19") && is.data.frame(results_setting19)) {
  cat("Setting 19: Misspecified Rho - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_setting19$power_theo_true - results_setting19$power_emp_true))))
}

if (exists("results_ols_regression") && is.data.frame(results_ols_regression)) {
  cat("Setting 20: OLS Regression Contrast - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_ols_regression$power_theo_ppi - results_ols_regression$power_emp_ppi))))
}

if (exists("results_glm_logistic") && is.data.frame(results_glm_logistic)) {
  cat("Setting 21: GLM Logistic Contrast - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_glm_logistic$power_theo_ppi - results_glm_logistic$power_emp_ppi))))
}

if (exists("results_2x2_or") && is.data.frame(results_2x2_or)) {
  cat("Setting 22 (OR): 2x2 Table Odds Ratio - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_2x2_or$power_theo_ppi - results_2x2_or$power_emp_ppi), na.rm = TRUE)))
}

if (exists("results_2x2_rr") && is.data.frame(results_2x2_rr)) {
  cat("Setting 22 (RR): 2x2 Table Relative Risk - Max |theo - emp| PPI++:\n")
  cat(sprintf("  %.4f\n\n",
              max(abs(results_2x2_rr$power_theo_ppi - results_2x2_rr$power_emp_ppi), na.rm = TRUE)))
}
