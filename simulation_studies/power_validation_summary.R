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

if (exists("results_eif_binary") && is.data.frame(results_eif_binary)) {
  cat("Setting 9: EIF Binary - Max |Theoretical - Empirical|:\n")
  max_diff_9 <- max(abs(results_eif_binary$power_theo - results_eif_binary$power_emp), na.rm = TRUE)
  cat(sprintf("  %.4f\n\n", max_diff_9))
}

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
  cat("Setting 13: Plugin Lambda Convergence (at n=2000):\n")
  subset_2k <- results_lambda_convergence[results_lambda_convergence$n == 2000, ]
  cat(sprintf("  Max RMSE: %.4f\n", max(subset_2k$rmse)))
  cat(sprintf("  Max |bias|: %.4f\n\n", max(abs(subset_2k$bias))))
}
