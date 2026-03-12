# =============================================================================
# GENERATE VALIDATION FIGURES
# =============================================================================
# Purpose: Validate that our analytical (theoretical) PPI++ power formulas
#          match Monte Carlo (empirical) power estimates.
# Style:   Each panel shows analytical lines + empirical points for PPI++.
# =============================================================================
cat("Generating validation figures...\n")
fig_start <- Sys.time()

library(ggplot2)
library(dplyr)
library(tidyr)

# Palette objects are exported by pppower; fallback for dev mode.
if (!exists("pppower_colors")) {
  tryCatch(source("R/palette.R", local = FALSE), error = function(e) {
    tryCatch(source("../R/palette.R", local = FALSE), error = function(e2) NULL)
  })
}

theme_validation <- theme_pppower

# Colors: PPI++ analytical/empirical + classical reference (mid-grey)
colors_val <- c(
  "Analytical" = pppower_colors$analytical,
  "Empirical"  = pppower_colors$empirical,
  "Classical"  = "gray55"
)

# =============================================================================
# Extract settings variables from cached results (for --figures mode)
# =============================================================================
if (exists("results_cont_mean") && is.data.frame(results_cont_mean)) {
  delta      <- results_cont_mean$delta[1]
  N_values   <- sort(unique(results_cont_mean$N))
  R          <- 1000
  if (!exists("settings_cont_mean")) {
    settings_cont_mean <- list(
      n_values = sort(unique(results_cont_mean$n)),
      N_values = N_values,
      rho_values = sort(unique(results_cont_mean$rho)),
      sigma_Y2 = 1, sigma_f2 = 1, alpha = 0.05, R = R
    )
  }
}
if (exists("results_bin_mean") && is.data.frame(results_bin_mean)) {
  p0 <- 0.3; delta_bin <- results_bin_mean$delta[1]
  N_bin_values <- sort(unique(results_bin_mean$N))
}
if (exists("results_cont_ttest") && is.data.frame(results_cont_ttest)) {
  delta_ttest <- results_cont_ttest$delta[1]
  N_ttest_values <- sort(unique(results_cont_ttest$N))
}
if (exists("results_bin_ttest") && is.data.frame(results_bin_ttest)) {
  delta_prop <- results_bin_ttest$delta[1]; p0_prop <- 0.3
  N_prop_values <- sort(unique(results_bin_ttest$N))
}
if (exists("results_paired") && is.data.frame(results_paired)) {
  delta_paired <- results_paired$delta[1]
  N_paired_values <- sort(unique(results_paired$N))
}
if (exists("results_paired_bin") && is.data.frame(results_paired_bin)) {
  delta_paired_bin <- results_paired_bin$delta[1]; p0_paired <- 0.3
  N_paired_bin_values <- sort(unique(results_paired_bin$N))
}
if (exists("results_lognorm") && is.data.frame(results_lognorm)) {
  delta_lognorm <- results_lognorm$delta[1]
  N_lognorm_values <- sort(unique(results_lognorm$N))
}
if (exists("results_tdist") && is.data.frame(results_tdist)) {
  delta_t <- results_tdist$delta[1]
  N_t_values <- sort(unique(results_tdist$N))
}

# Output directory
if (exists("script_dir") && nzchar(script_dir)) {
  output_dir <- script_dir
} else {
  output_dir <- "simulation_studies"
}
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("  Output directory: %s\n", normalizePath(output_dir, mustWork = FALSE)))

save_fig <- function(plot, name, width = 10, height = 4.5) {
  base <- file.path(output_dir, name)
  ggsave(paste0(base, ".pdf"), plot, width = width, height = height)
  ggsave(paste0(base, ".png"), plot, width = width, height = height, dpi = 400, bg = "white")
  cat(sprintf("    Saved %s\n", name))
}

# =============================================================================
# Panel A: Mean (Continuous)
# =============================================================================
plot_data_A <- results_cont_mean %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_A <- ggplot(plot_data_A, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A. Mean Estimation (Continuous)",
    subtitle = sprintf("PPI++ power: Delta = %.1f, R = %d replications",
                       delta, R),
    x = "Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel B: Mean (Binary)
# =============================================================================
plot_data_B <- results_bin_mean %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_B <- ggplot(plot_data_B, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "B. Mean Estimation (Binary/Prevalence)",
    subtitle = sprintf("PPI++ power: Delta = %.2f, p0 = %.1f",
                       delta_bin, p0),
    x = "Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel C: t-Test (Continuous)
# =============================================================================
plot_data_C <- results_cont_ttest %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_C <- ggplot(plot_data_C, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "C. Two-Sample t-Test (Continuous)",
    subtitle = sprintf("PPI++ power: Delta = %.1f", delta_ttest),
    x = "Per-Group Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel D: Proportion Test (Binary)
# =============================================================================
plot_data_D <- results_bin_ttest %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_D <- ggplot(plot_data_D, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "D. Two-Sample Proportion Test (Binary)",
    subtitle = sprintf("PPI++ power: Delta = %.2f, p0 = %.1f",
                       delta_prop, p0_prop),
    x = "Per-Group Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel E: Paired t-Test (Continuous)
# =============================================================================
plot_data_E <- results_paired %>%
  mutate(rho_label = paste0("rho[D] == ", rho_D),
         N_label = paste0("N == ", N))

fig_E <- ggplot(plot_data_E, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "E. Paired t-Test (Continuous)",
    subtitle = sprintf("PPI++ power: Delta = %.1f", delta_paired),
    x = "Number of Pairs (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel F: Paired Proportion Test (Binary)
# =============================================================================
plot_data_F <- results_paired_bin %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_F <- ggplot(plot_data_F, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "F. Paired Proportion Test (Binary)",
    subtitle = sprintf("PPI++ power: Delta = %.2f, p0 = %.1f",
                       delta_paired_bin, p0_paired),
    x = "Number of Pairs (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel G: Log-Normal (Skewed)
# =============================================================================
plot_data_G <- results_lognorm %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_G <- ggplot(plot_data_G, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "G. Mean Estimation (Log-Normal, Skewed)",
    subtitle = sprintf("PPI++ power: Delta = %.3f", delta_lognorm),
    x = "Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel H: t-distribution (Heavy-Tailed)
# =============================================================================
plot_data_H <- results_tdist %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_H <- ggplot(plot_data_H, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "H. Mean Estimation (t-dist df=5, Heavy-Tailed)",
    subtitle = sprintf("PPI++ power: Delta = %.1f", delta_t),
    x = "Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Save panels A-H
# =============================================================================
save_fig(fig_A, "fig_validation_A_mean_continuous")
save_fig(fig_B, "fig_validation_B_mean_binary")
save_fig(fig_C, "fig_validation_C_ttest_continuous")
save_fig(fig_D, "fig_validation_D_ttest_binary")
save_fig(fig_E, "fig_validation_E_paired_continuous")
save_fig(fig_F, "fig_validation_F_paired_binary")
save_fig(fig_G, "fig_validation_G_lognormal")
save_fig(fig_H, "fig_validation_H_tdist")

# =============================================================================
# Panel I: EIF Binary Validation (Setting 9)
# =============================================================================
if (exists("results_eif_binary") && is.data.frame(results_eif_binary)) {
  plot_data_I <- results_eif_binary
  plot_data_I$clf_label <- sprintf("Sens/Spec = %.0f%%", plot_data_I$sens * 100)
  plot_data_I$p_label <- paste0("p = ", plot_data_I$p)

  fig_I <- ggplot(plot_data_I, aes(x = m_cal)) +
    geom_line(aes(y = power_theo, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp, color = "Empirical"), size = 3.5) +
    facet_grid(p_label ~ clf_label) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "I. EIF Binary Surrogate Validation",
      subtitle = sprintf("PPI++ power: Delta = %.2f, R = %d", 0.05, 1000),
      x = "Calibration Sample Size (m_cal)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_I, "fig_validation_I_eif_binary", height = 5)
}

# =============================================================================
# Panel J: Sample Size Inversion (Setting 10)
# =============================================================================
if (exists("results_n_inversion") && is.data.frame(results_n_inversion)) {
  plot_data_J <- results_n_inversion[!is.na(results_n_inversion$analytical_achieved), ]

  fig_J <- ggplot(plot_data_J, aes(x = target_power, y = analytical_achieved)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = pppower_colors$reference) +
    geom_point(aes(color = design, shape = design), size = 3.5) +
    geom_point(aes(x = target_power, y = empirical_achieved, color = design),
               shape = 1, size = 3,
               data = plot_data_J[!is.na(plot_data_J$empirical_achieved), ]) +
    scale_color_manual(values = c(
      mean_continuous   = pppower_colors$ppi_pp,
      eif_binary        = pppower_colors$oracle,
      paired_continuous = pppower_colors$vanilla
    )) +
    scale_x_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
    scale_y_continuous(limits = c(0.5, 1), breaks = seq(0.5, 1, 0.1)) +
    labs(
      title = "J. Sample Size Inversion Validation",
      subtitle = "Achieved power at n* vs. target power (dashed = identity)",
      x = "Target Power", y = "Achieved Power at n*",
      color = "Design", shape = "Design"
    ) +
    theme_validation()

  save_fig(fig_J, "fig_validation_J_n_inversion", width = 7, height = 5)
}

# =============================================================================
# Panel K: Rule-of-Thumb (Setting 11)
# =============================================================================
if (exists("results_rule_of_thumb") && is.data.frame(results_rule_of_thumb)) {
  plot_data_K <- results_rule_of_thumb
  plot_data_K$N_label <- paste0("N = ", formatC(plot_data_K$N, big.mark = ","))

  fig_K <- ggplot(plot_data_K, aes(x = rho_sq)) +
    geom_line(aes(y = ratio, color = N_label), linewidth = 1.8) +
    geom_line(aes(y = theoretical_ratio), linetype = "dashed",
              color = pppower_colors$reference, linewidth = 1.2) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(
      title = expression(paste("K. Rule of Thumb: ", n[PPI] / n[classical], " vs. 1 - ", rho^2)),
      subtitle = "Dashed line = theoretical 1 - rho^2",
      x = expression(rho^2),
      y = expression(n[PPI] / n[classical]),
      color = "Unlabeled N"
    ) +
    theme_validation()

  save_fig(fig_K, "fig_validation_K_rule_of_thumb", width = 7, height = 5)
}

# =============================================================================
# Panel L: Type I Error (Setting 12)
# =============================================================================
if (exists("results_type1_error") && is.list(results_type1_error)) {
  ci_lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / 2000)
  ci_hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / 2000)

  if (!is.null(results_type1_error$continuous)) {
    plot_data_L <- results_type1_error$continuous
    plot_data_L$rho_label <- paste0("rho = ", plot_data_L$rho)
    plot_data_L$N_label <- paste0("N = ", plot_data_L$N)

    fig_L <- ggplot(plot_data_L, aes(x = n)) +
      geom_point(aes(y = rejection_ppi), color = pppower_colors$analytical, size = 3.5) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
      geom_hline(yintercept = ci_lo, linetype = "dotted", color = pppower_colors$reference) +
      geom_hline(yintercept = ci_hi, linetype = "dotted", color = pppower_colors$reference) +
      facet_grid(N_label ~ rho_label) +
      scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.02)) +
      labs(
        title = "L. Type I Error Calibration (PPI++, Gaussian Mean, H0)",
        subtitle = "Dashed = nominal 0.05; dotted = 95% MC CI",
        x = "Labeled Sample Size (n)", y = "Rejection Rate"
      ) +
      theme_validation()

    save_fig(fig_L, "fig_validation_L_type1_continuous", height = 5)
  }
}

# =============================================================================
# Panel M: Lambda Convergence (Setting 13)
# =============================================================================
if (exists("results_lambda_convergence") && is.data.frame(results_lambda_convergence)) {
  plot_data_M <- results_lambda_convergence
  plot_data_M$rho_label <- paste0("rho == ", plot_data_M$rho)

  fig_M <- ggplot(plot_data_M, aes(x = n)) +
    geom_ribbon(aes(ymin = lambda_hat_mean - lambda_hat_sd,
                    ymax = lambda_hat_mean + lambda_hat_sd),
                alpha = 0.2, fill = pppower_colors$ppi_pp) +
    geom_line(aes(y = lambda_hat_mean, color = "Plugin mean"), linewidth = 1.8) +
    geom_hline(aes(yintercept = lambda_oracle, color = "Oracle"),
               linetype = "dashed", linewidth = 0.8,
               data = plot_data_M[!duplicated(plot_data_M$rho), ]) +
    facet_wrap(~rho_label, scales = "free_y",
               labeller = labeller(rho_label = label_parsed)) +
    scale_color_manual(values = c(
      "Plugin mean" = pppower_colors$ppi_pp,
      "Oracle"      = pppower_colors$classical
    )) +
    scale_x_log10() +
    labs(
      title = "M. Plugin Lambda Convergence to Oracle",
      subtitle = "Ribbon = +/- 1 SD across R = 1000 replications",
      x = "Labeled Sample Size (n, log scale)",
      y = expression(lambda), color = ""
    ) +
    theme_validation()

  save_fig(fig_M, "fig_validation_M_lambda_convergence")
}

# =============================================================================
# Panel N: Plugin Lambda Power (Setting 14)
# =============================================================================
if (exists("results_plugin_lambda") && is.data.frame(results_plugin_lambda)) {
  plot_data_N <- results_plugin_lambda %>%
    mutate(rho_label = paste0("rho == ", rho),
           N_label = paste0("N == ", N))

  # Reshape to long format: oracle vs plugin empirical
  plot_data_N_long <- plot_data_N %>%
    tidyr::pivot_longer(
      cols = c(power_emp_oracle, power_emp_plugin),
      names_to = "lambda_type",
      values_to = "power_emp",
      names_prefix = "power_emp_"
    ) %>%
    mutate(lambda_label = ifelse(lambda_type == "oracle",
                                  "Oracle lambda*", "Plugin lambda"))

  fig_N <- ggplot(plot_data_N, aes(x = n)) +
    geom_line(aes(y = power_theo, color = "Analytical"), linewidth = 1.8) +
    geom_point(data = plot_data_N_long %>% filter(lambda_type == "oracle"),
               aes(y = power_emp, color = "Empirical (Oracle)"),
               size = 3.5, shape = 16) +
    geom_point(data = plot_data_N_long %>% filter(lambda_type == "plugin"),
               aes(y = power_emp, color = "Empirical (Plugin)"),
               size = 3.5, shape = 17) +
    facet_grid(N_label ~ rho_label,
               labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = c(
      "Analytical"        = pppower_colors$analytical,
      "Empirical (Oracle)" = pppower_colors$empirical,
      "Empirical (Plugin)" = pppower_colors$vanilla
    )) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "N. Plugin vs Oracle Lambda Power",
      subtitle = "Analytical power formula vs MC with oracle/plugin lambda",
      x = "Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_N, "fig_validation_N_plugin_lambda")
}

# =============================================================================
# Panel O: Power vs Effect Size (Setting 15)
# =============================================================================
if (exists("results_power_vs_delta") && is.data.frame(results_power_vs_delta)) {
  plot_data_O <- results_power_vs_delta %>%
    mutate(rho_label = paste0("rho == ", rho),
           n_label = paste0("n == ", n))

  fig_O <- ggplot(plot_data_O, aes(x = delta)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 2.5) +
    facet_grid(n_label ~ rho_label,
               labeller = labeller(rho_label = label_parsed, n_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "O. Power vs Effect Size (Delta)",
      subtitle = sprintf("PPI++ power: N = %d, varying delta",
                         results_power_vs_delta$N[1]),
      x = expression(Delta), y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_O, "fig_validation_O_power_vs_delta", height = 6)
}

# =============================================================================
# Panel P: Small n Regime (Setting 16)
# =============================================================================
if (exists("results_small_n") && is.data.frame(results_small_n)) {
  plot_data_P <- results_small_n %>%
    mutate(rho_label = paste0("rho == ", rho),
           N_label = paste0("N == ", N))

  fig_P <- ggplot(plot_data_P, aes(x = n)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_grid(N_label ~ rho_label,
               labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 30, linetype = "dotted", color = "gray70", linewidth = 0.5) +
    annotate("text", x = 32, y = 0.05, label = "n=30", hjust = 0,
             size = 3, color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "P. Small n Regime (CLT Breakdown)",
      subtitle = sprintf("PPI++ power: Delta = %.1f, R = 2000 reps",
                         results_small_n$delta[1]),
      x = "Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_P, "fig_validation_P_small_n")
}

# =============================================================================
# Panel Q: N/n Ratio Sensitivity (Setting 17)
# =============================================================================
if (exists("results_Nn_ratio") && is.data.frame(results_Nn_ratio)) {
  plot_data_Q <- results_Nn_ratio %>%
    mutate(rho_label = paste0("rho == ", rho))

  fig_Q <- ggplot(plot_data_Q, aes(x = Nn_ratio)) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical (PPI++)"),
              linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical (PPI++)"),
               size = 3.5) +
    geom_hline(aes(yintercept = power_theo_classical),
               linetype = "dashed", color = "gray55",
               linewidth = 0.8,
               data = plot_data_Q[!duplicated(plot_data_Q$rho), ]) +
    facet_wrap(~rho_label, labeller = labeller(rho_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 100)) +
    scale_color_manual(values = c(
      "Analytical (PPI++)" = pppower_colors$analytical,
      "Empirical (PPI++)"  = pppower_colors$empirical
    )) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Q. N/n Ratio Sensitivity",
      subtitle = sprintf("n = %d fixed; dashed = classical power (no PPI)",
                         results_Nn_ratio$n[1]),
      x = "N / n Ratio (log scale)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_Q, "fig_validation_Q_Nn_ratio", width = 10, height = 4)
}

# =============================================================================
# Panel R: Unequal Group Sizes (Setting 18)
# =============================================================================
if (exists("results_unequal") && is.data.frame(results_unequal)) {
  plot_data_R <- results_unequal %>%
    mutate(rho_label = paste0("rho == ", rho),
           alloc_label = sprintf("n_A:n_B = %d:%d", n_A, n_B))

  fig_R <- ggplot(plot_data_R, aes(x = alloc_ratio)) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_wrap(~rho_label, labeller = labeller(rho_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    scale_x_continuous(breaks = sort(unique(results_unequal$alloc_ratio)),
                       labels = sprintf("%.0f:1", sort(unique(results_unequal$alloc_ratio)))) +
    labs(
      title = "R. Unequal Group Sizes (Two-Sample t-Test)",
      subtitle = sprintf("n_total = %d, N_total = %d, Delta = %.1f",
                         sum(results_unequal$n_A[1], results_unequal$n_B[1]),
                         sum(results_unequal$N_A[1], results_unequal$N_B[1]),
                         results_unequal$delta[1]),
      x = expression(n[A]:n[B] ~ "Allocation Ratio"), y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_R, "fig_validation_R_unequal_groups", width = 10, height = 4)
}

# =============================================================================
# Panel S: OLS Regression Contrast (Setting 20)
# =============================================================================
if (exists("results_ols_regression") && is.data.frame(results_ols_regression)) {
  plot_data_S <- results_ols_regression %>%
    mutate(rho_label = paste0("rho == ", rho),
           N_label   = paste0("N == ", N))

  fig_S <- ggplot(plot_data_S, aes(x = n)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_grid(N_label ~ rho_label,
               labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "S. OLS Regression Contrast",
      subtitle = sprintf("PPI++ power: a = (1,-1), Delta = %.1f, R = %d",
                         results_ols_regression$delta[1], R),
      x = "Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_S, "fig_validation_S_ols_regression")
}

# =============================================================================
# Panel T: GLM Logistic Regression Contrast (Setting 21)
# =============================================================================
if (exists("results_glm_logistic") && is.data.frame(results_glm_logistic)) {
  plot_data_T <- results_glm_logistic %>%
    mutate(acc_label = sprintf("Accuracy = %.0f%%", accuracy * 100),
           N_label   = paste0("N == ", N))

  fig_T <- ggplot(plot_data_T, aes(x = n)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_grid(N_label ~ acc_label,
               labeller = labeller(N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "T. GLM Logistic Regression Contrast",
      subtitle = sprintf("PPI++ power: a = (1,-1), Delta = %.1f, R = %d",
                         results_glm_logistic$delta[1], R),
      x = "Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_T, "fig_validation_T_glm_logistic")
}

# =============================================================================
# Panel U: 2x2 Table OR (Setting 22)
# =============================================================================
if (exists("results_2x2_or") && is.data.frame(results_2x2_or)) {
  plot_data_U <- results_2x2_or %>%
    mutate(p_exp_label = sprintf("p[exp] == %.2f", p_exp),
           acc_label   = sprintf("Accuracy = %.0f%%", accuracy * 100),
           N_label     = paste0("N == ", N))

  fig_U <- ggplot(plot_data_U, aes(x = n)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_grid(N_label ~ acc_label,
               labeller = labeller(N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "U. 2x2 Table: Odds Ratio (PPI++ Logistic)",
      subtitle = sprintf("p_ctrl = %.2f, prev_exp = %.1f, R = %d",
                         results_2x2_or$p_ctrl[1], 0.5, 1000),
      x = "Total Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_U, "fig_validation_U_2x2_or")
}

# =============================================================================
# Panel V: 2x2 Table RR (Setting 22)
# =============================================================================
if (exists("results_2x2_rr") && is.data.frame(results_2x2_rr)) {
  plot_data_V <- results_2x2_rr %>%
    mutate(p_exp_label = sprintf("p[exp] == %.2f", p_exp),
           acc_label   = sprintf("Accuracy = %.0f%%", accuracy * 100),
           N_label     = paste0("N == ", N))

  fig_V <- ggplot(plot_data_V, aes(x = n)) +
    geom_line(aes(y = power_theo_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.8) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_grid(N_label ~ acc_label,
               labeller = labeller(N_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "V. 2x2 Table: Relative Risk (PPI++ Delta Method)",
      subtitle = sprintf("p_ctrl = %.2f, prev_exp = %.1f, R = %d",
                         results_2x2_rr$p_ctrl[1], 0.5, 1000),
      x = "Total Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_V, "fig_validation_V_2x2_rr")
}

cat(sprintf("Figures generated. (%.1f seconds)\n", difftime(Sys.time(), fig_start, units = "secs")))
cat(sprintf("Figures saved to: %s\n\n", normalizePath(output_dir, mustWork = FALSE)))
