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
theme_validation_focus <- function(base_size = 16) {
  theme_validation(base_size = base_size) +
    theme(
      axis.text = element_text(face = "bold", size = base_size - 1),
      axis.title = element_text(face = "bold", size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size - 1)
    )
}

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
    title = "Mean Estimation (Continuous)",
    x = "Labeled Sample Size (n)", y = "Power", color = ""
  ) +
  theme_validation()

# =============================================================================
# Panel B: Mean (Binary)
# =============================================================================
plot_data_B <- results_bin_mean
if (!"mc_se_emp_ppi" %in% names(plot_data_B)) {
  mc_reps_B <- if ("mc_reps" %in% names(plot_data_B)) plot_data_B$mc_reps else 1000
  plot_data_B$mc_se_emp_ppi <- sqrt(
    pmax(plot_data_B$power_emp_ppi * (1 - plot_data_B$power_emp_ppi), 0) / mc_reps_B
  )
}

plot_data_B <- plot_data_B %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N),
         power_emp_ppi_lo = pmax(0, power_emp_ppi - 1.96 * mc_se_emp_ppi),
         power_emp_ppi_hi = pmin(1, power_emp_ppi + 1.96 * mc_se_emp_ppi))

fig_B <- ggplot(plot_data_B, aes(x = n)) +
  geom_line(aes(y = power_theo_classical, color = "Classical"),
            linetype = "dashed", linewidth = 1.8) +
  geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
  geom_errorbar(aes(ymin = power_emp_ppi_lo, ymax = power_emp_ppi_hi, color = "Empirical"),
                width = 4, linewidth = 0.7, alpha = 0.8) +
  geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_val) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "Mean Estimation (Binary)",
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
    title = "Two-Sample t-Test (Continuous)",
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
    title = "Two-Sample Proportion Test (Binary)",
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
    title = "Paired t-Test (Continuous)",
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
    title = "Paired Proportion Test (Binary)",
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
    title = "Mean Estimation (Log-Normal)",
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
    title = "Mean Estimation (t-Distribution, df = 5)",
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
# Panel J: Sample Size Inversion (Setting 10)
# =============================================================================
if (exists("results_n_inversion") && is.data.frame(results_n_inversion)) {
  supported_designs_J <- c(
    "mean_continuous",
    "mean_binary",
    "ttest_continuous",
    "paired_continuous"
  )

  plot_data_J <- results_n_inversion %>%
    filter(!is.na(analytical_achieved), design %in% supported_designs_J) %>%
    mutate(
      design_label = recode(
        design,
        mean_continuous = "Mean (Continuous)",
        mean_binary = "Mean (Binary)",
        ttest_continuous = "Two-Sample (Continuous)",
        paired_continuous = "Paired (Continuous)"
      )
    )

  plot_data_J_long <- plot_data_J %>%
    select(design_label, target_power, analytical_achieved, empirical_achieved) %>%
    pivot_longer(
      cols = c(analytical_achieved, empirical_achieved),
      names_to = "method",
      values_to = "achieved_power"
    ) %>%
    filter(!is.na(achieved_power)) %>%
    mutate(
      method = recode(
        method,
        analytical_achieved = "Analytical",
        empirical_achieved = "Empirical"
      ),
      target_power_label = sprintf("%.2f", target_power),
      deviation = achieved_power - target_power
    )

  fig_J <- ggplot(plot_data_J_long, aes(x = target_power_label, y = deviation, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = pppower_colors$reference, linewidth = 0.9) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.45),
      size = 3.2,
      alpha = 0.9
    ) +
    facet_wrap(~design_label, ncol = 2) +
    scale_color_manual(values = c(
      "Analytical" = pppower_colors$analytical,
      "Empirical" = pppower_colors$empirical
    )) +
    scale_y_continuous(breaks = seq(-0.10, 0.10, 0.05)) +
    coord_cartesian(ylim = c(-0.10, 0.10)) +
    labs(
      title = "Sample Size Inversion",
      x = "Target Power",
      y = "Achieved Power - Target Power",
      color = ""
    ) +
    theme_validation_focus()

  save_fig(fig_J, "fig_validation_J_n_inversion", width = 9, height = 6)
}

# =============================================================================
# Panel K: Rule-of-Thumb (Setting 11)
# =============================================================================
if (exists("results_rule_of_thumb") && is.data.frame(results_rule_of_thumb)) {
  plot_data_K <- results_rule_of_thumb %>%
    mutate(
      N_label = factor(
        paste0("N = ", formatC(N, big.mark = ",")),
        levels = paste0("N = ", formatC(sort(unique(N)), big.mark = ","))
      ),
      rule_gap = ratio - theoretical_ratio
    )

  n_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

  fig_K <- ggplot(plot_data_K, aes(x = rho_sq, y = rule_gap, color = N_label)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = pppower_colors$reference, linewidth = 0.9) +
    geom_line(linewidth = 1.3) +
    geom_point(size = 2.6) +
    scale_color_manual(values = n_colors) +
    scale_y_continuous(breaks = seq(-0.12, 0.12, 0.04)) +
    scale_x_continuous(breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Rule-of-Thumb Error",
      x = expression(rho^2),
      y = expression(frac(n[PPI], n[classical]) - (1 - rho^2)),
      color = "Unlabeled N"
    ) +
    theme_validation_focus()

  save_fig(fig_K, "fig_validation_K_rule_of_thumb", width = 8, height = 5.5)
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
        title = "Type I Error Calibration",
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
  plot_data_M <- results_lambda_convergence %>%
    mutate(rho_label = factor(paste0("rho = ", rho),
                              levels = paste0("rho = ", sort(unique(rho)))))

  fig_M <- ggplot(plot_data_M, aes(x = n, y = rmse, color = rho_label, group = rho_label)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = pppower_colors$reference, linewidth = 0.8) +
    geom_line(linewidth = 1.4) +
    geom_point(size = 3.6) +
    scale_color_manual(
      values = c(pppower_colors$ppi_pp, pppower_colors$classical, pppower_colors$vanilla)
    ) +
    scale_x_continuous(breaks = sort(unique(plot_data_M$n))) +
    labs(
      title = "Plugin Lambda Convergence",
      x = "Labeled Sample Size (n)",
      y = "RMSE Relative to Oracle Lambda",
      color = ""
    ) +
    theme_validation_focus()

  save_fig(fig_M, "fig_validation_M_lambda_convergence", width = 8, height = 5)
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
      title = "Plugin vs Oracle Lambda Power",
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
    geom_line(aes(y = power_emp_ppi, color = "Empirical"), linewidth = 1.0, alpha = 0.85) +
    facet_grid(n_label ~ rho_label,
               labeller = labeller(rho_label = label_parsed, n_label = label_parsed)) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray70") +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Power vs Effect Size",
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
      title = "Small n Regime",
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
      title = "N/n Ratio Sensitivity",
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
    arrange(rho, alloc_ratio) %>%
    mutate(
      rho_label = paste0("rho == ", rho),
      alloc_label = factor(
        paste0(n_A, ":", n_B),
        levels = unique(paste0(n_A, ":", n_B))
      )
    )

  fig_R <- ggplot(plot_data_R, aes(x = alloc_label, group = 1)) +
    geom_line(aes(y = power_emp_classical, color = "Classical"),
              linetype = "dashed", linewidth = 1.4) +
    geom_line(aes(y = power_theo_ppi, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_ppi, color = "Empirical"), size = 3.5) +
    facet_wrap(~rho_label, labeller = labeller(rho_label = label_parsed)) +
    scale_color_manual(values = colors_val) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "Unequal Group Sizes",
      x = expression(n[A]:n[B] ~ "allocation"),
      y = "Power",
      color = ""
    ) +
    theme_validation_focus()

  save_fig(fig_R, "fig_validation_R_unequal_groups", width = 10, height = 4)
}

# =============================================================================
# Panel Setting 19: Misspecified Prediction Quality
# =============================================================================
if (exists("results_setting19") && is.data.frame(results_setting19)) {
  plot_data_19 <- results_setting19 %>%
    mutate(plan_label = paste0("rho[plan] == ", sprintf("%.2f", rho_planning)))

  fig_19 <- ggplot(plot_data_19, aes(x = rho_true)) +
    geom_line(aes(y = power_theo_true, color = "Analytical"), linewidth = 1.8) +
    geom_point(aes(y = power_emp_true, color = "Empirical"), size = 3.4) +
    geom_hline(yintercept = unique(plot_data_19$power_target),
               linetype = "dashed", color = pppower_colors$reference, linewidth = 0.9) +
    facet_wrap(~plan_label, labeller = labeller(plan_label = label_parsed), scales = "free_x") +
    scale_color_manual(values = c(
      "Analytical" = pppower_colors$analytical,
      "Empirical" = pppower_colors$empirical
    )) +
    scale_x_continuous(
      breaks = sort(unique(plot_data_19$rho_true)),
      labels = function(x) sprintf("%.2f", x)
    ) +
    scale_y_continuous(limits = c(0.45, 1), breaks = seq(0.5, 1, 0.1)) +
    labs(
      title = "Misspecified Prediction Quality",
      x = expression(rho[true]),
      y = "Achieved Power at Planned n",
      color = ""
    ) +
    theme_validation_focus() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
    )

  save_fig(fig_19, "fig_validation_setting19_misspecified_rho", width = 10.5, height = 5.2)
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
      title = "OLS Regression Contrast",
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
      title = "GLM Logistic Regression Contrast",
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
    filter(p_exp == 0.40) %>%
    mutate(acc_label = sprintf("Accuracy = %.0f%%", accuracy * 100),
           N_label  = paste0("N == ", N))

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
      title = "2x2 Table: Odds Ratio",
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
    filter(p_exp == 0.40) %>%
    mutate(acc_label = sprintf("Accuracy = %.0f%%", accuracy * 100),
           N_label  = paste0("N == ", N))

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
      title = "2x2 Table: Relative Risk",
      x = "Total Labeled Sample Size (n)", y = "Power", color = ""
    ) +
    theme_validation()

  save_fig(fig_V, "fig_validation_V_2x2_rr")
}

cat(sprintf("Figures generated. (%.1f seconds)\n", difftime(Sys.time(), fig_start, units = "secs")))
cat(sprintf("Figures saved to: %s\n\n", normalizePath(output_dir, mustWork = FALSE)))
