# =============================================================================
# GENERATE VALIDATION FIGURES
# =============================================================================
cat("Generating validation figures...\n")
fig_start <- Sys.time()

library(ggplot2)
library(dplyr)
library(tidyr)

# Use centralized palette from the package
if (!exists("pppower_method_colors_full")) {
  tryCatch(source("R/palette.R", local = FALSE), error = function(e) {
    tryCatch(source(file.path(dirname(getwd()), "R", "palette.R"), local = FALSE),
             error = function(e2) NULL)
  })
}

# Aliases for backward compatibility
theme_validation <- theme_pppower
colors_method <- pppower_method_colors_full

# -----------------------------------------------------------------------------
# Panel A: Mean (Continuous)
# -----------------------------------------------------------------------------
plot_data_A <- results_cont_mean %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_A <- ggplot(plot_data_A, aes(x = n)) +
  # Theoretical lines
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  # Empirical points
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A. Mean Estimation (Continuous)",
    subtitle = sprintf("Delta = %.1f, N in {%s}, R = %d replications",
                       delta, paste(N_values, collapse = ", "), R),
    x = "Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel A2: Mean (Continuous) Power vs Effect Size
# -----------------------------------------------------------------------------
delta_grid <- seq(0, 0.6, by = 0.05)
plot_data_A_delta <- expand.grid(
  delta = delta_grid,
  n = settings_cont_mean$n_values,
  N = settings_cont_mean$N_values,
  rho = settings_cont_mean$rho_values
)
plot_data_A_delta <- plot_data_A_delta %>%
  mutate(
    cov_Yf = rho * sqrt(settings_cont_mean$sigma_Y2 * settings_cont_mean$sigma_f2)
  )

plot_data_A_delta$power_ppi <- mapply(
  theo_power_ppi_mean,
  delta = plot_data_A_delta$delta,
  n = plot_data_A_delta$n,
  N = plot_data_A_delta$N,
  sigma_Y2 = settings_cont_mean$sigma_Y2,
  sigma_f2 = settings_cont_mean$sigma_f2,
  cov_y_f = plot_data_A_delta$cov_Yf,
  MoreArgs = list(alpha = settings_cont_mean$alpha)
)

plot_data_A_delta$power_classical <- mapply(
  theo_power_classical_onesample,
  delta = plot_data_A_delta$delta,
  n = plot_data_A_delta$n,
  sigma_Y2 = settings_cont_mean$sigma_Y2,
  MoreArgs = list(alpha = settings_cont_mean$alpha)
)

# Empirical power for delta grid (Gaussian outcomes)
total_empirical <- nrow(plot_data_A_delta)
if (total_empirical > 0) {
  cat("Computing empirical power for delta grid (Setting A2)...\n")
}
plot_data_A_delta$power_emp_ppi <- NA_real_
plot_data_A_delta$power_emp_classical <- NA_real_
for (i in seq_len(nrow(plot_data_A_delta))) {
  row <- plot_data_A_delta[i, ]
  progress_bar(i, total_empirical, "Setting A2",
               sprintf("rho=%.1f, N=%d, n=%d, delta=%.2f: ",
                       row$rho, row$N, row$n, row$delta))

  sigma_f2 <- settings_cont_mean$sigma_f2
  cov_Yf <- row$cov_Yf
  Sigma <- matrix(c(settings_cont_mean$sigma_Y2, cov_Yf, cov_Yf, sigma_f2), 2, 2)
  L <- chol(Sigma)
  lambda_oracle <- cov_Yf / ((1 + row$n / row$N) * sigma_f2)

  rejections_ppi <- 0
  rejections_classical <- 0
  for (r in 1:settings_cont_mean$R) {
    Z_L <- matrix(rnorm(row$n * 2), row$n, 2) %*% L
    Y_L <- Z_L[, 1] + row$delta
    f_L <- Z_L[, 2] + row$delta

    Z_U <- matrix(rnorm(row$N * 2), row$N, 2) %*% L
    f_U <- Z_U[, 2] + row$delta

    test_ppi <- run_ppi_mean_test(Y_L, f_L, f_U, theta0 = 0, alpha = settings_cont_mean$alpha,
                                  lambda_oracle = lambda_oracle)
    rejections_ppi <- rejections_ppi + test_ppi$reject

    test_classical <- classical_mean_test(Y_L, theta0 = 0, alpha = settings_cont_mean$alpha)
    rejections_classical <- rejections_classical + test_classical$reject
  }
  plot_data_A_delta$power_emp_ppi[i] <- rejections_ppi / settings_cont_mean$R
  plot_data_A_delta$power_emp_classical[i] <- rejections_classical / settings_cont_mean$R
}

plot_data_A_delta_long <- plot_data_A_delta %>%
  select(delta, n, N, rho, power_ppi, power_classical,
         power_emp_ppi, power_emp_classical) %>%
  mutate(
    n_label = paste0("n = ", n),
    rho_label = paste0("rho == ", rho),
    N_label = paste0("N == ", N)
  ) %>%
  pivot_longer(
    cols = c(power_ppi, power_classical, power_emp_ppi, power_emp_classical),
    names_to = "method",
    values_to = "power"
  ) %>%
  mutate(
    method = recode(
      method,
      power_ppi = "PPI++ (Theoretical)",
      power_classical = "Classical (Theoretical)",
      power_emp_ppi = "PPI++ (Empirical)",
      power_emp_classical = "Classical (Empirical)"
    )
  )

fig_A_delta <- ggplot(plot_data_A_delta_long, aes(x = delta, y = power)) +
  geom_line(
    data = subset(plot_data_A_delta_long, grepl("Theoretical", method)),
    aes(color = method, linetype = n_label),
    linewidth = 1
  ) +
  geom_point(
    data = subset(plot_data_A_delta_long, grepl("Empirical", method)),
    aes(color = method, shape = n_label),
    size = 2.2
  ) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "A2. Mean Estimation (Continuous): Power vs Effect Size",
    subtitle = sprintf("Gaussian outcomes; N in {%s}, n in {%s}",
                       paste(settings_cont_mean$N_values, collapse = ", "),
                       paste(settings_cont_mean$n_values, collapse = ", ")),
    x = "Effect Size (Delta)",
    y = "Power",
    color = "Method",
    linetype = "Labeled n",
    shape = "Labeled n"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel B: Mean (Binary)
# -----------------------------------------------------------------------------
plot_data_B <- results_bin_mean %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_B <- ggplot(plot_data_B, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "B. Mean Estimation (Binary/Prevalence)",
    subtitle = sprintf("Delta = %.2f, p0 = %.1f, N in {%s}",
                       delta_bin, p0, paste(N_bin_values, collapse = ", ")),
    x = "Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel C: t-Test (Continuous)
# -----------------------------------------------------------------------------
plot_data_C <- results_cont_ttest %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_C <- ggplot(plot_data_C, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "C. Two-Sample t-Test (Continuous)",
    subtitle = sprintf("Delta = %.1f, N in {%s} per group",
                       delta_ttest, paste(N_ttest_values, collapse = ", ")),
    x = "Per-Group Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel D: Proportion Test (Binary)
# -----------------------------------------------------------------------------
plot_data_D <- results_bin_ttest %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_D <- ggplot(plot_data_D, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "D. Two-Sample Proportion Test (Binary)",
    subtitle = sprintf("Delta = %.2f, p0 = %.1f, N in {%s} per group",
                       delta_prop, p0_prop, paste(N_prop_values, collapse = ", ")),
    x = "Per-Group Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel E: Paired t-Test (Continuous)
# -----------------------------------------------------------------------------
plot_data_E <- results_paired %>%
  mutate(rho_label = paste0("rho[D] == ", rho_D),
         N_label = paste0("N == ", N))

fig_E <- ggplot(plot_data_E, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "E. Paired t-Test (Continuous)",
    subtitle = sprintf("Delta = %.1f, N in {%s}",
                       delta_paired, paste(N_paired_values, collapse = ", ")),
    x = "Number of Pairs (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel F: Paired Proportion Test (Binary)
# -----------------------------------------------------------------------------
plot_data_F <- results_paired_bin %>%
  mutate(classifier_label = sprintf("Sens/Spec = %.0f%%", sens * 100),
         N_label = paste0("N == ", N))

fig_F <- ggplot(plot_data_F, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ classifier_label, labeller = labeller(N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "F. Paired Proportion Test (Binary)",
    subtitle = sprintf("Delta = %.2f, p0 = %.1f, N in {%s}",
                       delta_paired_bin, p0_paired, paste(N_paired_bin_values, collapse = ", ")),
    x = "Number of Pairs (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel G: Log-Normal (Skewed)
# -----------------------------------------------------------------------------
plot_data_G <- results_lognorm %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_G <- ggplot(plot_data_G, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "G. Mean Estimation (Log-Normal, Skewed)",
    subtitle = sprintf("Delta = %.3f, N in {%s}",
                       delta_lognorm, paste(N_lognorm_values, collapse = ", ")),
    x = "Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Panel H: t-distribution (Heavy-Tailed)
# -----------------------------------------------------------------------------
plot_data_H <- results_tdist %>%
  mutate(rho_label = paste0("rho == ", rho),
         N_label = paste0("N == ", N))

fig_H <- ggplot(plot_data_H, aes(x = n)) +
  geom_line(aes(y = power_theo_ppi, color = "PPI++ (Theoretical)"), linewidth = 1) +
  geom_line(aes(y = power_theo_classical, color = "Classical (Theoretical)"),
            linewidth = 1, linetype = "dashed") +
  geom_point(aes(y = power_emp_ppi, color = "PPI++ (Empirical)"), size = 3) +
  geom_point(aes(y = power_emp_classical, color = "Classical (Empirical)"),
             size = 3, shape = 17) +
  facet_grid(N_label ~ rho_label,
             labeller = labeller(rho_label = label_parsed, N_label = label_parsed)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray50") +
  scale_color_manual(values = colors_method) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    title = "H. Mean Estimation (t-dist df=5, Heavy-Tailed)",
    subtitle = sprintf("Delta = %.1f, N in {%s}",
                       delta_t, paste(N_t_values, collapse = ", ")),
    x = "Labeled Sample Size (n)",
    y = "Power",
    color = "Method"
  ) +
  theme_validation()

# -----------------------------------------------------------------------------
# Save individual figures as PDFs to write-up folder
# -----------------------------------------------------------------------------
output_dir <- "../write-up/PPI_Power_Analysis/figures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(file.path(output_dir, "fig_validation_A_mean_continuous.pdf"), fig_A,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_A2_mean_continuous_delta.pdf"), fig_A_delta,
       width = 10, height = 4.8, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_B_mean_binary.pdf"), fig_B,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_C_ttest_continuous.pdf"), fig_C,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_D_ttest_binary.pdf"), fig_D,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_E_paired_continuous.pdf"), fig_E,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_F_paired_binary.pdf"), fig_F,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_G_lognormal.pdf"), fig_G,
       width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "fig_validation_H_tdist.pdf"), fig_H,
       width = 10, height = 4, dpi = 300)

# -----------------------------------------------------------------------------
# Panel I: EIF Binary Validation (Setting 9)
# -----------------------------------------------------------------------------
if (exists("results_eif_binary") && is.data.frame(results_eif_binary)) {
  plot_data_I <- results_eif_binary
  plot_data_I$clf_label <- sprintf("Sens/Spec = %.0f%%", plot_data_I$sens * 100)
  plot_data_I$p_label <- paste0("p = ", plot_data_I$p)

  fig_I <- ggplot(plot_data_I, aes(x = m_cal)) +
    geom_line(aes(y = power_theo, color = "PPI++ (Theoretical)"), linewidth = 1) +
    geom_point(aes(y = power_emp, color = "PPI++ (Empirical)"), size = 2.5) +
    geom_line(aes(y = power_classical, color = "Classical (Theoretical)"),
              linewidth = 1, linetype = "dashed") +
    facet_grid(p_label ~ clf_label) +
    geom_hline(yintercept = 0.8, linetype = "dotted", color = pppower_colors$power_target) +
    scale_color_manual(values = colors_method) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(
      title = "I. EIF Binary Surrogate Validation",
      subtitle = sprintf("Delta = %.2f, R = %d replications", 0.05, 1000),
      x = "Calibration Sample Size (m_cal)",
      y = "Power",
      color = "Method"
    ) +
    theme_validation()

  ggsave(file.path(output_dir, "fig_validation_I_eif_binary.pdf"), fig_I,
         width = 10, height = 5, dpi = 300)
}

# -----------------------------------------------------------------------------
# Panel J: Sample Size Inversion (Setting 10)
# -----------------------------------------------------------------------------
if (exists("results_n_inversion") && is.data.frame(results_n_inversion)) {
  plot_data_J <- results_n_inversion[!is.na(results_n_inversion$analytical_achieved), ]

  fig_J <- ggplot(plot_data_J, aes(x = target_power, y = analytical_achieved)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = pppower_colors$reference) +
    geom_point(aes(color = design, shape = design), size = 3) +
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
      x = "Target Power",
      y = "Achieved Power at n*",
      color = "Design",
      shape = "Design"
    ) +
    theme_validation()

  ggsave(file.path(output_dir, "fig_validation_J_n_inversion.pdf"), fig_J,
         width = 7, height = 5, dpi = 300)
}

# -----------------------------------------------------------------------------
# Panel K: Rule-of-Thumb (Setting 11)
# -----------------------------------------------------------------------------
if (exists("results_rule_of_thumb") && is.data.frame(results_rule_of_thumb)) {
  plot_data_K <- results_rule_of_thumb
  plot_data_K$N_label <- paste0("N = ", formatC(plot_data_K$N, big.mark = ","))

  fig_K <- ggplot(plot_data_K, aes(x = rho_sq)) +
    geom_line(aes(y = ratio, color = N_label), linewidth = 1) +
    geom_line(aes(y = theoretical_ratio), linetype = "dashed",
              color = pppower_colors$reference, linewidth = 0.8) +
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

  ggsave(file.path(output_dir, "fig_validation_K_rule_of_thumb.pdf"), fig_K,
         width = 7, height = 5, dpi = 300)
}

# -----------------------------------------------------------------------------
# Panel L: Type I Error (Setting 12)
# -----------------------------------------------------------------------------
if (exists("results_type1_error") && is.list(results_type1_error)) {
  ci_lo <- 0.05 - 1.96 * sqrt(0.05 * 0.95 / 2000)
  ci_hi <- 0.05 + 1.96 * sqrt(0.05 * 0.95 / 2000)

  if (!is.null(results_type1_error$continuous)) {
    plot_data_L <- results_type1_error$continuous
    plot_data_L$rho_label <- paste0("rho = ", plot_data_L$rho)
    plot_data_L$N_label <- paste0("N = ", plot_data_L$N)

    fig_L <- ggplot(plot_data_L, aes(x = n)) +
      geom_point(aes(y = rejection_ppi, color = "PPI++"), size = 2.5) +
      geom_point(aes(y = rejection_classical, color = "Classical"),
                 size = 2.5, shape = 17) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
      geom_hline(yintercept = ci_lo, linetype = "dotted", color = pppower_colors$reference) +
      geom_hline(yintercept = ci_hi, linetype = "dotted", color = pppower_colors$reference) +
      facet_grid(N_label ~ rho_label) +
      scale_color_manual(values = pppower_method_colors) +
      scale_y_continuous(limits = c(0, 0.12), breaks = seq(0, 0.12, 0.02)) +
      labs(
        title = "L. Type I Error Calibration (Gaussian Mean, H0)",
        subtitle = "Dashed = nominal 0.05; dotted = 95% MC CI",
        x = "Labeled Sample Size (n)",
        y = "Rejection Rate",
        color = "Method"
      ) +
      theme_validation()

    ggsave(file.path(output_dir, "fig_validation_L_type1_continuous.pdf"), fig_L,
           width = 10, height = 5, dpi = 300)
  }
}

# -----------------------------------------------------------------------------
# Panel M: Lambda Convergence (Setting 13)
# -----------------------------------------------------------------------------
if (exists("results_lambda_convergence") && is.data.frame(results_lambda_convergence)) {
  plot_data_M <- results_lambda_convergence
  plot_data_M$rho_label <- paste0("rho == ", plot_data_M$rho)

  fig_M <- ggplot(plot_data_M, aes(x = n)) +
    geom_ribbon(aes(ymin = lambda_hat_mean - lambda_hat_sd,
                    ymax = lambda_hat_mean + lambda_hat_sd),
                alpha = 0.2, fill = pppower_colors$ppi_pp) +
    geom_line(aes(y = lambda_hat_mean, color = "Plugin mean"), linewidth = 1) +
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
      y = expression(lambda),
      color = ""
    ) +
    theme_validation()

  ggsave(file.path(output_dir, "fig_validation_M_lambda_convergence.pdf"), fig_M,
         width = 10, height = 4, dpi = 300)
}

cat(sprintf("Figures generated. (%.1f seconds)\n", difftime(Sys.time(), fig_start, units = "secs")))
cat("Figures saved to write-up/PPI_Power_Analysis/figures/\n\n")
