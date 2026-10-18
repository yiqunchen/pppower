#!/usr/bin/env Rscript
# =============================================================================
# Figures and tables for Setting 23 (estimated tuning parameter)
# =============================================================================
# Usage: Rscript simulation_studies/setting23_figures.R
# Reads rds_cache/power_validation_setting23_plugin_lambda.rds and writes
#   fig_validation_W_plugin_type1.{pdf,png}
#   fig_validation_X_plugin_power.{pdf,png}
#   table_setting23_variance.tex
#   table_setting23_planning.tex

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
script_dir <- if (dir.exists("simulation_studies")) "simulation_studies" else "."
if (!exists("theme_pppower")) {
  if (requireNamespace("pppower", quietly = TRUE)) {
    theme_pppower <- pppower::theme_pppower
  } else source(file.path(script_dir, "..", "R", "palette.R"))
}

x <- readRDS(file.path(script_dir, "rds_cache", "power_validation_setting23_plugin_lambda.rds"))
r <- x$results_setting23
p <- x$plan_setting23

method_labels <- c(
  oracle     = "Oracle lambda",
  plugin     = "Plugin lambda",
  plugin_adj = "Plugin, corrected variance",
  crossfit   = "Cross-fitted lambda",
  pilot      = "Lambda from pilot",
  classical  = "Classical"
)
# Okabe-Ito colors, distinct for every method; classical in gray
method_colors <- c(
  "Oracle lambda"              = "#000000",
  "Plugin lambda"              = "#E69F00",
  "Plugin, corrected variance" = "#D55E00",
  "Cross-fitted lambda"        = "#CC79A7",
  "Lambda from pilot"          = "#0072B2",
  "Classical"                  = "gray55"
)
method_shapes <- c(
  "Oracle lambda" = 16, "Plugin lambda" = 17, "Plugin, corrected variance" = 15,
  "Cross-fitted lambda" = 18, "Lambda from pilot" = 1, "Classical" = 4
)
family_levels <- c("Gaussian", "t5", "Log-normal (sd 0.5)", "Log-normal (sd 1)",
                   "Binary (p = 0.35)", "Binary (p = 0.05)")

# One representative cell per family: rho = 0.5 (sens = spec = 0.85 for
# binary), and the larger unlabeled pool.
pick <- function(d) {
  d %>%
    group_by(family) %>%
    filter((is.na(rho_target) & sens == 0.85) | (!is.na(rho_target) & rho_target == 0.5),
           N == max(N)) %>%
    ungroup() %>%
    mutate(family = factor(family, levels = family_levels),
           method_label = factor(method_labels[method], levels = method_labels),
           design = factor(design, levels = c("one-sample", "two-sample"),
                           labels = c("One-sample mean", "Two-sample difference")),
           lo = pmax(0, reject - 1.96 * mc_se), hi = pmin(1, reject + 1.96 * mc_se))
}

theme_fig <- theme_pppower(base_size = 13) +
  theme(legend.position = "bottom", strip.text = element_text(face = "bold"),
        axis.text = element_text(face = "bold"), panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------
# Type I error
# ---------------------------------------------------------------------------
d0 <- pick(r %>% filter(hyp == "H0"))
fig_W <- ggplot(d0, aes(n, reject, color = method_label, shape = method_label)) +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "gray40") +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5, alpha = 0.7) +
  geom_point(size = 2.2) +
  facet_grid(design ~ family, scales = "free_x") +
  scale_x_continuous(trans = "log10", breaks = c(15, 30, 50, 100, 200, 400)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = "Labeled sample size n (per group)", y = "Type I error", color = NULL, shape = NULL) +
  theme_fig + guides(color = guide_legend(nrow = 1))

# ---------------------------------------------------------------------------
# Power: empirical (points) against analytical oracle (solid) and plug-in
# corrected (dashed) curves
# ---------------------------------------------------------------------------
d1 <- pick(r %>% filter(hyp == "H1"))
curves <- d1 %>% distinct(family, design, n, power_theo_oracle, power_theo_plugin, power_theo_classical)
fig_X <- ggplot(d1 %>% filter(method %in% c("plugin", "plugin_adj", "pilot", "classical")),
                aes(n, reject, color = method_label, shape = method_label)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray40") +
  geom_line(data = curves, aes(n, power_theo_oracle), inherit.aes = FALSE,
            color = "black", linewidth = 0.9) +
  geom_line(data = curves, aes(n, power_theo_plugin), inherit.aes = FALSE,
            color = "black", linewidth = 0.9, linetype = "dashed") +
  geom_line(data = curves, aes(n, power_theo_classical), inherit.aes = FALSE,
            color = "gray55", linewidth = 0.9, linetype = "dashed") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, linewidth = 0.5, alpha = 0.7,
                position = position_dodge(width = 0.06)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.06)) +
  facet_grid(design ~ family, scales = "free_x") +
  scale_x_continuous(trans = "log10", breaks = c(15, 30, 50, 100, 200, 400)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = "Labeled sample size n (per group)", y = "Power", color = NULL, shape = NULL) +
  theme_fig + guides(color = guide_legend(nrow = 1))

ggsave(file.path(script_dir, "fig_validation_W_plugin_type1.pdf"), fig_W, width = 14, height = 5.6)
ggsave(file.path(script_dir, "fig_validation_W_plugin_type1.png"), fig_W, width = 14, height = 5.6, dpi = 300, bg = "white")
ggsave(file.path(script_dir, "fig_validation_X_plugin_power.pdf"), fig_X, width = 14, height = 5.6)
ggsave(file.path(script_dir, "fig_validation_X_plugin_power.png"), fig_X, width = 14, height = 5.6, dpi = 300, bg = "white")

# ---------------------------------------------------------------------------
# Variance table (one-sample design, representative cell)
# ---------------------------------------------------------------------------
v <- r %>%
  filter(hyp == "H1", design == "one-sample") %>%
  group_by(family) %>%
  filter((is.na(rho_target) & sens == 0.85) | (!is.na(rho_target) & rho_target == 0.5),
         N == max(N), n %in% c(min(n), sort(unique(n))[3], max(n))) %>%
  ungroup() %>%
  filter(method %in% c("plugin", "plugin_adj", "pilot")) %>%
  select(family, n, N, method, mc_var, V_opt, V_plugin, mean_vhat, coverage) %>%
  pivot_wider(names_from = method, values_from = c(mc_var, mean_vhat, coverage)) %>%
  mutate(family = factor(family, levels = family_levels)) %>%
  arrange(family, n) %>%
  transmute(family, n,
            ratio_opt = mc_var_plugin / V_opt, ratio_pl = mc_var_plugin / V_plugin,
            vhat_plugin = mean_vhat_plugin / mc_var_plugin,
            vhat_adj = mean_vhat_plugin_adj / mc_var_plugin_adj,
            vhat_pilot = mean_vhat_pilot / mc_var_pilot,
            cov_plugin = coverage_plugin, cov_adj = coverage_plugin_adj, cov_pilot = coverage_pilot)

fmt <- function(z) sprintf("%.2f", z)
tex <- c("\\begin{tabular}{llrrrrrrrr}", "\\toprule",
         " & & \\multicolumn{2}{c}{Variance of plugin estimator} & \\multicolumn{3}{c}{Mean estimated / true variance} & \\multicolumn{3}{c}{Coverage of 95\\% interval} \\\\",
         "\\cmidrule(lr){3-4}\\cmidrule(lr){5-7}\\cmidrule(lr){8-10}",
         "Outcome & $n$ & / asymptotic & / corrected & plugin & corrected & pilot & plugin & corrected & pilot \\\\", "\\midrule")
for (i in seq_len(nrow(v))) {
  w <- v[i, ]
  tex <- c(tex, paste(c(as.character(w$family), w$n, fmt(w$ratio_opt), fmt(w$ratio_pl), fmt(w$vhat_plugin),
                        fmt(w$vhat_adj), fmt(w$vhat_pilot), fmt(w$cov_plugin), fmt(w$cov_adj), fmt(w$cov_pilot)),
                      collapse = " & "), "\\\\")
  tex[length(tex) - 1] <- paste(tex[length(tex) - 1], tex[length(tex)]); tex <- tex[-length(tex)]
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_setting23_variance.tex"))

# ---------------------------------------------------------------------------
# Planning table: planned n and achieved power at the oracle-planned n
# ---------------------------------------------------------------------------
pl <- p %>%
  filter(plan == "oracle", method %in% c("plugin", "plugin_adj", "pilot")) %>%
  group_by(family, design) %>%
  summarise(n_range = sprintf("%d to %d", min(n_oracle_plan), max(n_oracle_plan)),
            extra = sprintf("%d to %d", min(n_plugin_plan - n_oracle_plan), max(n_plugin_plan - n_oracle_plan)),
            plugin = sprintf("%.2f (%.2f)", mean(achieved_power[method == "plugin"]), min(achieved_power[method == "plugin"])),
            adj = sprintf("%.2f (%.2f)", mean(achieved_power[method == "plugin_adj"]), min(achieved_power[method == "plugin_adj"])),
            pilot = sprintf("%.2f (%.2f)", mean(achieved_power[method == "pilot"]), min(achieved_power[method == "pilot"])),
            .groups = "drop") %>%
  mutate(family = factor(family, levels = family_levels)) %>% arrange(design, family)
tex <- c("\\begin{tabular}{llllll}", "\\toprule",
         "Outcome & Design & Planned $n^\\star$ & Extra labels & \\multicolumn{2}{c}{Achieved power, mean (min)} \\\\",
         "\\cmidrule(lr){5-6}", " & & & (plugin correction) & plugin & pilot $\\lambda$ \\\\", "\\midrule")
for (i in seq_len(nrow(pl))) {
  w <- pl[i, ]
  tex <- c(tex, paste0(paste(as.character(w$family), w$design, w$n_range, w$extra, w$plugin, w$pilot, sep = " & "), " \\\\"))
}
tex <- c(tex, "\\bottomrule", "\\end{tabular}")
writeLines(tex, file.path(script_dir, "table_setting23_planning.tex"))

# ---------------------------------------------------------------------------
# Figures for the response letter (Reviewer 1, comments 1(a) and 3), without the
# strongly skewed log-normal outcome (sd 1), which is kept in the Supplement.
# ---------------------------------------------------------------------------
resp_families <- setdiff(family_levels, "Log-normal (sd 1)")
fig_R1_W <- fig_W %+% (d0 %>% filter(family != "Log-normal (sd 1)") %>% droplevels())

# Achieved power at the planned n (closed form, oracle lambda in the plan)
pa <- p %>%
  filter(plan == "oracle", method %in% c("plugin", "plugin_adj", "pilot"), family %in% resp_families) %>%
  mutate(family = factor(family, levels = resp_families, labels = sub("^t5$", "t (5 df)", resp_families)),
         method_label = factor(method_labels[method], levels = method_labels),
         design = factor(design, levels = c("one-sample", "two-sample"),
                         labels = c("One-sample mean", "Two-sample difference")))
fig_R1_power <- ggplot(pa, aes(family, achieved_power, color = method_label, shape = method_label)) +
  geom_hline(yintercept = 0.8, linetype = "dotted", color = "gray40") +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6, seed = 1),
             size = 1.8, alpha = 0.75) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.4,
               position = position_dodge(width = 0.6), show.legend = FALSE) +
  facet_wrap(~ design, nrow = 1) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  scale_y_continuous(limits = c(0.6, 1), breaks = seq(0.6, 1, 0.1)) +
  labs(x = NULL, y = "Achieved power at the planned n", color = NULL, shape = NULL) +
  theme_fig + theme(axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(color = guide_legend(nrow = 1))

# Mean estimated variance over Monte Carlo variance, one-sample design
vr <- r %>%
  filter(hyp == "H1", design == "one-sample", family %in% resp_families,
         method %in% c("plugin", "plugin_adj", "pilot")) %>%
  group_by(family) %>%
  filter((is.na(rho_target) & sens == 0.85) | (!is.na(rho_target) & rho_target == 0.5), N == max(N)) %>%
  ungroup() %>%
  mutate(family = factor(family, levels = resp_families),
         method_label = factor(method_labels[method], levels = method_labels),
         ratio = mean_vhat / mc_var)
fig_R1_var <- ggplot(vr, aes(n, ratio, color = method_label, shape = method_label)) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "gray40") +
  geom_line(linewidth = 1) + geom_point(size = 2.2) +
  facet_wrap(~ family, nrow = 1, scales = "free_x") +
  scale_x_continuous(trans = "log10", breaks = c(15, 30, 50, 100, 200, 400)) +
  scale_color_manual(values = method_colors) +
  scale_shape_manual(values = method_shapes) +
  labs(x = "Labeled sample size n", y = "Estimated / true variance", color = NULL, shape = NULL) +
  theme_fig + guides(color = guide_legend(nrow = 1))

ggsave(file.path(script_dir, "fig_response_R1_type1.pdf"), fig_R1_W, width = 12, height = 5.6)
ggsave(file.path(script_dir, "fig_response_R1_power.pdf"), fig_R1_power, width = 12, height = 4.6)
ggsave(file.path(script_dir, "fig_response_R1_variance.pdf"), fig_R1_var, width = 12, height = 3.8)
cat("Setting 23 figures and tables written\n")
