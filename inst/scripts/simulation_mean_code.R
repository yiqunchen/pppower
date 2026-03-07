library(furrr)
library(purrr)
library(tidyr)
library(utils)
library(future)
library(dplyr)
library(ggplot2)
library(pppower)

# True parameter
X_sampler_L <- function(n) data.frame(x1 = rnorm(n), x2 = rnorm(n))
f_true <- function(X) with(X, x1 + x1 * x2 + x2^3)


########################################################### Estimating lambda on an external train set, then do power grid test

n_grid <- c(50, 100, 150, 200, 300)
N <- 2000
R <- 500  
deltas <- seq(0, 0.6, by = 0.05)
n_external <- c(50, 100, 200, 300, 500, 1000)
models <- c("glm_correct", "glm_mis", "glm_wrong")
models <- c("glm_correct", "glm_mis", "glm_wrong", "rf")
ppi_types <- c("PPI", "PPI++")

grid <- expand.grid(
  delta = deltas,
  model_type = models,
  ppi_type = ppi_types,
  n = n_grid,
  lambda_mode = "plugin",
  lambda_external = TRUE,
  stringsAsFactors = FALSE
)

## Compare with the default, X-fit lambda estimation
# grid <- expand.grid(
#   delta = deltas,
#   model_type = models,
#   ppi_type = ppi_types,
#   n = n_grid,
#   lambda_mode = "plugin",
#   lambda_external = FALSE,
#   stringsAsFactors = FALSE
# )

# Parallelization
plan(multisession, workers = 16) 

power_grid <- grid %>%
  mutate(
    result = furrr::future_pmap(
      list(model_type, delta, ppi_type, n, lambda_mode, lambda_external),
      ~ ppi_empirical_power(
        R = R,                              # Monte Carlo replicates
        n = ..4,                            # labeled sample size
        N = N,                              # unlabeled sample size
        X_sampler_L = X_sampler_L,          # covariate sampler (DGP)
        f_generator = f_true,               # true f(X)
        model_type = ..1,                   # model
        delta = ..2,                        # mean shift
        alpha = 0.05,                       # test level
        ppi_type = ..3,                     # PPI or PPI++
        lambda_mode = ..5,                  # plugin / theory
        lambda_external = ..6,              # FALSE: normal; TRUE: external lambda
        seed = 42                           # reproducible seed
      ),
      .options = furrr::furrr_options(seed = TRUE),
      .progress = TRUE
    )
  )


power_df <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  select(model_type, ppi_type, n, delta, empirical_power, mc_se, avg_SE, lambda_external) %>%
  mutate(
    ppi_type = factor(ppi_type, levels = c("naive", "PPI", "PPI++")),
    n        = factor(n, levels = sort(unique(n))),  # columns by n
    model_type = factor(model_type, levels = c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  )

lambda_summary <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  select(model_type, ppi_type, n, delta, avg_lambda_hat, sd_lambda_hat, lambda_external) %>%
  group_by(model_type, ppi_type, n, lambda_external) %>%
  summarise(
    mean_lambda = mean(avg_lambda_hat, na.rm = TRUE),
    sd_lambda   = mean(sd_lambda_hat, na.rm = TRUE),
    .groups = "drop"
  )

# # Normal X-fit lambda by changing `lambda_external`
# default_lambda_summary <- lambda_summary
# default_lambda_summary %>% filter(ppi_type == "PPI++")

power_df_ext <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  dplyr::select(
    all_of(c(
      "model_type", "ppi_type", "n", "delta", "empirical_power",
      "mc_se", "avg_lambda_hat", "avg_SE", "lambda_external", "n_external"
    ))
  ) %>%
  mutate(
    n_external_label = ifelse(lambda_external,
                              paste0("ext=", n_external),
                              "plugin (no external λ)")
  )

p_lambda <- ggplot(
  plot_df,
  aes(
    x = delta,
    y = empirical_power,
    color = as.character(n),
    group = interaction(n, n_external)  # ensure lines separate by n and external size
  )
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  facet_grid(model_type ~ interaction(ppi_type, n_external), scales = "free_y") +
  labs(
    x = expression(Delta),
    y = "Empirical Power",
    color = "Labeled sample size (n)",
    title = "Empirical Power vs Effect Size",
    subtitle = "Within each facet: curves for different n; columns = (Estimator * External λ size)"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

# ggsave("external_lambda_verification_Xfit_n_corrected_th.png", p_lambda, width = 13, height = 7, dpi = 300)

####################################### normal X-fit, power curve comparison

# Rows: Powers of PPI, PPI++, and PPI++ - PPI; Columns: sample sizes

power_df_ext <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  dplyr::select(
    all_of(c(
      "model_type", "ppi_type", "n", "delta", "empirical_power",
      "mc_se", "avg_lambda_hat", "avg_SE", "lambda_external"
    ))
  )


plot_df <- power_df_ext %>%
  filter(ppi_type %in% c("PPI", "PPI++")) %>%
  select(model_type, n, ppi_type, delta, empirical_power, mc_se)

# Compute difference row (PPI++ − PPI)

plot_diff <- plot_df %>%
  pivot_wider(names_from = ppi_type, values_from = c(empirical_power, mc_se)) %>%
  mutate(
    empirical_power = `empirical_power_PPI++` - empirical_power_PPI,
    mc_se = sqrt(`mc_se_PPI++`^2 + `mc_se_PPI`^2),  # conservative SE combination
    ppi_type = "PPI++ - PPI"
  ) %>%
  select(model_type, n, delta, ppi_type, empirical_power, mc_se)

# Combine all three (PPI, PPI++, Difference)
plot_full <- bind_rows(plot_df, plot_diff) %>%
  mutate(
    ppi_type = factor(ppi_type, levels = c("PPI", "PPI++", "PPI++ - PPI")),
    model_type = factor(model_type,
                        levels = c("glm_correct", "glm_mis", "glm_wrong", "rf"),
                        labels = c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  )

# Figure
p_power_cf <- ggplot(
  plot_full,
  aes(x = delta, y = empirical_power,
      color = model_type, group = model_type)
) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = empirical_power - 1.96 * mc_se, ymax = empirical_power + 1.96 * mc_se),
    width = 0.015, alpha = 0.5
  ) +
  facet_grid(ppi_type ~ n, scales = "free_y") +
  labs(
    x = expression("Effect size " * delta),
    y = "Empirical Power (or difference in power)",
    color = "Model",
    title = "Empirical Power estimated via Monte Carlo",
    subtitle = "Rows: estimator type; last row shows PPI++ - PPI; Columns: labeled sample size n",
    caption = "labmda for PPI++ estimator were estimated using an independent external set of labeled data, size n/2."
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10)
  )

p_power_cf

#ggsave("Power_curve_external_lambda_sample_estimates.png", p_power_cf, width = 13, height = 7, dpi = 300)

####################################### Empirical power vs. Theoretical power

alpha <- 0.05
z_alpha <- qnorm(1 - alpha/2)
N_big <- 2000   # same unlabeled size as used in simulation

power_compare <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  select(model_type, ppi_type, n, delta,
         empirical_power, theoretical_power, mc_se, avg_SE, lambda_external) %>%
  mutate(
    ppi_type = factor(ppi_type, levels = c("naive", "PPI", "PPI++")),
    n        = factor(n, levels = sort(unique(n))),
    model_type = factor(
      model_type,
      levels = c("glm_correct", "glm_mis", "glm_wrong", "rf")
    )
  )

emp_theoretical_power <- ggplot(
  power_compare %>% filter(ppi_type != "naive"),
  aes(x = theoretical_power, y = empirical_power, color = model_type)
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 2, alpha = 0.6) +
  facet_grid(ppi_type ~ n, scales = "free") +
  labs(
    title = "Empirical vs. Analytical (Sample-based) Power",
    subtitle = "Rows: estimator type; Columns: labeled sample size n",
    x = "Closed-form (sample-based) power",
    y = "Empirical (simulated) power",
    color = "Model",
    caption = "λ estimated from an independent external labeled set of size n/2."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0.5, face = "italic", size = 10)
  )


ggsave("Emp_theoretical_power_compare_internal_lambda_sample_estimates_with_rf.png", emp_theoretical_power, width = 11, height = 7, dpi = 300)


####################################### for testing all scenarios, X-fit/Non-Xfit/Independent-Split for PPI/PPI++

# Grid Test
n_grid   <- c(50, 100, 150, 200, 300)
N <- 2000
R <- 500  
deltas <- seq(0, 0.6, by = 0.05)
models <- c("glm_correct", "glm_mis", "glm_wrong", "rf")
#models <- c("glm_correct", "glm_mis", "glm_wrong")
ppi_types <- c("PPI", "PPI++")

############################# for testing all scenarios

design_flags <- expand.grid(
  PPI_crossfit = c(TRUE, FALSE),
  PPIpp_crossfit = c(TRUE, FALSE),
  PPIpp_independent_split = c(TRUE, FALSE),
  stringsAsFactors = FALSE
) %>%
  ## For PPI++, we only keep valid combinations:
  ##  - if PPIpp_crossfit == TRUE → PPIpp_independent_split must be FALSE
  ##  - if PPIpp_independent_split == TRUE → PPIpp_crossfit must be FALSE
  ##  - for PPI, PPIpp flags are ignored (but we’ll keep them to match structure)
  filter(!(PPIpp_crossfit & PPIpp_independent_split))


grid <- expand.grid(
  delta = deltas,
  model_type = models,
  ppi_type = ppi_types,
  n = n_grid,
  lambda_mode = "plugin",
  lambda_external = TRUE,
  stringsAsFactors = FALSE
) %>%
  crossing(design_flags)

plan(multicore, workers = 16)

############################# grid test for all scenarios

power_grid <- grid %>%
  mutate(result = furrr::future_pmap(
    list(model_type, delta, ppi_type, n, lambda_mode, lambda_external,
         PPI_crossfit, PPIpp_crossfit, PPIpp_independent_split),
    ~ ppi_empirical_power(
      R = R,
      n = ..4,
      N = N,
      X_sampler_L = X_sampler_L,
      f_generator = f_true,
      model_type = ..1,
      delta = ..2,
      alpha = 0.05,
      ppi_type = ..3,
      lambda_mode = ..5,
      lambda_external = ..6,
      seed = 42,
      PPI_crossfit = ..7,
      PPIpp_crossfit = ..8,
      PPIpp_independent_split = ..9
    ),
    .options = furrr_options(seed = TRUE, globals = FALSE),
    .progress = TRUE
  ))

power_grid_all_comb <- power_grid


# Raw file
# saveRDS(power_grid_all_comb, file = "power_grid_all_comb.rds")

# Read file
# power_all <- readRDS("./power_grid_all_comb.rds")

power_summary <- power_grid_all_comb %>%
  mutate(
    emp_power = map_dbl(result, "empirical_power"),
    mc_se     = map_dbl(result, "mc_se")
  ) %>%
  select(model_type, ppi_type, n, delta, lambda_mode,
         PPI_crossfit, PPIpp_crossfit, PPIpp_independent_split,
         emp_power, mc_se) %>%
  mutate(
    design = case_when(
      ppi_type == "PPI"   &  PPI_crossfit ~ "PPI_crossfit",
      ppi_type == "PPI"   & !PPI_crossfit ~ "PPI_noncrossfit",
      ppi_type == "PPI++" &  PPIpp_crossfit ~ "PPI++_crossfit",
      ppi_type == "PPI++" &  PPIpp_independent_split ~ "PPI++_independent_split",
      ppi_type == "PPI++" & !PPIpp_crossfit & !PPIpp_independent_split ~ "PPI++_noncrossfit",
      TRUE ~ "unknown"
    ),
    default = case_when(
    # PPI default: noncrossfit
    ppi_type == "PPI"   & !PPI_crossfit               ~ TRUE,
    # PPI++ default: crossfit
    ppi_type == "PPI++" &  PPIpp_crossfit             ~ TRUE,
    # everything else is not the default design
    TRUE ~ FALSE
  )
  )

# All variance-covariance statistics
cov_var_all <- power_grid_all_comb %>%
  mutate(cov_tab = map(result, "cov_var_table")) %>%
  select(
    model_type, ppi_type, n, delta, lambda_mode,
    PPI_crossfit, PPIpp_crossfit, PPIpp_independent_split,
    cov_tab
  ) %>%
  unnest(cov_tab, names_sep = "_") %>%
  select(
    # keep the "outer" identifiers
    model_type, ppi_type, n, delta, lambda_mode,
    PPI_crossfit, PPIpp_crossfit, PPIpp_independent_split,
    # pull in only the useful inner summaries
    avg_lambda_hat      = cov_tab_avg_lambda_hat,
    sd_lambda_hat       = cov_tab_sd_lambda_hat,
    lambda_oracle       = cov_tab_lambda_oracle,
    var_fhat            = cov_tab_var_fhat,
    var_res             = cov_tab_var_res,
    cor_Y_fhat          = cov_tab_cor_Y_fhat,
    ratio_cov_over_var  = cov_tab_ratio_cov_over_var,
    N                   = cov_tab_N
  ) %>%
  mutate(
    design = case_when(
      ppi_type == "PPI"   &  PPI_crossfit ~ "PPI_crossfit",
      ppi_type == "PPI"   & !PPI_crossfit ~ "PPI_noncrossfit",
      ppi_type == "PPI++" &  PPIpp_crossfit ~ "PPI++_crossfit",
      ppi_type == "PPI++" &  PPIpp_independent_split ~ "PPI++_independent_split",
      ppi_type == "PPI++" & !PPIpp_crossfit & !PPIpp_independent_split ~ "PPI++_noncrossfit",
      TRUE ~ "unknown"
    ),
    default = case_when(
    # PPI default: noncrossfit
    ppi_type == "PPI"   & !PPI_crossfit               ~ TRUE,
    # PPI++ default: crossfit
    ppi_type == "PPI++" &  PPIpp_crossfit             ~ TRUE,
    # everything else is not the default design
    TRUE ~ FALSE
  )
  )

# # Raw file
# saveRDS(cov_var_all, file = "cov_var_all.rds")


############################################

# Main empirical power table
power_summary %>%
  select(model_type, ppi_type, design, n, delta, emp_power, mc_se) %>%
  arrange(model_type, ppi_type, design, n, delta) -> power_summary_table

# Detailed cov/var summary table
cov_var_all %>%
  select(model_type, ppi_type, design, n, delta,
         avg_lambda_hat, sd_lambda_hat, lambda_oracle,
         var_fhat, var_res, cor_Y_fhat) %>%
  arrange(model_type, ppi_type, design, n, delta) -> cov_var_summary_table

write.csv(power_summary_table, "power_summary.csv", row.names = FALSE)