library(tidyverse)
library(furrr)
library(dplyr)
library(future)

devtools::load_all()

X_sampler_L <- function(n) data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n),
  x3 = rnorm(n)
)

# Linear predictor: η(x)
f_true <- function(X) {
  X <- as.data.frame(X)
  with(X, 1 + 2 * x1 + 3.5 * x2 - x3)
}

# Contrasts for GLM coefficients β = (β1, β2, β3)
contrasts_list <- list(
  c1 = c(1, 0, 0),      # test β1
  c2 = c(0, 1, 0),      # test β2
  c3 = c(1, 0, -1)      # test β1 – β3
)

contrast_names <- names(contrasts_list)

# Compute θ₀ = cᵀ β*
theta0_list <- lapply(
  contrasts_list, compute_theta0,
  X_sampler_L = X_sampler_L,
  f_generator = f_true
)
names(theta0_list) <- contrast_names

## Simulation grid

n_grid <- c(50, 100, 150, 200, 300)
N      <- 2000
R      <- 500

deltas <- seq(0, 0.6, by = 0.05)

models    <- c("glm_correct", "glm_mis", "glm_wrong")
ppi_types <- c("PPI", "PPI++")

grid <- expand.grid(
  contrast = contrast_names,
  delta = deltas,
  model_type = models,
  ppi_type = ppi_types,
  n = n_grid,
  lambda_mode = "plugin",
  lambda_external = TRUE,
  stringsAsFactors = FALSE
) %>%
  mutate(
    contrast   = as.character(contrast),
    model_type = as.character(model_type),
    ppi_type   = as.character(ppi_type)
  )


## Parallel computation setup

plan(multicore, workers = 8)


## Run GLM simulation

power_grid <- grid %>%
  mutate(
    result = future_pmap(
      list(contrast, model_type, delta, ppi_type, n, lambda_mode, lambda_external),
      ~ ppi_glm_empirical_power(
          R = R,
          n = ..5,
          N = N,
          X_sampler_L = X_sampler_L,
          f_generator = f_true,
          model_type   = ..2,
          a            = contrasts_list[[as.character(..1)]],
          delta        = ..3,
          theta0       = theta0_list[[as.character(..1)]],
          ppi_type     = ..4,
          lambda_mode  = ..6,
          lambda_external = ..7,
          n_external   = ceiling(..5/2),
          PPIpp_crossfit = TRUE,
          alpha = 0.05,
          seed = 42
      ),
      .options = furrr_options(seed = TRUE, packages = "pppower"),
      .progress = TRUE
    )
  )


## Expand results into a dataframe

power_df <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  select(
    contrast, model_type, ppi_type, n, delta,
    empirical_power, mc_se,
    avg_se_hat,       # empirical SE
    se_theory,        # theoretical SE
    theoretical_power,
    lambda_external, avg_lambda_hat
  ) %>%
  mutate(
    ppi_type   = factor(ppi_type, levels = c("PPI", "PPI++")),
    model_type = factor(model_type, levels = models),
    n          = factor(n, levels = sort(unique(n)))
  )


lambda_summary <- power_df %>%
  group_by(contrast, model_type, ppi_type, n, lambda_external) %>%
  summarise(
    mean_lambda = mean(avg_lambda_hat, na.rm = TRUE),
    sd_lambda = sd(avg_lambda_hat, na.rm = TRUE),
    .groups = "drop"
  )

plot_contrast_power_glm <- function(ct, plot_df, models) {
  # Filter for selected contrast
  dfc <- plot_df %>%
    filter(contrast == ct) %>%
    mutate(
      ppi_type = factor(ppi_type, levels = c("PPI", "PPI++")),
      model_type = factor(model_type, levels = models),
      n = factor(n, levels = sort(unique(n)))
    )
  
  # Compute PPI++ - PPI differences 
  dfc_diff <- dfc %>%
    select(delta, model_type, n, ppi_type, empirical_power) %>%
    tidyr::pivot_wider(
      names_from = ppi_type,
      values_from = empirical_power
    ) %>%
    mutate(power_diff = `PPI++` - PPI)
  
  # Convert to long format with 3 panels
  dfc_long <- bind_rows(
    dfc %>% mutate(panel = ppi_type),
    dfc_diff %>% mutate(ppi_type = NA, empirical_power = NA, panel = "Difference")
  )
  
  dfc_long$panel <- factor(
    dfc_long$panel,
    levels = c("PPI", "PPI++", "Difference")
  )
  
  # Panel-specific y-axis range
  dfc_long <- dfc_long %>%
    mutate(
      yplot = case_when(
        panel == "Difference" ~ power_diff,
        TRUE ~ empirical_power
      ),
      ymin_fake = case_when(
        panel == "Difference" ~ -0.5,
        TRUE ~ 0
      ),
      ymax_fake = case_when(
        panel == "Difference" ~ 0.5,
        TRUE ~ 1
      )
    )
  
  p <- ggplot(
    dfc_long,
    aes(
      x = delta,
      y = yplot,
      color = model_type,
      group = model_type
    )
  ) +
    geom_blank(aes(y = ymin_fake)) +
    geom_blank(aes(y = ymax_fake)) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.8) +
    facet_grid(panel ~ n, scales = "free_y") +
    labs(
      title = paste0("GLM Logistic Contrast Power Curve — Contrast ", ct),
      subtitle = "Rows: PPI, PPI++, and (PPI++ - PPI)",
      x = expression(Delta),
      y = "Empirical Power / Difference",
      color = "Model"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
  
  return(p)
}

plot_df_glm <- power_df_glm %>%
  filter(ppi_type %in% c("PPI", "PPI++"))


plot_contrast_power_glm("c1", plot_df_glm, models)
plot_contrast_power_glm("c2", plot_df_glm, models)
plot_contrast_power_glm("c3", plot_df_glm, models)


emp_theoretical_power <- ggplot(
  power_df,
  aes(
    x = theoretical_power,
    y = empirical_power,
    color = model_type
  )
) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  geom_point(size = 2, alpha = 0.6) +
  facet_grid(ppi_type ~ n + contrast, scales = "free") +
  labs(
    title = "Empirical vs Analytical Power (GLM Contrast Tests)",
    subtitle = "Rows: Estimator type (PPI / PPI++), Columns: n × contrast",
    x = "Theoretical Power",
    y = "Empirical Power",
    color = "Model",
    caption = "Closed-form GLM variance; λ from plugin/external as specified."
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

print(emp_theoretical_power)
