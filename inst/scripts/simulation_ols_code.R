library(tidyverse)
library(furrr)
library(dplyr)


# Data-generating process

X_sampler_L <- function(n) data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
)

# Nonlinear truth for OLS tests
f_true <- function(X) {
  X <- as.data.frame(X)
  x1 <- X[["x1"]]
  x2 <- X[["x2"]]
  x1 + x1 * x2 + x2^3
}

contrasts_list <- list(
  c1 = c(1, 0),     # β1
  c2 = c(0, 1),     # β2
  c3 = c(1, -1)     # β1 - β2
)

contrast_names <- names(contrasts_list)

# Simulation grid

n_grid <- c(50, 100, 150, 200, 300)
N <- 2000
R <- 500

deltas <- seq(0, 0.6, by = 0.05)
n_external <- c(50, 100, 200, 300, 500, 1000)

models <- c("glm_correct", "glm_mis", "glm_wrong", "rf")
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
)

plan(multisession, workers = 16, scheduling = TRUE)

power_grid <- grid %>%
  mutate(
    result = future_pmap(
      list(contrast, model_type, delta, ppi_type, n, lambda_mode, lambda_external),
      ~ ppi_ols_empirical_power(
          R = R,
          n = ..5,                       # labeled
          N = N,                         # unlabeled
          X_sampler_L = X_sampler_L,
          f_generator = f_true,
          model_type = ..2,
          a = contrasts_list[[..1]],
          delta = ..3,
          ppi_type = ..4,
          lambda_mode = ..6,
          lambda_external = ..7,
          n_external = ceiling(..5/2),   # consistent default
          PPIpp_crossfit = TRUE,
          alpha = 0.05,
          seed = 42,
          scheduling = Inf
      ),
      .options = furrr_options(seed = TRUE),
      .progress = TRUE
    )
  )

power_df <- power_grid %>%
  tidyr::unnest_wider(result) %>%
  select(
    contrast, model_type, ppi_type, n, delta,
    empirical_power, mc_se, avg_SE, lambda_external, avg_lambda_hat
  ) %>%
  mutate(
    ppi_type = factor(ppi_type, levels = c("PPI", "PPI++")),
    n = factor(n, levels = sort(unique(n))),
    model_type = factor(model_type, levels = models)
  )


lambda_summary <- power_df %>%
  group_by(contrast, model_type, ppi_type, n, lambda_external) %>%
  summarise(
    mean_lambda = mean(avg_lambda_hat, na.rm = TRUE),
    sd_lambda = sd(avg_lambda_hat, na.rm = TRUE),
    .groups = "drop"
  )


plot_df <- power_df %>%
  filter(ppi_type %in% c("PPI", "PPI++"))

for (ct in contrast_names) {
  dfc <- plot_df %>% filter(contrast == ct)

  p <- ggplot(
    dfc,
    aes(x = delta, y = empirical_power,
        color = model_type, group = model_type)
  ) +
    geom_line(linewidth = 1.1) +
    geom_point(size = 1.8) +
    facet_grid(ppi_type ~ n, scales = "free_y") +
    labs(
      title = paste("OLS Contrast Power Curve — contrast:", ct),
      x = expression(Delta),
      y = "Empirical Power"
    ) +
    theme_bw(13)

  print(p)

}