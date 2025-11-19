library(tidyverse)
library(future)
library(furrr)
library(pppower)


devtools::load_all()


cl <- parallel::makeCluster(8)
plan(cluster, workers = cl)

parallel::clusterCall(cl, function() {
  library(pppower)
  library(tidyverse)
})

compute_theta0 <- function(a, X_sampler_L, f_generator, n_ref = 5e5) {
  Xref <- X_sampler_L(n_ref)
  Xref_m <- as.matrix(Xref)
  yref <- f_generator(Xref_m)

  H <- crossprod(Xref_m) / n_ref
  G <- crossprod(Xref_m, yref) / n_ref
  beta_star <- solve(H, G)

  as.numeric(sum(a * beta_star))
}


X_sampler_L <- function(n) data.frame(
  x1 = rnorm(n),
  x2 = rnorm(n)
)

f_true <- function(X) {
  X <- as.data.frame(X)
  with(X, 1 + 2*x1 - 1*x2 + 0.3*x1*x2 + 0.4*x2^2)
}

contrasts_list <- list(
  c1 = c(1, 0),     # β1
  c2 = c(0, 1),     # β2
  c3 = c(1, -1)     # β1 - β2
)

contrast_names <- names(contrasts_list)

theta0_list <- lapply(contrasts_list, compute_theta0,
                      X_sampler_L = X_sampler_L,
                      f_generator = f_true)
names(theta0_list) <- names(contrasts_list)
# Simulation grid

n_grid <- c(50, 100, 150, 200, 300)
N <- 2000
R <- 500

deltas <- seq(0, 0.6, by = 0.05)
#n_external <- c(50, 100, 200, 300, 500, 1000)

models <- c("glm_correct", "glm_mis", "glm_wrong")
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
    contrast     = as.character(contrast),
    model_type   = as.character(model_type),
    ppi_type     = as.character(ppi_type),
    lambda_mode  = as.character(lambda_mode)
  )

power_grid <- grid %>%
  mutate(
    result = future_pmap(
      list(contrast, model_type, delta, ppi_type, n, lambda_mode, lambda_external),
      ~ ppi_ols_empirical_power(
          R = R,
          n = ..5,
          N = N,
          X_sampler_L = X_sampler_L,
          f_generator = f_true,
          model_type  = ..2,
          a           = contrasts_list[[ as.character(..1) ]],
          delta       = ..3,
          theta0      = theta0_list[[ as.character(..1) ]],
          ppi_type    = ..4,
          lambda_mode = ..6,
          lambda_external = ..7,
          n_external  = ceiling(..5/2),
          PPIpp_crossfit = TRUE,
          alpha       = 0.05,
          seed        = 42,
          keep_reps = FALSE
      ),
      .options = furrr_options(
        seed = TRUE,
        packages = "pppower"
      ),
      .progress = TRUE
    )
  )

print("Finished gird loop. Saving .rds file...")

saveRDS(power_grid, "grid_results_ols.rds")

parallel::stopCluster(cl)