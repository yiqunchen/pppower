library(furrr)     # Required
library(purrr)     # Required
library(tidyr)     # Required
library(tibble)    # Required
library(dplyr)     # Not required if rewrite code
library(ggplot2)   # Suggested for tutorial & vignettes
library(ggforce)   # Suggested for tutorial & vignettes
library(patchwork) # Suggested for tutorial & vignettes


################################### Simulation code from start to end

# Initialize simulated data for all models
simulate_one_draw <- function(n_labeled, 
                              n_unlabeled,
                              X_sampler_L, 
                              X_sampler_U = NULL,
                              f_generator,
                              eps_sampler = function(n) rnorm(n, 0, 1),
                              delta = 0) {
  if (is.null(X_sampler_U)) X_sampler_U <- X_sampler_L
 
  X_L <- X_sampler_L(n_labeled)
  X_U <- X_sampler_U(n_unlabeled)

  f_L_base <- f_generator(X_L)
  f_U_base <- f_generator(X_U)

  # mean shift by delta
  f_L <- f_L_base + delta
  f_U <- f_U_base + delta

  y_L <- f_L + eps_sampler(n_labeled)

  list(
    labeled   = list(X = X_L, y = y_L, f = f_L),
    unlabeled = list(X = X_U, f = f_U),
    theta_true = mean(f_L),   # population mean of Y is ~ mean(f) here
    delta = delta
  )
}

# Wrapper for predicting the outcome for all models
.predict_any <- function(fit, newdata) {
  if (inherits(fit, "ranger")) {
    predict(fit, data = newdata)$predictions
  } else {
    stats::predict(fit, newdata = newdata)
  }
}

# Train and predict
fit_predict_model <- function(model_type, X_L, y_L, X_U,
                              mtry = NULL,
                              rf_engine = c("ranger", "randomForest"),
                              rf_trees = 200,
                              rf_min_node_size = 5,
                              rf_num_threads = NULL,
                              rf_seed = NULL) {
  model_type <- match.arg(model_type, c("glm_correct", "glm_mis", "glm_wrong", "rf"))
  rf_engine  <- match.arg(rf_engine, c("ranger", "randomForest"))
  dat_L <- data.frame(y = y_L, X_L)

  if (model_type == "glm_correct") {
    fit <- stats::glm(y ~ x1 + I(x1 * x2) + I(x2^3), data = dat_L)
    fhat_U <- stats::predict(fit, newdata = X_U)

  } else if (model_type == "glm_mis") {
    fit <- stats::glm(y ~ x1 + x2, data = dat_L)
    fhat_U <- stats::predict(fit, newdata = X_U)

  } else if (model_type == "glm_wrong") {
    fit <- stats::glm(y ~ I(x1 * x2), data = dat_L)
    fhat_U <- stats::predict(fit, newdata = X_U)

  } else if (model_type == "rf") {
    mtry_eff  <- if (is.null(mtry)) floor(sqrt(ncol(X_L))) else mtry
    use_rgr   <- (rf_engine == "ranger") && requireNamespace("ranger", quietly = TRUE)
    if (use_rgr) {
      fit <- ranger::ranger(
        y ~ ., data = dat_L,
        num.trees     = rf_trees,
        mtry          = mtry_eff,
        min.node.size = rf_min_node_size,
        importance    = "none",
        seed          = rf_seed,
        num.threads   = rf_num_threads
      )
      fhat_U <- predict(fit, data = X_U)$predictions
    } else {
      fit <- randomForest::randomForest(y ~ ., data = dat_L, ntree = rf_trees, mtry = mtry_eff)
      fhat_U <- stats::predict(fit, newdata = X_U)
    }
  }

  list(fit = fit, fhat_U = as.numeric(fhat_U))
}

# One Monte Carlo replicate of the inference procedure for a given configuration
ppi_mean_rep_one <- function(
  n, N,
  X_sampler_L, X_sampler_U = NULL,
  f_generator,
  eps_sampler = function(n) rnorm(n, 0, 1),
  model_type = c("glm_correct", "glm_mis", "glm_wrong", "rf"),
  ppi_type = c("naive", "PPI", "PPI++"),
  lambda_mode = c("plugin", "theory"),  # theory: uses oracle lambda
  lambda_external = FALSE,              # lambda estimated on an external labeled set of data
  n_external = ceiling(n/2),            # if above is TRUE, size of that labeled data
  alpha = 0.05,
  delta = 0,
  seed = NULL,
  pred_noise = 0,
  PPI_crossfit = FALSE,                 # Default flag, controls Xfit or not for PPI
  PPIpp_crossfit = TRUE,                # Default flag, controls Xfit or not for PPI++
  PPIpp_independent_split = FALSE       # Disabled
) {
  model_type <- match.arg(model_type)
  ppi_type   <- match.arg(ppi_type)
  if (!is.null(seed)) set.seed(seed)
  if (is.null(X_sampler_U)) X_sampler_U <- X_sampler_L

  sim <- simulate_one_draw(
    n_labeled   = n,
    n_unlabeled = N,
    X_sampler_L = X_sampler_L,
    X_sampler_U = X_sampler_U,
    f_generator = f_generator,
    eps_sampler = eps_sampler,
    delta       = delta
  )
  X_L <- sim$labeled$X
  y_L <- sim$labeled$y
  X_U <- sim$unlabeled$X

  ## Initialize
  fhat_l <- fhat_u <- rep(0, 0)

  ## CASE 1: NAIVE ESTIMATOR
  if (ppi_type == "naive") {
    fhat_l <- rep(0, n)
    fhat_u <- rep(0, N)

  ## CASE 2: PPI
  } else if (ppi_type == "PPI") {
    if (PPI_crossfit) {
      ## Cross-fit version
      fold_id <- sample(rep(1:2, length.out = n))
      fhat_l <- numeric(n)
      fhat_u_fold1 <- numeric(N)
      fhat_u_fold2 <- numeric(N)
      for (fold in 1:2) {
        tr <- which(fold_id != fold)
        te <- which(fold_id == fold)
        m <- fit_predict_model(model_type, X_L[tr, , drop = FALSE], y_L[tr], X_U)
        fhat_l[te] <- as.numeric(.predict_any(m$fit, X_L[te, , drop = FALSE]))
        assign(paste0("fhat_u_fold", fold), m$fhat_U)
      }
      fhat_u <- (fhat_u_fold1 + fhat_u_fold2) / 2
    } else {
      ## Independent training for PPI (default)
      n_train <- n / 2
      sim_train <- simulate_one_draw(
        n_labeled = n_train,
        n_unlabeled = 0,
        X_sampler_L = X_sampler_L,
        f_generator = f_generator,
        eps_sampler = eps_sampler,
        delta = delta
      )
      m_full <- fit_predict_model(model_type, sim_train$labeled$X, sim_train$labeled$y, X_U)
      fhat_l <- as.numeric(.predict_any(m_full$fit, X_L))
      fhat_u <- m_full$fhat_U
    }

  ## CASE 3: PPI++
  } else if (ppi_type == "PPI++") {
    if (PPIpp_crossfit) {
      ## Cross-fitting version (default)
      fold_id <- sample(rep(1:2, length.out = n))
      fhat_l <- numeric(n)
      fhat_u_fold1 <- numeric(N)
      fhat_u_fold2 <- numeric(N)
      for (fold in 1:2) {
        tr <- which(fold_id != fold)
        te <- which(fold_id == fold)
        m <- fit_predict_model(model_type, X_L[tr, , drop = FALSE], y_L[tr], X_U)
        fhat_l[te] <- as.numeric(.predict_any(m$fit, X_L[te, , drop = FALSE]))
        assign(paste0("fhat_u_fold", fold), m$fhat_U)
      }
      fhat_u <- (fhat_u_fold1 + fhat_u_fold2) / 2

    } else if (PPIpp_independent_split) {
      ## Independent-split version (3 disjoint parts)
      n_fit <- ceiling(n / 2)  # size of training set for f
      sim_fit <- simulate_one_draw(
        n_labeled = n_fit,
        n_unlabeled = 0,
        X_sampler_L = X_sampler_L,
        f_generator = f_generator,
        eps_sampler = eps_sampler,
        delta = delta
      )

      # Fit f on independent training data (X_f, Y_f)
      m_fit <- fit_predict_model(model_type, sim_fit$labeled$X, sim_fit$labeled$y, X_U)

      # Predictions on labeled and unlabeled sets (f fixed)
      fhat_l <- as.numeric(.predict_any(m_fit$fit, X_L))
      fhat_u <- m_fit$fhat_U

    } else {
      ## Non-crossfit, shared model (trained on labeled data itself)
      m_full <- fit_predict_model(model_type, X_L, y_L, X_U)
      fhat_l <- as.numeric(.predict_any(m_full$fit, X_L))
      fhat_u <- m_full$fhat_U
    }
  }

  ## Optionally add prediction noise, disabled currently
  if (pred_noise > 0) {
    fhat_l <- fhat_l + rnorm(length(fhat_l), 0, pred_noise)
    fhat_u <- fhat_u + rnorm(length(fhat_u), 0, pred_noise)
  }

  ## Estimation & Inference
  res_labeled <- y_L - fhat_l
  sigma_eff_sq <- mean(res_labeled^2)
  sigma_eff <- sqrt(sigma_eff_sq)
  var_f_hat <- var_res_hat <- NA_real_

  if (ppi_type == "naive") {
    theta_hat <- mean(y_L)
    var_res_hat <- var(y_L)
    se_hat <- sqrt(var_res_hat / n)
    lambda_hat <- NA

  } else if (ppi_type == "PPI") {
    theta_hat <- mean(fhat_u) - mean(fhat_l - y_L)
    var_f_hat   <- var(fhat_u)
    var_res_hat <- var(res_labeled)
    se_hat <- sqrt(var_f_hat / N + var_res_hat / n)
    lambda_hat <- 1

  } else if (ppi_type == "PPI++") {
    cov_yf <- cov(y_L, fhat_l)
    var_f  <- var(c(fhat_l, fhat_u)) # Computed from the pooled data set (labeled & unlabeled)
    if (lambda_mode == "theory") {
      lambda_hat <- 1 / (1 + n / N)
    } else if (exists("lambda_external") && isTRUE(lambda_external)) {
    # External lambda sub-branch:
    # Draw an independent external labeled dataset of size n/2
    n_ext <- n_external
    sim_ext <- simulate_one_draw(
      n_labeled   = n_ext,
      n_unlabeled = 0,
      X_sampler_L = X_sampler_L,
      f_generator = f_generator,
      eps_sampler = eps_sampler,
      delta       = delta
    )

    X_ext <- sim_ext$labeled$X
    y_ext <- sim_ext$labeled$y

    # Predict fhat for the external set using model trained on the main labeled data
    fit_ext <- fit_predict_model(
      model_type = model_type,
      X_L = X_L,
      y_L = y_L,
      X_U = X_ext
    )
    fhat_ext <- fit_ext$fhat_U

    # Estimate lambda from external set
    cov_yf_ext <- cov(y_ext, fhat_ext)
    var_f_ext  <- var(fhat_ext) # Computed only based on the external labeled set
    lambda_raw <- cov_yf_ext / ((1 + n / N) * var_f_ext)
    lambda_hat <- min(1, max(0, lambda_raw))
    } else {
      # Normal estimation of lambda
      lambda_raw <- cov_yf / ((1 + n / N) * var_f)
      lambda_hat <- min(1, max(0, lambda_raw))
    }
    theta_hat <- mean(y_L) + lambda_hat * (mean(fhat_u) - mean(fhat_l))
    var_res_lambda <- var(y_L - lambda_hat * fhat_l)
    se_hat <- sqrt(var_res_lambda / n + (lambda_hat^2) * var(fhat_u) / N)
    var_f_hat <- var(fhat_u)
    var_res_hat <- var_res_lambda
  }

  ## Inference test
  z_alpha <- qnorm(1 - alpha / 2)
  z_stat  <- (theta_hat - 0) / se_hat
  reject  <- is.finite(z_stat) && abs(z_stat) > z_alpha   # Rejection status

  mse_l <- mean((y_L - fhat_l)^2)
  r2_l  <- 1 - mse_l / var(y_L)
  cor_l <- suppressWarnings(cor(y_L, fhat_l))

  list(
    theta_hat = theta_hat,
    se = se_hat,
    z = z_stat,
    reject = reject,
    lambda_hat = if (ppi_type == "PPI++") lambda_hat else 1,
    sigma_eff = sigma_eff,
    components = list(var_f = var_f_hat, var_res = var_res_hat),
    mse = mse_l,
    r2 = r2_l,
    cor = cor_l
  )
}

# Run Monte Carlo for a fixed number of times and output 
ppi_empirical_power <- function(R = 2000,
                                n, N,
                                X_sampler_L, X_sampler_U = NULL,
                                f_generator,
                                eps_sampler = function(n) rnorm(n, 0, 1),
                                model_type = c("glm_correct", "glm_mis", "glm_wrong", "rf"),
                                ppi_type = c("naive", "PPI", "PPI++"),
                                lambda_mode = c("plugin", "theory"),
                                lambda_external = FALSE,
                                n_external = ceiling(n/2),
                                alpha = 0.05,
                                delta = 0,
                                seed = NULL,
                                pred_noise = 0,
                                PPI_crossfit = FALSE,
                                PPIpp_crossfit = TRUE,         # Default settings
                                PPIpp_independent_split = FALSE) {
  model_type  <- match.arg(model_type)
  ppi_type    <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  if (!is.null(seed)) set.seed(seed)

  ## Monte Carlo repetitions 
  reps <- replicate(
    R,
    ppi_mean_rep_one(
      n = n, N = N,
      X_sampler_L = X_sampler_L,
      X_sampler_U = X_sampler_U,
      f_generator = f_generator,
      eps_sampler = eps_sampler,
      model_type  = model_type,
      ppi_type    = ppi_type,
      alpha       = alpha,
      delta       = delta,
      lambda_mode = lambda_mode,
      lambda_external = lambda_external,
      n_external = n_external,
      pred_noise  = pred_noise,
      PPI_crossfit = PPI_crossfit,
      PPIpp_crossfit = PPIpp_crossfit,
      PPIpp_independent_split = PPIpp_independent_split
    ),
    simplify = FALSE
  )

  ## Extract main summaries 
  rej         <- vapply(reps, `[[`, logical(1), "reject")
  se_avg      <- mean(vapply(reps, function(x) x$se, numeric(1)))
  theta_avg   <- mean(vapply(reps, function(x) x$theta_hat, numeric(1)))
  emp_power   <- mean(rej)
  mc_se       <- sqrt(emp_power * (1 - emp_power) / R)

  ## Lambda summaries
  lambda_hats <- vapply(reps, function(x) x$lambda_hat, numeric(1))
  lambda_avg  <- mean(lambda_hats, na.rm = TRUE)
  lambda_sd   <- sd(lambda_hats, na.rm = TRUE)

  ## Cov–Var summaries based on MC runs, later plug into the close-form equation for power calculation
  var_fhat   <- mean(vapply(reps, function(x) x$components$var_f,   numeric(1)), na.rm = TRUE)
  var_res    <- mean(vapply(reps, function(x) x$components$var_res, numeric(1)), na.rm = TRUE)
  cor_y_fhat <- mean(vapply(reps, function(x) x$cor, numeric(1)), na.rm = TRUE)

  ## Approximate Cov(Y,f)/Var(f)
  ratio_cov_over_var <- mean(vapply(reps, function(x) {
    cval <- x$cor
    vf   <- x$components$var_f
    vr   <- x$components$var_res
    if (is.na(vf) || is.na(cval)) return(NA_real_)
    sd_y <- sqrt(vf + vr)
    cval * sd_y / sqrt(vf)
  }, numeric(1)), na.rm = TRUE)

  ## Oracle lambda based on sample-derived moments
  lambda_oracle <- ratio_cov_over_var * (1 / (1 + n / N))

  ## Identify the inference design 
  design_label <- dplyr::case_when(
    ppi_type == "naive" ~ "naive",
    ppi_type == "PPI" &  PPI_crossfit                 ~ "PPI_crossfit",
    ppi_type == "PPI" & !PPI_crossfit                 ~ "PPI_independent",
    ppi_type == "PPI++" & PPIpp_crossfit              ~ "PPI++_crossfit",
    ppi_type == "PPI++" & PPIpp_independent_split     ~ "PPI++_independent_split",
    ppi_type == "PPI++" & !PPIpp_crossfit & !PPIpp_independent_split ~ "PPI++_shared",
    TRUE ~ "unknown"
  )

  z_alpha <- qnorm(1 - alpha / 2)
  se_theory <- dplyr::case_when(
    ppi_type == "PPI"   ~ sqrt(var_fhat / N + var_res / n),
    ppi_type == "PPI++" ~ sqrt(var_res / n + (lambda_oracle^2) * var_fhat / N),
    TRUE ~ NA_real_
  )

  th_power <- 1 - pnorm(z_alpha - delta / se_theory) +
                  pnorm(-z_alpha - delta / se_theory)

  ## Tidy table
  cov_var_table <- tibble::tibble(
    model_type         = model_type,
    ppi_type           = ppi_type,
    lambda_mode        = lambda_mode,
    design             = design_label,
    n                  = n,
    N                  = N,
    var_fhat           = var_fhat,
    var_res            = var_res,
    cor_Y_fhat         = cor_y_fhat,
    ratio_cov_over_var = ratio_cov_over_var,
    avg_lambda_hat     = lambda_avg,
    sd_lambda_hat      = lambda_sd,
    lambda_oracle      = lambda_oracle,
    se_theory          = se_theory,
    power_theory       = th_power
  )


  ## Return results
  list(
    empirical_power = emp_power,
    theoretical_power = th_power,
    mc_se           = mc_se,
    avg_SE          = se_avg,
    avg_theta_hat   = theta_avg,
    avg_lambda_hat  = lambda_avg,
    sd_lambda_hat   = lambda_sd,
    n_external      = n_external,
    cov_var_table   = cov_var_table,
    details         = reps
  )
}


#############################

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

# write.csv(power_summary_table, "power_summary.csv", row.names = FALSE)

############################################ Check for lambda hats, oracle lambda, and theoretical lambda

