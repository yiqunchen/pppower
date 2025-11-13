
#' One Monte Carlo replicate for PPI / PPI++
#'
#' @description
#' Internal engine for computing one repetition of (theta_hat, lambda_hat, SE, reject).
#'
#' @keywords internal
#'
#' @return A list with theta_hat, se, z, reject, lambda_hat, and variance components.
#'
#' @importFrom stats qnorm var cov rnorm
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

#' Monte Carlo Estimation of Empirical and Close-Form Power for PPI / PPI++ mean estimation
#'
#' @description
#' Runs a full Monte Carlo simulation for estimating power for the
#' naive, PPI, or PPI++ estimators. Returns empirical rejection probability,
#' Monte Carlo standard error, and theoretical (closed-form) power based
#' on estimated population-level covariance components.
#'
#' @param R Number of Monte Carlo repetitions.
#' @param n Size of labeled sample.
#' @param N Size of unlabeled sample.
#' @param X_sampler_L Function generating labeled covariates.
#' @param X_sampler_U Function generating unlabeled covariates (defaults to `X_sampler_L`).
#' @param f_generator Function generating signal mean function f(x).
#' @param eps_sampler Error distribution sampler, default standard normal.
#' @param model_type Choice of predictive model: `"glm_correct"`, `"glm_mis"`, `"glm_wrong"`, `"rf"`.
#' @param ppi_type `"naive"`, `"PPI"`, or `"PPI++"`.
#' @param lambda_mode `"plugin"` or `"theory"`.
#' @param lambda_external Whether to estimate lambda on an external set.
#' @param n_external Size of external set if used.
#' @param alpha Significance level.
#' @param delta True mean shift.
#' @param seed Optional integer random seed for reproducibility.
#' @param pred_noise Optional added prediction noise.
#' @param PPI_crossfit Logical, use cross-fitting for PPI.
#' @param PPIpp_crossfit Logical, use cross-fitting for PPI++.
#' @param PPIpp_independent_split Logical, use 3-way independent split.
#'
#' @return A list containing:
#' \item{empirical_power}{Monte Carlo rejection rate}
#' \item{mc_se}{Monte Carlo SE of empirical power}
#' \item{theoretical_power}{closed-form population power}
#' \item{avg_lambda_hat}{average estimated \eqn{\hat\lambda}}
#' \item{cov_var_table}{tibble with covariance summaries used for power formula}
#'
#' @export
#'
#' @examples
#' ppi_empirical_power(
#'   R = 500,
#'   n = 100,
#'   N = 10000,
#'   X_sampler_L = function(n) data.frame(x1 = rnorm(n), x2 = rnorm(n)),
#'   f_generator = function(X) 2 * X$x1,
#'   model_type = "glm_correct",
#'   ppi_type = "PPI"
#' )
#' @importFrom stats qnorm rnorm cov var predict
#' @importFrom tibble tibble
#' 
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
