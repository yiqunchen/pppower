#' One Monte Carlo Replicate for OLS PPI / PPI++
#'
#' @description
#' Generates one Monte Carlo replication of the OLS–based PPI or PPI++
#' estimator for a scalar linear contrast \eqn{\\theta = a^\top \\beta}.
#'
#' This function returns all intermediate objects needed for empirical and
#' closed-form power calculations and is intended for internal Monte Carlo use.
#'
#'
#' @keywords internal
#'
#' @param n Number of labeled observations.
#' @param N Number of unlabeled observations.
#' @param X_sampler_L Function generating covariate matrices (\code{n x p}).
#' @param f_generator Function generating the regression signal \eqn{f(X)}.
#' @param eps_sampler Noise generator, default i.i.d. standard normal.
#' @param model_type Predictive model type ("glm_correct", "glm_mis",
#'   "glm_wrong", "rf").
#' @param a Contrast vector of length \eqn{p}.
#' @param delta Effect size injected via \eqn{f(X) + \delta (X^\top a)}.
#' @param theta0 Null value of the contrast \eqn{\theta=a^\top\beta^\star}.
#' @param ppi_type Either "PPI" ({\eqn{\lambda = 1}}) or "PPI++".
#' @param lambda_mode For PPI++: "plugin" (estimate {\eqn{\lambda}})
#'   or "oracle" (not implemented internally).
#' @param lambda_external Logical; if TRUE, estimate {\eqn{\lambda}}
#'   using an external labeled sample.
#' @param n_external Size of the external labeled set if
#'   \code{lambda_external = TRUE}.
#' @param PPIpp_crossfit Logical; if TRUE, use 2-fold cross-fitting
#'   for PPI++.
#' @param alpha Wald test significance level.
#' @param seed Optional random seed.
#'
#' @return A list with:
#' \item{theta_hat}{Estimated contrast \eqn{a^\top \hat\beta}.}
#' \item{se_hat}{Estimated standard error.}
#' \item{reject}{Logical Wald rejection indicator.}
#' \item{lambda_hat}{Shrinkage parameter used.}
#' \item{beta_hat}{Estimated regression coefficients.}
#' \item{H_L}{Labeled second-moment matrix.}
#' \item{H_U}{Unlabeled second-moment matrix.}
#' \item{Sigma_YY_hat}{Covariance of residual-score terms.}
#' \item{Sigma_ff_hat}{Covariance of prediction-error terms.}
#' \item{Sigma_Yf_hat}{Covariance of residual/prediction interactions.}
#' 
ppi_ols_rep_one <- function(
  n, N,
  X_sampler_L,
  f_generator,
  eps_sampler = function(n) rnorm(n, 0, 1),
  model_type = c("glm_correct", "glm_mis", "glm_wrong", "rf"),
  a,
  delta = 0,
  ppi_type = c("PPI", "PPI++"),
  theta0 = 0,
  lambda_mode = c("plugin", "oracle"),
  lambda_external = FALSE,
  n_external = ceiling(n/2),
  PPIpp_crossfit = TRUE,
  alpha = 0.05,
  seed = NULL
) {
  
  model_type <- match.arg(model_type)
  ppi_type   <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  if (!is.null(seed)) set.seed(seed)

  a <- as.numeric(a)

  X_L_df <- X_sampler_L(n)  # data.frame with x1, x2
  X_U_df <- X_sampler_L(N)

  X_L <- as.matrix(as.data.frame(X_L_df))
  X_U <- as.matrix(as.data.frame(X_U_df))
  mode(X_L) <- "numeric"
  mode(X_U) <- "numeric"

  p <- ncol(X_L)
  stopifnot(length(a) == p)

  f_L_true <- f_generator(X_L_df)
  f_U_true <- f_generator(X_U_df)

  # inject contrast-aligned signal shift
  f_L <- f_L_true + delta * as.numeric(X_L %*% a)
  #f_U <- f_U_true + delta * as.numeric(X_U %*% a)

  y_L <- f_L + eps_sampler(n)

  if (ppi_type == "PPI++" && PPIpp_crossfit) {

    fold_id <- sample(rep(1:2, length.out = n))
    fhat_L <- numeric(n)
    fhat_U_fold1 <- numeric(N); fhat_U_fold2 <- numeric(N)

    fhat_U_folds <- vector("list", 2)

    for (fold in 1:2) {
      tr <- which(fold_id != fold)
      te <- which(fold_id == fold)

      m <- fit_predict_model_ols(
        model_type,
        X_L_df[tr, , drop = FALSE], y_L[tr],
        X_U_df
      )

      fhat_L[te] <- as.numeric(.predict_any(m$fit, X_L[te, , drop=FALSE]))
      fhat_U_folds[[fold]] <- m$fhat_U
    }

    # average U predictions
    fhat_U <- (fhat_U_folds[[1]] + fhat_U_folds[[2]]) / 2

  } else {
    ## Non-cross-fit model (PPI or PPI++)
    m <- fit_predict_model_ols(model_type, X_L_df, y_L, X_U_df)
    fhat_L <- as.numeric(.predict_any(m$fit, X_L_df))
    fhat_U <- m$fhat_U
  }

  ## Compute moment matrices H_L, H_U, G_L, G_U, W (correct PPI++ form)
  H_L <- crossprod(X_L, X_L) / n
  H_U <- crossprod(X_U, X_U) / N

  G_L <- crossprod(X_L, y_L) / n
  G_U <- crossprod(X_U, fhat_U) / N
  W_L <- crossprod(X_L, f_L_true) / n   # population "W" analog

  ## Compute lambda for PPI++
  if (ppi_type == "PPI") {

    lambda_hat <- 1

  } else {

    if (lambda_external) {

      ## External labeled sample
      sim_ext <- simulate_one_draw_contrast(
        n_labeled = n_external,
        X_sampler_L = X_sampler_L,
        f_generator = f_generator,
        eps_sampler = eps_sampler,
        a = a,
        delta = delta
      )
      
      X_ext_df <- sim_ext$X
      y_ext    <- sim_ext$y

      X_ext <- as.matrix(as.data.frame(X_ext_df))
      mode(X_ext) <- "numeric"

      ## Predict f on external data using model trained on main X_L,y_L
      m_ext <- fit_predict_model_ols(model_type, X_L_df, y_L, X_ext_df)
      fhat_ext <- m_ext$fhat_U

      ## External score covariances
      resid_ext <- y_ext - fhat_ext
      XY_ext <- sweep(X_ext, 1, resid_ext, `*`)
      Sigma_YY_ext <- cov(XY_ext)

      lin_ext <- X_ext %*% solve(H_U, G_U)
      ff_ext <- (fhat_ext - lin_ext)
      X_ff_ext <- sweep(X_ext, 1, ff_ext, `*`)
      Sigma_ff_ext <- cov(X_ff_ext)
      Sigma_Yf_ext <- cov(XY_ext, X_ff_ext)

      HinvU <- solve(H_U)
      B_ext <- t(a) %*% HinvU %*% Sigma_ff_ext %*% HinvU %*% a
      C_ext <- t(a) %*% HinvU %*% Sigma_Yf_ext %*% HinvU %*% a

      r <- n / N
      lambda_hat <- as.numeric(C_ext / ((1+r) * B_ext))
      lambda_hat <- max(0, min(1, lambda_hat))

    } else if (lambda_mode == "plugin") {

      ## Plugin lambda using main-labeled covariances
      resid_L <- y_L - fhat_L
      XY_resid <- sweep(X_L, 1, resid_L, `*`)
      Sigma_YY <- cov(XY_resid)

      lin_L <- X_L %*% solve(H_U, G_U)
      ff_L <- (fhat_L - lin_L)
      X_ff_L <- sweep(X_L, 1, ff_L, `*`)
      Sigma_ff <- cov(X_ff_L)
      Sigma_Yf <- cov(XY_resid, X_ff_L)

      HinvU <- solve(H_U)
      B <- t(a) %*% HinvU %*% Sigma_ff %*% HinvU %*% a
      C <- t(a) %*% HinvU %*% Sigma_Yf %*% HinvU %*% a

      r <- n / N
      lambda_hat <- as.numeric(C / ((1+r)*B))
      lambda_hat <- max(0, min(1, lambda_hat))

    } else {
      stop("Oracle lambda not implemented internally.")
    }
  }

  # PPI uses only H_U and (G_U + correction)
  if (ppi_type == "PPI") {
    G_p <- G_U + (crossprod(X_L, (y_L - fhat_L)) / n)
    H_p <- H_U
  } else {
    lambda <- lambda_hat
    # Correct PPI++ estimating equation:
    # betahat_lambda = ((1-λ)H_L + λH_U)^(-1) [G_L + λ(G_U - W_L)]
    H_p <- (1 - lambda) * H_L + lambda * H_U
    G_p <- G_L + lambda * (G_U - W_L)
  }

  beta_hat <- solve(H_p, G_p)
  theta_hat <- sum(a * beta_hat)

  ## Variance components (sample analogues)
  resid_L <- y_L - as.numeric(X_L %*% beta_hat)
  XY_resid <- sweep(X_L, 1, resid_L, `*`)
  Sigma_YY <- cov(XY_resid)

  ff_L <- (fhat_L - as.numeric(X_L %*% beta_hat))
  X_ff <- sweep(X_L, 1, ff_L, `*`)
  Sigma_ff <- cov(X_ff)

  Sigma_Yf <- cov(XY_resid, X_ff)

  ## Sandwich SE for PPI / PPI++
  Hinv <- solve(H_p)
  A_scalar <- as.numeric(t(a) %*% Hinv %*% Sigma_YY %*% Hinv %*% a)
  B_scalar <- as.numeric(t(a) %*% Hinv %*% Sigma_ff %*% Hinv %*% a)
  C_scalar <- as.numeric(t(a) %*% Hinv %*% Sigma_Yf %*% Hinv %*% a)

  Vhat <- A_scalar / n +
    lambda_hat^2 * (B_scalar / N + B_scalar / n) -
    2 * lambda_hat * C_scalar / n

  se_hat <- sqrt(as.numeric(Vhat))

  z_stat <- (theta_hat - theta0) / se_hat
  z_alpha <- qnorm(1 - alpha/2)
  reject <- abs(z_stat) > z_alpha

  ## Return all components
  list(
    theta_hat   = theta_hat,
    se_hat      = se_hat,
    reject      = reject,
    lambda_hat  = lambda_hat,
    beta_hat    = beta_hat,
    H_L = H_L, 
    H_U = H_U,
    Sigma_YY_hat = Sigma_YY,
    Sigma_ff_hat = Sigma_ff,
    Sigma_Yf_hat = Sigma_Yf,
    A_hat        = A_scalar,
    B_hat        = B_scalar,
    C_hat        = C_scalar
  )
}

#' Monte Carlo Estimation of Empirical and Closed-Form Power for OLS PPI / PPI++
#'
#' @description
#' Runs a full Monte Carlo simulation for estimating power for the
#' OLS-based naive, PPI, or PPI++ estimators. This function repeatedly
#' calls \code{ppi_ols_rep_one()} and aggregates:
#' \enumerate{
#'   \item Empirical rejection rate (empirical power),
#'   \item Monte Carlo SE,
#'   \item Estimated quadratic forms
#'     \eqn{A = a^\top H^{-1}\Sigma_{YY}H^{-1}a},
#'     \eqn{B = a^\top H^{-1}\Sigma_{ff}H^{-1}a},
#'     \eqn{C = a^\top H^{-1}\Sigma_{Yf}H^{-1}a},
#'   \item Closed-form variance and Wald power,
#'   \item Diagnostic summary table.
#' }
#'
#' @param R Number of Monte Carlo repetitions.
#' @param n Number of labeled observations.
#' @param N Number of unlabeled observations.
#' @param X_sampler_L Covariate generator.
#' @param f_generator Signal generator.
#' @param model_type Predictive model class for \eqn{\hat f}.
#' @param a Contrast vector of length \eqn{p}.
#' @param delta Effect size.
#' @param theta0 True theta.
#' @param ppi_type `"PPI"` or `"PPI++"`.
#' @param lambda_mode `"plugin"` or `"oracle"`.
#' @param lambda_external Whether to estimate lambda from an external labeled dataset.
#' @param n_external Size of external labeled dataset (if used).
#' @param PPIpp_crossfit Whether to use cross-fitting in PPI++.
#' @param alpha Wald test significance level.
#' @param seed Optional seed.
#' @param keep_reps Logical; whether or not output intermediate statistics. Default FALSE.
#'
#' @return A list with:
#'   \item{empirical_power}{Empirical rejection rate over R reps.}
#'   \item{mc_se}{Monte Carlo SE of power.}
#'   \item{avg_theta_hat}{mean(\eqn{\hat\theta}).}
#'   \item{avg_se_hat}{mean estimated SE.}
#'   \item{avg_lambda_hat}{mean(\eqn{\hat\lambda}).}
#'   \item{A_hat, B_hat, C_hat}{Monte Carlo means of quadratic forms.}
#'   \item{lambda_star_hat}{Closed-form minimizing λ.}
#'   \item{se_theory}{Closed-form SE.}
#'   \item{theoretical_power}{Closed-form Wald power.}
#'   \item{cov_var_table}{Tidy diagnostic table.}
#'   \item{reps}{(Optional) Raw replicate-level results, if \code{keep_reps = TRUE}.}
#'
#' @export
#' 
ppi_ols_empirical_power <- function(
  R,
  n, N,
  X_sampler_L,
  f_generator,
  model_type,
  a,
  delta,
  theta0,
  ppi_type = c("PPI", "PPI++"),
  lambda_mode = c("plugin", "oracle"),
  lambda_external = FALSE,
  n_external = ceiling(n/2),
  PPIpp_crossfit = TRUE,
  alpha = 0.05,
  seed = 1,
  keep_reps = FALSE
) {

  ppi_type    <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  set.seed(seed)
  a <- as.numeric(a)

  ## Run R Monte Carlo replicates
  reps <- replicate(
    R,
    ppi_ols_rep_one(
      n = n, N = N,
      X_sampler_L = X_sampler_L,
      f_generator = f_generator,
      model_type   = model_type,
      a            = a,
      delta        = delta,
      theta0       = theta0,
      ppi_type     = ppi_type,
      lambda_mode  = lambda_mode,
      lambda_external = lambda_external,
      n_external   = n_external,
      PPIpp_crossfit = PPIpp_crossfit,
      alpha        = alpha
    ),
    simplify = FALSE
  )

  ## Extract replicate-level results
  reject_vec <- vapply(reps, `[[`, logical(1), "reject")
  theta_vec  <- vapply(reps, `[[`, numeric(1), "theta_hat")
  se_vec     <- vapply(reps, `[[`, numeric(1), "se_hat")
  lambda_vec <- vapply(reps, `[[`, numeric(1), "lambda_hat")

  empirical_power <- mean(reject_vec)
  mc_se_power     <- sqrt(empirical_power * (1 - empirical_power) / R)

  A_vals <- vapply(reps, `[[`, numeric(1), "A_hat")
  B_vals <- vapply(reps, `[[`, numeric(1), "B_hat")
  C_vals <- vapply(reps, `[[`, numeric(1), "C_hat")

  z_alpha <- qnorm(1 - alpha/2)

  ## analytic contrast shift:  θ* - θ0 = δ‖a‖^2
  theta_shift <- delta * sum(a * a)

  se_theory_vec <- se_vec
  power_theory_vec <- numeric(R)

  for (r in seq_len(R)) {
    se_theory_r <- se_theory_vec[r]
    mu_r <- theta_shift / se_theory_r

    power_theory_vec[r] <-
      pnorm(-z_alpha + mu_r) + (1 - pnorm(z_alpha + mu_r))
  }

  ## Averages across R
  A_hat <- mean(A_vals)
  B_hat <- mean(B_vals)
  C_hat <- mean(C_vals)

  se_theory         <- mean(se_theory_vec)
  theoretical_power <- mean(power_theory_vec)

  ## lambda summary
  lambda_avg <- mean(lambda_vec)
  lambda_sd  <- sd(lambda_vec)

  ## Oracle lambda based on sample-averaged B,C
  r_ratio <- n / N
  lambda_star_hat <-
    if (ppi_type == "PPI") 1 else C_hat / ((1 + r_ratio) * B_hat)
  lambda_star_hat <- max(0, min(1, lambda_star_hat))

  ## Condition number
  H_condition_num <- mean(vapply(
    reps,
    function(x) {
      if (ppi_type == "PPI") {
        Hp <- x$H_U
      } else {
        Hp <- (1 - x$lambda_hat) * x$H_L + x$lambda_hat * x$H_U
      }
      kappa(Hp)
    },
    numeric(1)
  ))

  ## Covariance summary table 
  cov_var_table <- tibble::tibble(
    model_type       = model_type,
    ppi_type         = ppi_type,
    lambda_mode      = lambda_mode,
    lambda_external  = lambda_external,
    design           = ppi_type,
    n                = n,
    N                = N,
    A_hat            = A_hat,
    B_hat            = B_hat,
    C_hat            = C_hat,
    avg_lambda_hat   = lambda_avg,
    sd_lambda_hat    = lambda_sd,
    lambda_star_hat  = lambda_star_hat,
    H_condition_num  = H_condition_num,
    se_theory        = se_theory,
    power_theory     = theoretical_power
  )

  out <- list(
    empirical_power  = empirical_power,
    theoretical_power = theoretical_power,
    mc_se            = mc_se_power,
    avg_theta_hat    = mean(theta_vec),
    avg_se_hat       = mean(se_vec),
    avg_lambda_hat   = lambda_avg,
    lambda_star_hat = lambda_star_hat,
    se_theory         = se_theory,
    cov_var_table = cov_var_table
  )

  if (keep_reps) {
    out$reps <- reps
  }
  
  out
}
