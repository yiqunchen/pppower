#' One Monte Carlo Replicate for GLM–Based PPI / PPI++ Inference
#'
#' @description
#' Performs one Monte Carlo replication of the Prediction-Powered Inference
#' (PPI) or PPI++ estimator for a scalar linear contrast
#' \eqn{\theta = a^\top \beta^\star} in a **generalized linear model (GLM)**
#' with canonical link (e.g., logistic regression).
#'
#' This function implements:
#' \itemize{
#'   \item the canonical score formulation \eqn{U(\beta)} from Section 4,
#'   \item the PPI and PPI++ estimating equations,
#'   \item the contrast-aligned alternative 
#'         \eqn{\eta_\delta(x) = \eta(x) + \delta (x^\top a)},
#'   \item cross-fitted prediction models for \eqn{\mu_f(x)},
#'   \item sandwich variance formulas for GLMs
#'         (equations (35)–(38) of the write-up),
#'   \item λ-estimation for PPI++ using plug-in covariance estimators.
#' }
#'
#' The returned object contains all components needed to compute **empirical
#' power**, **theoretical (closed-form) power**, and **λ-optimality analyses**
#' for GLM PPI/PPI++.
#'
#'
#' @keywords internal
#' 
#' @param n Number of labeled observations.
#' @param N Number of unlabeled observations.
#' @param X_sampler_L Function generating covariate matrices (\code{n x p}).
#' @param f_generator Function generating the regression signal \eqn{f(X)}.
#' @param eps_sampler Noise generator, default NULL. Ignored for Bernoulli outcomes.
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
#' @param use_true_J Logical; if TRUE, use true Fisher J instead of estimated.
#' 
#' @return A list with:
#'   \item{theta_hat}{Estimated contrast \eqn{a^\top \hat\beta}.}
#'   \item{se_hat}{Estimated standard error of \eqn{a^\top \hat\beta}.}
#'   \item{reject}{Logical indicator for Wald test rejection of
#'     \eqn{H_0 : \theta = \theta_0}.}
#'   \item{lambda_hat}{Estimated PPI++ shrinkage parameter \eqn{\hat\lambda}.}
#'   \item{beta_hat}{Estimated GLM coefficient vector \eqn{\hat\beta}.}
#'   \item{J_L}{Estimated labeled Fisher information matrix
#'     \eqn{J_L = \frac{1}{n} \sum_{i=1}^n w(X_i)\,X_i X_i^\top}.}
#'   \item{J_U}{Estimated unlabeled Fisher information matrix
#'     \eqn{J_U = \frac{1}{N} \sum_{j=1}^N w(\tilde X_j)\,\tilde X_j \tilde X_j^\top}.}
#'   \item{Sigma_YY_hat}{Sample analogue of the score–residual covariance
#'     \eqn{\Sigma_{YY} = \mathrm{Var}\big(X\{Y - \mu_{\beta^\star}(X)\}\big)}.}
#'   \item{Sigma_ff_hat}{Sample analogue of the prediction–error covariance
#'     \eqn{\Sigma_{ff} = \mathrm{Var}\big(X\{\mu_f(X) - \mu_{\beta^\star}(X)\}\big)}.}
#'   \item{Sigma_Yf_hat}{Sample analogue of the cross–covariance
#'     \eqn{\Sigma_{Yf} = \mathrm{Cov}\big(
#'       X\{Y - \mu_{\beta^\star}(X)\},\,
#'       X\{\mu_f(X) - \mu_{\beta^\star}(X)\}
#'     \big)}.}
#' 
#' @importFrom stats plogis rbinom
#' 
ppi_glm_rep_one <- function(
  n, N,
  X_sampler_L,
  f_generator,                 # returns linear predictor η(x)
  eps_sampler = NULL,          # unused for Bernoulli
  model_type = c("glm_correct", "glm_mis", "glm_wrong", "rf", "oracle"),
  a,
  delta = 0,
  ppi_type = c("PPI", "PPI++"),
  theta0 = 0,
  lambda_mode = c("plugin", "oracle"),
  lambda_external = FALSE,
  n_external = ceiling(n/2),
  PPIpp_crossfit = TRUE,
  alpha = 0.05,
  seed = NULL,
  use_true_J = FALSE
) {

  if (!is.null(seed)) set.seed(seed)

  model_type <- match.arg(model_type)
  ppi_type   <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  a <- as.numeric(a)

  ## Generate X
  X_L_df <- X_sampler_L(n)
  X_U_df <- X_sampler_L(N)

  X_L <- as.matrix(X_L_df)
  X_U <- as.matrix(X_U_df)
  p   <- ncol(X_L)
  stopifnot(length(a) == p)

  # baseline linear predictor for β_true
  eta_L_base <- f_generator(X_L_df)
  eta_U_base <- f_generator(X_U_df)

  ## Apply delta shift
  eta_L_star <- eta_L_base + delta * as.numeric(X_L %*% a)
  eta_U_star <- eta_U_base + delta * as.numeric(X_U %*% a)

  mu_L_star <- plogis(eta_L_star)   # true mean for labeled
  mu_U_star <- plogis(eta_U_star)   # true mean for unlabeled

  # generate Bernoulli response
  y_L <- rbinom(n, 1, mu_L_star)

  ## Fit predictive model: μ_f(x)
  if (ppi_type == "PPI++" && PPIpp_crossfit) {
    fold_id <- sample(rep(1:2, length.out = n))
    muhat_L <- numeric(n)
    muhat_U_folds <- vector("list", 2)

    for (fold in 1:2) {
      tr <- which(fold_id != fold)
      te <- which(fold_id == fold)

      m <- fit_predict_model_glm(
        model_type,
        X_L_df[tr, , drop=FALSE],
        y_L[tr],
        X_U_df,
        f_generator = f_generator
      )

      muhat_L[te] <- m$fhat_L
      muhat_U_folds[[fold]] <- m$fhat_U
    }

    muhat_U <- (muhat_U_folds[[1]] + muhat_U_folds[[2]]) / 2

  } else {
    m <- fit_predict_model_glm(model_type, X_L_df, y_L, X_U_df, f_generator = f_generator)
    muhat_L <- m$fhat_L
    muhat_U <- m$fhat_U
  }

  ## Canonical GLM matrices (Fisher information / Hessian)
  ## J := E[w(X) XXᵀ],  w = μ(1−μ) for logistic link
  w_L <- mu_L_star * (1 - mu_L_star)
  w_U <- mu_U_star * (1 - mu_U_star)

  # H_L and H_U correspond to J_L and J_U
  J_L <- crossprod(X_L * w_L, X_L) / n
  J_U <- crossprod(X_U * w_U, X_U) / N

  ## PPI λ-estimation
  if (ppi_type == "PPI") {
    lambda_hat <- 1 # NOT USED
  } else {

    if (lambda_external) {

      ##  External Labeled Sample 
      sim_ext <- simulate_one_draw_contrast_glm(
        n_labeled   = n_external,
        X_sampler_L = X_sampler_L,
        f_generator = f_generator,
        a           = a,
        delta       = delta
      )

      X_ext_df <- sim_ext$X
      y_ext    <- sim_ext$y
      X_ext <- as.matrix(X_ext_df)

      ## Predict f on external data using model from main sample
      m_ext <- fit_predict_model_glm(model_type, X_L_df, y_L, X_ext_df, f_generator = f_generator)
      muhat_ext <- m_ext$fhat_U

      eta_ext <- f_generator(X_ext_df) + delta * as.numeric(X_ext %*% a)
      mu_ext <- plogis(eta_ext)

      resid_ext <- y_ext - mu_ext
      XY_ext <- sweep(X_ext, 1, resid_ext, `*`)
      Sigma_YY_ext <- cov(XY_ext)

      ## prediction error: (μ_f - μ_star)
      ff_ext <- muhat_ext - mu_ext
      X_ff_ext <- sweep(X_ext, 1, ff_ext, `*`)
      Sigma_ff_ext <- cov(X_ff_ext)

      Sigma_Yf_ext <- cov(XY_ext, X_ff_ext)

      ## External-based B, C
      Jinv_U <- solve(J_U)

      B_ext <- t(a) %*% Jinv_U %*% Sigma_ff_ext %*% Jinv_U %*% a
      C_ext <- t(a) %*% Jinv_U %*% Sigma_Yf_ext %*% Jinv_U %*% a

      r <- n / N
      lambda_hat <- as.numeric(C_ext / ((1+r) * B_ext)) 
      lambda_hat <- max(0, min(1, lambda_hat))

    } else if (lambda_mode == "plugin") {

      resid_L <- y_L - mu_L_star
      XY_resid <- sweep(X_L, 1, resid_L, `*`)
      Sigma_YY <- cov(XY_resid)

      ff_L <- muhat_L - mu_L_star
      X_ff <- sweep(X_L, 1, ff_L, `*`)
      Sigma_ff <- cov(X_ff)

      Sigma_Yf <- cov(XY_resid, X_ff)

      Jinv_U <- solve(J_U)

      B <- t(a) %*% Jinv_U %*% Sigma_ff %*% Jinv_U %*% a
      C <- t(a) %*% Jinv_U %*% Sigma_Yf %*% Jinv_U %*% a

      r <- n / N
      lambda_hat <- as.numeric(C / ((1+r)*B))  ## Eq (29)
      lambda_hat <- max(0, min(1, lambda_hat))

    } else {
      stop("Oracle lambda not implemented for GLM")
    }

  } ## end PPI++ branch

  ## Solve one-step PPI / PPI++ estimating equation by Newton iterations
  # score: U(β) = X (Y - μ_β)

  maxit <- 25L
  tol   <- 1e-8
  ridge <- 1e-4   # larger ridge for numerical stability

  glm_init <- glm.fit(X_L, y_L, family = binomial())
  beta     <- glm_init$coefficients

  for (iter in seq_len(maxit)) {

    eta_L_hat <- as.numeric(X_L %*% beta)
    eta_U_hat <- as.numeric(X_U %*% beta)
    mu_L_hat  <- plogis(eta_L_hat)
    mu_U_hat  <- plogis(eta_U_hat)

    ## clip logits to avoid mu = 0 or 1
    eps <- 1e-8
    mu_L_hat <- pmin(pmax(mu_L_hat, eps), 1 - eps)
    mu_U_hat <- pmin(pmax(mu_U_hat, eps), 1 - eps)

    ## scores
    U_L <- crossprod(X_L, (y_L      - mu_L_hat)) / n
    U_U <- crossprod(X_U, (muhat_U  - mu_U_hat)) / N
    U_f <- crossprod(X_L, (muhat_L  - mu_L_hat)) / n

    ## Hessians
    w_L_hat <- mu_L_hat * (1 - mu_L_hat)
    w_U_hat <- mu_U_hat * (1 - mu_U_hat)

    J_L_hat <- crossprod(X_L * w_L_hat, X_L) / n
    J_U_hat <- crossprod(X_U * w_U_hat, X_U) / N

    if (ppi_type == "PPI") {
      score <- U_L + U_U
      J_eff <- J_L_hat + J_U_hat
    } else {
      score <- U_L + lambda_hat*(U_U - U_f)
      J_eff <- (1 - lambda_hat)*J_L_hat + lambda_hat*J_U_hat
    }

    ## regularize and guard against non-finite entries
    J_eff_reg <- J_eff + diag(ridge, p)

    if (!all(is.finite(J_eff_reg)) || !all(is.finite(score))) {
      ## fall back to current beta and stop iterating
      warning("Non-finite J_eff or score encountered; using current beta")
      break
    }

    step <- tryCatch(
      solve(J_eff_reg, score),
      error = function(e) rep(NA_real_, p)
    )

    ## if step failed or has NA, break and keep current beta
    if (!all(is.finite(step))) {
      warning("Newton step non-finite; using current beta")
      break
    }

    beta_new <- beta + step

    ## SAFE convergence check
    if (max(abs(step)) < tol) {
      beta <- beta_new
      break
    }

    beta <- beta_new
  }

  beta_hat <- as.numeric(beta)

  ## Compute variance components at β_hat
  eta_L_hat <- as.numeric(X_L %*% beta_hat)
  eta_U_hat <- as.numeric(X_U %*% beta_hat)
  mu_L_hat  <- plogis(eta_L_hat)
  mu_U_hat  <- plogis(eta_U_hat)

  w_L_hat <- mu_L_hat * (1 - mu_L_hat)
  w_U_hat <- mu_U_hat * (1 - mu_U_hat)

  J_L_hat <- crossprod(X_L * w_L_hat, X_L) / n
  J_U_hat <- crossprod(X_U * w_U_hat, X_U) / N

  resid_L <- y_L - mu_L_star         ### FIXED
  XY_resid <- sweep(X_L, 1, resid_L, `*`)
  Sigma_YY <- cov(XY_resid)

  ff_L <- muhat_L - mu_L_star        ### FIXED
  X_ff <- sweep(X_L, 1, ff_L, `*`)
  Sigma_ff <- cov(X_ff)

  Sigma_Yf <- cov(XY_resid, X_ff)

  ## Sandwich variance (eq 25)
  if (use_true_J) {
      ## TRUE Fisher: evaluated at μ_star, not μ_hat
      J_L_eff <- J_L   # from mu_L_star
      J_U_eff <- J_U   # from mu_U_star
  } else {
      ## ESTIMATED Fisher: evaluated at mu_hat(beta_hat)
      J_L_eff <- J_L_hat
      J_U_eff <- J_U_hat
  }


if (ppi_type == "PPI") {
    ## PPI Fisher + variance 
    J_p_r <- J_L_eff + J_U_eff
    lambda_eff <- 1

    Mid <- 
      Sigma_YY / n +
      Sigma_ff / N

} else {
    ## PPI++ 
    J_p_r <- (1 - lambda_hat) * J_L_eff + lambda_hat * J_U_eff
    lambda_eff <- lambda_hat

    Mid <-
      Sigma_YY / n +
      lambda_eff^2 * (Sigma_ff / N + Sigma_ff / n) -
      2 * lambda_eff * Sigma_Yf / n
}

  ridge <- 1e-8
  Jinv_r <- solve(J_p_r + diag(ridge, ncol(J_p_r)))

  Vhat <- t(a) %*% Jinv_r %*% Mid %*% Jinv_r %*% a
  se_hat <- sqrt(as.numeric(Vhat))

  ## Wald test
  theta_hat <- sum(a * beta_hat)
  z_stat <- (theta_hat - theta0) / se_hat
  reject <- abs(z_stat) > qnorm(1 - alpha/2)

  ## Return
  list(
      theta_hat   = theta_hat,
      se_hat      = se_hat,
      reject      = reject,
      lambda_hat  = lambda_hat,
      beta_hat    = beta_hat,
      # TRUE Fisher info (from true μ*)
      J_L_true = J_L,
      J_U_true = J_U,
      # ESTIMATED Fisher info (from μ_hat(β_hat))
      J_L_hat  = J_L_hat,
      J_U_hat  = J_U_hat,
      # EFFECTIVE Fisher info for variance calculation
      J_L_eff = J_L_eff,
      J_U_eff = J_U_eff,
      J_p_r = J_p_r,
      Sigma_YY_hat = Sigma_YY,
      Sigma_ff_hat = Sigma_ff,
      Sigma_Yf_hat = Sigma_Yf
  )
}

#' Monte Carlo Estimation of Empirical and Closed-Form Power for GLM PPI / PPI++
#'
#' @description
#' Runs a full Monte Carlo study for inference on a scalar linear contrast
#' \eqn{\theta = a^\top\beta^\star} under a generalized linear model (GLM),
#' using either PPI (\eqn{\lambda=1}) or PPI++ (\eqn{\lambda\in[0,1]}).
#' 
#' This function repeatedly calls \code{ppi_glm_rep_one()} to generate
#' replicate estimates \eqn{\hat\theta}, \eqn{\widehat{\mathrm{se}}}, and
#' \eqn{\hat\lambda}, and aggregates:
#' \enumerate{
#'   \item Empirical rejection probability (empirical power),
#'   \item Monte Carlo standard error of power,
#'   \item Estimated quadratic forms
#'     \eqn{
#'     A = a^\top J^{-1}\Sigma_{YY}J^{-1}a,\quad
#'     B = a^\top J^{-1}\Sigma_{ff}J^{-1}a,\quad
#'     C = a^\top J^{-1}\Sigma_{Yf}J^{-1}a,
#'     }
#'     where \eqn{J} is the Fisher information matrix,
#'   \item Closed-form asymptotic variance and Wald power for the GLM,
#'   \item A tidy diagnostic summary table of all variance components.
#' }
#'
#' The GLM formulation follows the canonical-score expansions described in
#' Section~4 of the accompanying manuscript.  For logistic regression,
#' \eqn{J = \mathbb{E}[\,\mu(1-\mu)XX^\top\,]}.
#'
#'
#' @param R Number of Monte Carlo repetitions.
#' @param n Number of labeled observations.
#' @param N Number of unlabeled observations.
#' @param X_sampler_L Function that generates labeled covariates \eqn{X_L}.
#' @param f_generator Function generating the true linear predictor
#'   \eqn{\eta(X)}, before applying the logistic (or general GLM) link.
#' @param model_type Predictive model class for the nuisance predictor
#'   \eqn{\mu_f}, e.g.\ \code{"glm_correct"}, \code{"glm_wrong"}, \code{"rf"}.
#' @param a Contrast vector \eqn{a\in\mathbb{R}^p}.
#' @param delta Effect size injected via \eqn{\eta_\delta(x)
#'   = \eta(x) + \delta\,x^\top a}.
#' @param theta0 Null value of the contrast \eqn{\theta}.
#' @param ppi_type `"PPI"` or `"PPI++"`.
#' @param lambda_mode `"plugin"` or `"oracle"`.  Only `"plugin"` is
#'   implemented internally.
#' @param lambda_external Logical; whether to estimate \eqn{\lambda} from an
#'   external labeled sample.
#' @param n_external Size of the external labeled sample (if used).
#' @param PPIpp_crossfit Logical; whether to use 2-fold cross-fitting for PPI++.
#' @param alpha Wald test significance level.
#' @param seed Optional random seed.
#' @param use_true_J Logical; if TRUE, use true Fisher J instead of estimated.
#' @param keep_reps Logical; if \code{TRUE}, return the full list of replicate
#'   outputs from \code{ppi_glm_rep_one()}.
#' @param theta_shift_mode Character string specifying how the effect shift is applied.
#'   Must be one of:
#'   \describe{
#'     \item{`"population"`}{Shift is applied to the population-level
#'       contrast \eqn{c^\top \beta}.}
#'     \item{`"empirical_global"`}{Shift is applied once to a globally estimated
#'       contrast, shared by all Monte Carlo replicates.}
#'     \item{`"empirical_rep"`}{Shift is applied separately **within each** Monte Carlo
#'       replicate based on its own empirical estimate.}
#'   }
#'
#' @return A list with:
#'   \item{empirical_power}{Empirical rejection probability over \eqn{R} reps.}
#'   \item{theoretical_power}{Closed-form Wald power (GLM).}
#'   \item{mc_se}{Monte Carlo standard error of empirical power.}
#'   \item{avg_theta_hat}{Mean of \eqn{\hat\theta}.}
#'   \item{avg_se_hat}{Mean estimated standard error.}
#'   \item{avg_lambda_hat}{Mean of the fitted \eqn{\hat\lambda}.}
#'   \item{lambda_star_hat}{Closed-form minimizer \eqn{\lambda^\star(a)} of
#'       the GLM variance (oracle value).}
#'   \item{se_theory}{Closed-form asymptotic standard error for the GLM contrast.}
#'   \item{cov_var_table}{A tidy \code{tibble} summarizing
#'     \eqn{A}, \eqn{B}, \eqn{C}, \eqn{\lambda}, condition number diagnostics,
#'     and theoretical variance/power.}
#'   \item{reps}{(Optional) Raw replicate-level results, if
#'     \code{keep_reps = TRUE}.}
#'
#' @details
#' The GLM variance decomposition uses:
#' \deqn{
#' \mathrm{Var}\big(a^\top\hat\beta_\lambda\big)
#'   \approx
#'   a^\top J^{-1}\!\left(
#'     \frac{\Sigma_{YY}}{n}
#'     + \lambda^2\Big(\frac{\Sigma_{ff}}{N}+\frac{\Sigma_{ff}}{n}\Big)
#'     - \frac{2\lambda}{n}\Sigma_{Yf}
#'   \right)J^{-1}a,
#' }
#' where \eqn{J} is the Fisher information matrix.
#'
#' The effect size shift used for power calculations is
#' \eqn{\theta_\delta - \theta_0 = \delta\,\|a\|_2^2}.
#'
#' @export
#' 
ppi_glm_empirical_power <- function(
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
  use_true_J = FALSE,
  keep_reps = FALSE,
  theta_shift_mode = c("population", "empirical_global", "empirical_rep")
) {

  ppi_type    <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  theta_shift_mode <- match.arg(theta_shift_mode)
  set.seed(seed)
  a <- as.numeric(a)
  z_alpha <- qnorm(1 - alpha/2)

  ## population-level theoretical shift (delta * ||a||^2)
  theta_shift_pop <- delta * sum(a * a)

  ## Monte Carlo Replicates
  reps <- replicate(
    R,
    ppi_glm_rep_one(
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
      alpha        = alpha,
      use_true_J   = use_true_J
    ),
    simplify = FALSE
  )

  ## extract result vectors
  reject_vec <- vapply(reps, `[[`, logical(1), "reject")
  theta_vec  <- vapply(reps, `[[`, numeric(1), "theta_hat")
  se_vec     <- vapply(reps, `[[`, numeric(1), "se_hat")
  lambda_vec <- vapply(reps, `[[`, numeric(1), "lambda_hat")

  empirical_power <- mean(reject_vec)
  mc_se_power     <- sqrt(empirical_power * (1 - empirical_power) / R)

  ## θ-shift corrections:
  if (theta_shift_mode == "population") {
    theta_shift_global <- theta_shift_pop
  } 
  if (theta_shift_mode == "empirical_global") {
    theta_shift_global <- mean(theta_vec) - theta0
  }

  ## Compute theoretical power
  A_vals <- numeric(R)
  B_vals <- numeric(R)
  C_vals <- numeric(R)
  se_theory_vec    <- numeric(R)
  power_theory_vec <- numeric(R)

  for (r in seq_len(R)) {
    x <- reps[[r]]
    lambda_r <- x$lambda_hat

    ## effective Fisher information
    J_p_r <- x$J_p_r
    SYY_r <- x$Sigma_YY_hat
    Sff_r <- x$Sigma_ff_hat
    SYf_r <- x$Sigma_Yf_hat

    ridge <- 1e-8
    Jinv_r <- solve(J_p_r + diag(ridge, ncol(J_p_r)))

    ## Quadratic forms
    A_r <- as.numeric(t(a) %*% Jinv_r %*% SYY_r %*% Jinv_r %*% a)
    B_r <- as.numeric(t(a) %*% Jinv_r %*% Sff_r %*% Jinv_r %*% a)
    C_r <- as.numeric(t(a) %*% Jinv_r %*% SYf_r %*% Jinv_r %*% a)

    A_vals[r] <- A_r
    B_vals[r] <- B_r
    C_vals[r] <- C_r

    ## lambda_eff
    lambda_eff <- if (ppi_type == "PPI") 1 else lambda_r

    ## variance formula
    if (ppi_type == "PPI") {
      V_r <- A_r / n + B_r / N
    } else {
      V_r <- A_r / n +
        lambda_eff^2 * (B_r / N + B_r / n) -
        2 * lambda_eff * C_r / n
    }

    se_theory_r <- sqrt(V_r)
    se_theory_vec[r] <- se_theory_r

    ## μ_r: the noncentrality parameter

    if (theta_shift_mode == "population") {
      theta_shift_r <- theta_shift_pop
    }
    if (theta_shift_mode == "empirical_global") {
      theta_shift_r <- theta_shift_global
    }
    if (theta_shift_mode == "empirical_rep") {
      theta_shift_r <- theta_vec[r] - theta0
    }

    mu_r <- theta_shift_r / se_theory_r

    ## Wald Power
    power_theory_vec[r] <-
      pnorm(-z_alpha + mu_r) + (1 - pnorm(z_alpha + mu_r))
  }

  ## summary
  theoretical_power <- mean(power_theory_vec)
  A_hat <- mean(A_vals)
  B_hat <- mean(B_vals)
  C_hat <- mean(C_vals)

  r_ratio <- n/N
  lambda_star_hat <- if (ppi_type == "PPI") 1 else C_hat / ((1+r_ratio)*B_hat)
  lambda_star_hat <- max(0, min(1, lambda_star_hat))

  cov_var_table <- tibble::tibble(
    model_type       = model_type,
    ppi_type         = ppi_type,
    lambda_mode      = lambda_mode,
    lambda_external  = lambda_external,
    n                = n,
    N                = N,
    A_hat            = A_hat,
    B_hat            = B_hat,
    C_hat            = C_hat,
    avg_lambda_hat   = mean(lambda_vec),
    sd_lambda_hat    = sd(lambda_vec),
    lambda_star_hat  = lambda_star_hat,
    se_theory        = mean(se_theory_vec),
    power_theory     = theoretical_power,
    theta_shift_mode = theta_shift_mode
  )

  out <- list(
    empirical_power   = empirical_power,
    theoretical_power = theoretical_power,
    mc_se             = mc_se_power,
    avg_theta_hat     = mean(theta_vec),
    avg_se_hat        = mean(se_vec),
    avg_lambda_hat    = mean(lambda_vec),
    lambda_star_hat   = lambda_star_hat,
    se_theory         = mean(se_theory_vec),
    cov_var_table     = cov_var_table
  )

  if (keep_reps) out$reps <- reps

  out
}