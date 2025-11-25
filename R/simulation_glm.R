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
#' @param ppi_type Either "PPI" (\\eqn{\\lambda=1}) or "PPI++".
#' @param lambda_mode For PPI++: "plugin" (estimate \\eqn{\\lambda})
#'   or "oracle" (not implemented internally).
#' @param lambda_external Logical; if TRUE, estimate \\eqn{\\lambda}
#'   using an external labeled sample.
#' @param n_external Size of the external labeled set if
#'   \code{lambda_external = TRUE}.
#' @param PPIpp_crossfit Logical; if TRUE, use 2-fold cross-fitting
#'   for PPI++.
#' @param alpha Wald test significance level.
#' @param seed Optional random seed.
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

  ## True GLM signal & alternative perturbation
  eta_L <- f_generator(X_L_df)
  eta_U <- f_generator(X_U_df)

  # contrast-aligned shift: η_δ(x) = η(x) + δ x^T a
  eta_L_delta <- eta_L + delta * as.numeric(X_L %*% a)
  eta_U_delta <- eta_U + delta * as.numeric(X_U %*% a)

  # logistic mean μ = logistic(η)
  mu_L <- plogis(eta_L_delta)
  mu_U <- plogis(eta_U_delta)

  # generate Bernoulli response
  y_L <- rbinom(n, 1, mu_L)

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
        X_U_df
      )

      muhat_L[te] <- m$fhat_L
      muhat_U_folds[[fold]] <- m$fhat_U
    }

    muhat_U <- (muhat_U_folds[[1]] + muhat_U_folds[[2]]) / 2

  } else {
    m <- fit_predict_model_glm(model_type, X_L_df, y_L, X_U_df)
    muhat_L <- m$fhat_L
    muhat_U <- m$fhat_U
  }

  ## Canonical GLM matrices (Fisher information / Hessian)
  ## J := E[w(X) XXᵀ],  w = μ(1−μ) for logistic link
  w_L <- mu_L * (1 - mu_L)
  w_U <- mu_U * (1 - mu_U)

  # H_L and H_U correspond to J_L and J_U
  J_L <- crossprod(X_L * w_L, X_L) / n
  J_U <- crossprod(X_U * w_U, X_U) / N

  ## PPI λ-estimation
  if (ppi_type == "PPI") {
    lambda_hat <- 1
  } else {

    if (lambda_external) {

      ##  External Labeled Sample 
      sim_ext <- simulate_one_draw_contrast(
        n_labeled   = n_external,
        X_sampler_L = X_sampler_L,
        f_generator = f_generator,
        eps_sampler = NULL,
        a           = a,
        delta       = delta
      )

      X_ext_df <- sim_ext$X
      y_ext    <- sim_ext$y
      X_ext <- as.matrix(X_ext_df)

      ## Predict f on external data using model from main sample
      m_ext <- fit_predict_model_glm(
        model_type,
        X_L_df, y_L,
        X_ext_df
      )
      muhat_ext <- m_ext$fhat_U

      ## External moment matrices
      mu_ext_true <- plogis(f_generator(X_ext_df))
      resid_ext   <- y_ext - muhat_ext
      XY_ext      <- sweep(X_ext, 1, resid_ext, `*`)
      Sigma_YY_ext <- cov(XY_ext)

      ff_ext <- muhat_ext - mu_ext_true
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

      resid_L <- y_L - muhat_L
      XY_resid <- sweep(X_L, 1, resid_L, `*`)
      Sigma_YY <- cov(XY_resid)

      ff_L <- (muhat_L - mu_L)
      X_ff <- sweep(X_L, 1, ff_L, `*`)
      Sigma_ff <- cov(X_ff)

      Sigma_Yf <- cov(XY_resid, X_ff)

      Jinv_U <- solve(J_U)
      B <- t(a) %*% Jinv_U %*% Sigma_ff %*% Jinv_U %*% a
      C <- t(a) %*% Jinv_U %*% Sigma_Yf %*% Jinv_U %*% a

      r <- n / N
      lambda_hat <- as.numeric(C / ((1+r)*B))
      lambda_hat <- max(0, min(1, lambda_hat))

    } else {
      stop("Oracle lambda not implemented for GLM")
    }

  } ## end PPI++ branch

  ## Solve PPI / PPI++ estimating equation
  # score: U(β) = X (Y - μ_β)
  # PPI++ score: U_λ as in eq (31) of write-up

  # We use 1-step Newton update starting at β=0.
  beta_hat <- rep(0, p)

  for (iter in 1:3) {
    eta_hat_L <- as.numeric(X_L %*% beta_hat)
    mu_hat_L  <- plogis(eta_hat_L)

    # Score U_λ(β)
    U_L <- crossprod(X_L, (y_L - mu_hat_L)) / n
    U_U <- crossprod(X_U, (muhat_U - plogis(X_U %*% beta_hat))) / N
    U_f <- crossprod(X_L, (muhat_L - plogis(X_L %*% beta_hat))) / n

    if (ppi_type == "PPI") {
      U_lambda <- U_L + U_U
      Jeff <- J_U
    } else {
      U_lambda <- U_L + lambda_hat * (U_U - U_f)
      Jeff <- (1 - lambda_hat) * J_L + lambda_hat * J_U
    }

    beta_hat <- beta_hat + solve(Jeff, U_lambda)
  }

  ## Compute variance components at β_hat
  eta_hat <- as.numeric(X_L %*% beta_hat)
  mu_hat  <- plogis(eta_hat)

  resid_L <- y_L - mu_hat
  XY_resid <- sweep(X_L, 1, resid_L, `*`)
  Sigma_YY <- cov(XY_resid)

  ff_L <- muhat_L - mu_hat
  X_ff <- sweep(X_L, 1, ff_L, `*`)
  Sigma_ff <- cov(X_ff)

  Sigma_Yf <- cov(XY_resid, X_ff)

  ## Sandwich variance (eq 37)
  J_p <- if (ppi_type == "PPI") J_U else (1-lambda_hat)*J_L + lambda_hat*J_U
  Jinv <- solve(J_p)

  Mid <-
    Sigma_YY / n +
    lambda_hat^2 * (Sigma_ff / N + Sigma_ff / n) -
    2 * lambda_hat * Sigma_Yf / n

  Vhat <- t(a) %*% Jinv %*% Mid %*% Jinv %*% a
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
    J_L         = J_L,
    J_U         = J_U,
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
#' @param keep_reps Logical; if \code{TRUE}, return the full list of replicate
#'   outputs from \code{ppi_glm_rep_one()}.
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
  keep_reps = FALSE
) {

  ppi_type    <- match.arg(ppi_type)
  lambda_mode <- match.arg(lambda_mode)
  set.seed(seed)
  a <- as.numeric(a)

  ## Run R Monte Carlo replicates
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

  A_vals <- numeric(R)
  B_vals <- numeric(R)
  C_vals <- numeric(R)

  se_theory_vec    <- numeric(R)
  power_theory_vec <- numeric(R)

  z_alpha <- qnorm(1 - alpha/2)

  ## analytic contrast shift:  θ* - θ0 = δ‖a‖²
  theta_shift <- delta * sum(a * a)

  for (r in seq_len(R)) {
    x <- reps[[r]]

    lambda_r <- x$lambda_hat
    J_L_r    <- x$J_L
    J_U_r    <- x$J_U

    SYY_r <- x$Sigma_YY_hat
    Sff_r <- x$Sigma_ff_hat
    SYf_r <- x$Sigma_Yf_hat

    ## Fisher information for this replicate
    if (ppi_type == "PPI") {
      J_p_r <- J_U_r
    } else {
      J_p_r <- (1 - lambda_r) * J_L_r + lambda_r * J_U_r
    }
    Jinv_r <- solve(J_p_r)

    ## Sandwich components (GLM variance: eq 37)
    A_r <- as.numeric(t(a) %*% Jinv_r %*% SYY_r %*% Jinv_r %*% a)
    B_r <- as.numeric(t(a) %*% Jinv_r %*% Sff_r %*% Jinv_r %*% a)
    C_r <- as.numeric(t(a) %*% Jinv_r %*% SYf_r %*% Jinv_r %*% a)

    A_vals[r] <- A_r
    B_vals[r] <- B_r
    C_vals[r] <- C_r

    ## Variance for replicate r
    V_r <- A_r / n +
      lambda_r^2 * (B_r / N + B_r / n) -
      2 * lambda_r * C_r / n

    se_theory_r <- sqrt(V_r)
    se_theory_vec[r] <- se_theory_r

    ## Noncentral Z mean
    mu_r <- theta_shift / se_theory_r

    ## Theoretical power for replicate r
    p_r <- pnorm(-z_alpha + mu_r) + (1 - pnorm(z_alpha + mu_r))
    power_theory_vec[r] <- p_r
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

  ## Oracle lambda analytic (GLM version identical structure)
  r_ratio <- n / N
  lambda_star_hat <-
    if (ppi_type == "PPI") 1 else C_hat / ((1 + r_ratio) * B_hat)
  lambda_star_hat <- max(0, min(1, lambda_star_hat))

  ## Condition number
  H_condition_num <- mean(vapply(
    reps,
    function(x) {
      if (ppi_type == "PPI") {
        Jp <- x$J_U
      } else {
        Jp <- (1 - x$lambda_hat) * x$J_L + x$lambda_hat * x$J_U
      }
      kappa(Jp)
    },
    numeric(1)
  ))

  ## Covariance summary table 
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
    avg_lambda_hat   = lambda_avg,
    sd_lambda_hat    = lambda_sd,
    lambda_star_hat  = lambda_star_hat,
    H_condition_num  = H_condition_num,
    se_theory        = se_theory,
    power_theory     = theoretical_power
  )

  out <- list(
    empirical_power   = empirical_power,
    theoretical_power = theoretical_power,
    mc_se             = mc_se_power,
    avg_theta_hat     = mean(theta_vec),
    avg_se_hat        = mean(se_vec),
    avg_lambda_hat    = lambda_avg,
    lambda_star_hat   = lambda_star_hat,
    se_theory         = se_theory,
    cov_var_table     = cov_var_table
  )

  if (keep_reps)
    out$reps <- reps

  out
}