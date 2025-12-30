# One Monte Carlo replicate of given configurations
# Simulate labeled + unlabeled data, fits working model, builds matrices, compute the test statistics
ppi_ols_rep_one <- function(
  n, N,
  X_sampler_L,
  f_generator,
  eps_sampler = function(n) rnorm(n, 0, 1),
  model_type = c("glm_correct", "glm_mis", "glm_wrong", "rf"),
  a,                  # contrast vector length p
  delta = 0,          # effect size along that contrast
  alpha = 0.05,
  seed = NULL
) {
  model_type <- match.arg(model_type)
  if (!is.null(seed)) set.seed(seed)

  # simulate labeled + unlabeled X
  X_L <- X_sampler_L(n)        # n x p
  X_U <- X_sampler_L(N)        # N x p (same sampler)
  p   <- ncol(X_L)

  stopifnot(length(a) == p)

  # True signal
  f_L_true <- f_generator(X_L)
  f_U_true <- f_generator(X_U)

  # Inject effect along direction 'a' (same as mean shift):
  # we interpret delta as a shift in the contrast a^T beta*, i.e. bump the component
  # of f(X) aligned with 'a'. Simplest way: add delta * (X a) to f_true.
  f_L_shifted <- f_L_true + delta * as.numeric(X_L %*% a)
  f_U_shifted <- f_U_true + delta * as.numeric(X_U %*% a)

  # labeled outcomes Y = shifted f(X) + noise
  y_L <- f_L_shifted + eps_sampler(n)

  # fit working model on labeled data to get \hat f
  fit_obj <- fit_predict_model(
    model_type = model_type,
    X_L = X_L,
    y_L = y_L,
    X_U = X_U
  )
  # predictions
  fhat_L <- as.numeric(.predict_any(fit_obj$fit, X_L))
  fhat_U <- as.numeric(fit_obj$fhat_U)

  # build \hat H and \hat G
  # Hhat = (1/N) sum_j X_U[j,] X_U[j,]^T  (p x p)
  Hhat <- crossprod(X_U, X_U) / N

  # Ghat = (1/N) sum_j X_U[j,] fhat_U[j]  + (1/n) sum_i X_L[i,] * (y_L[i] - fhat_L[i])
  term_unl <- crossprod(X_U, fhat_U) / N         # p-vector
  term_lab <- crossprod(X_L, (y_L - fhat_L)) / n # p-vector
  Ghat <- as.numeric(term_unl + term_lab)

  # beta_hat_PPI = solve(Hhat, Ghat)
  beta_hat <- solve(Hhat, Ghat)

  # target contrast
  theta_hat <- sum(a * beta_hat)   # a^T beta_hat

  # estimate plug-in variance for Wald test
  # We need Sigma_ff and Sigma_YY (population-style) using sample analogues:
  #
  # Sigma_YY = Var( X * (Y - fhat_L) )  [p x p]
  resid_L <- y_L - fhat_L
  XY_resid <- sweep(X_L, 1, resid_L, `*`)  # n x p: row i = X_i * (Y_i - fhat_i)
  Sigma_YY_hat <- cov(XY_resid)            # sample cov -> p x p

  # Sigma_ff = Var( X * (fhat - X beta_hat) )  [p x p]
  lin_fit_L <- as.numeric(X_L %*% beta_hat)
  ff_L <- fhat_L - lin_fit_L
  X_ff <- sweep(X_L, 1, ff_L, `*`)         # n x p approx;
  Sigma_ff_hat <- cov(X_ff)                # p x p

  # sandwich variance:
  # Var(a^T beta_hat) \approx a^T H^{-1} [ Sigma_ff/N + Sigma_YY/n ] H^{-1} a
  Hinv <- solve(Hhat)
  Mid  <- (Sigma_ff_hat / N) + (Sigma_YY_hat / n)
  Vhat <- t(a) %*% Hinv %*% Mid %*% Hinv %*% a
  Vhat <- as.numeric(Vhat)
  se_hat <- sqrt(Vhat)

  # Wald test H0: a^T beta = 0
  z_stat <- theta_hat / se_hat
  z_alpha <- qnorm(1 - alpha/2)
  reject <- abs(z_stat) > z_alpha

  list(
    theta_hat   = theta_hat,
    se_hat      = se_hat,
    z_stat      = z_stat,
    reject      = reject,
    beta_hat    = beta_hat,
    Hhat        = Hhat,
    Sigma_YY_hat = Sigma_YY_hat,
    Sigma_ff_hat = Sigma_ff_hat
  )
}

ppi_ols_empirical_power <- function(
  R,
  n, N,
  X_sampler_L,
  f_generator,
  model_type,
  a,
  delta,
  alpha = 0.05,
  seed = 1
) {
  set.seed(seed)

  reps <- replicate(R, {
    out <- ppi_ols_rep_one(
      n = n, N = N,
      X_sampler_L = X_sampler_L,
      f_generator = f_generator,
      model_type = model_type,
      a = a,
      delta = delta,
      alpha = alpha
    )
    c(
      reject   = out$reject,
      theta    = out$theta_hat,
      se       = out$se_hat
    )
  }, simplify = TRUE)

  reject_vec <- reps["reject", ]
  theta_vec  <- reps["theta", ]
  se_vec     <- reps["se", ]

  empirical_power <- mean(reject_vec)
  mc_se_power     <- sqrt(empirical_power * (1 - empirical_power) / R)

  list(
    empirical_power = empirical_power,
    mc_se           = mc_se_power,
    avg_theta_hat   = mean(theta_vec),
    avg_se_hat      = mean(se_vec)
  )
}

ppi_ols_population_moments <- function(
  n_big = 100000,
  X_sampler_L,
  f_generator,
  model_type,
  eps_sampler = function(n) rnorm(n, 0, 1),
  a,
  delta = 0
) {
  # Simulate a huge population
  X_all <- X_sampler_L(n_big)          # [n_big x p]
  p <- ncol(X_all)

  # signal with shift along direction a
  f_all <- f_generator(X_all) + delta * as.numeric(X_all %*% a)
  y_all <- f_all + eps_sampler(n_big)

  # two-fold cross-fit to get out-of-fold predictions fhat_all
  fold_id <- sample(rep(1:2, length.out = n_big))
  fhat_all <- numeric(n_big)

  for (fold in 1:2) {
    train_idx <- which(fold_id != fold)
    test_idx  <- which(fold_id == fold)

    fit_obj <- fit_predict_model(
      model_type = model_type,
      X_L = X_all[train_idx, , drop = FALSE],
      y_L = y_all[train_idx],
      X_U = X_all[test_idx, , drop = FALSE]
    )
    fhat_all[test_idx] <- fit_obj$fhat_U
  }

  # Population H = E[XX^T]
  H_hat <- crossprod(X_all, X_all) / n_big  # p x p

  # Estimate beta_star like PPI OLS:
  #    G_hat = E[ X fhat ] + E[ X (Y - fhat) ]
  term_unl <- crossprod(X_all, fhat_all) / n_big
  term_lab <- crossprod(X_all, (y_all - fhat_all)) / n_big
  G_hat    <- as.numeric(term_unl + term_lab)         # length p
  beta_hat <- solve(H_hat, G_hat)                     # p-vector ≈ β^*

  # residuals around β^*
  uY <- as.numeric(y_all - X_all %*% beta_hat)        # Y - Xβ*
  uf <- as.numeric(f_all - X_all %*% beta_hat)        # f(X) - Xβ*
  # NOTE: use f_all not fhat_all

  # build X * residual terms
  X_uY <- sweep(X_all, 1, uY, `*`)  # each row i: X_i * uY_i
  X_uf <- sweep(X_all, 1, uf, `*`)  # each row i: X_i * uf_i

  # Σ matrices
  Sigma_YY <- stats::cov(X_uY)
  Sigma_ff <- stats::cov(X_uf)
  # Σ_Yf is cross-cov between X_uY and X_uf
  # cov(A,B) = E[ (A-meanA)(B-meanB)^T ]
  mu_XuY <- colMeans(X_uY)
  mu_Xuf <- colMeans(X_uf)
  Sigma_Yf <- crossprod( sweep(X_uY, 2, mu_XuY) ,
                         sweep(X_uf, 2, mu_Xuf) ) / (n_big - 1)

  # Precompute the pieces:
  H_inv <- solve(H_hat)
  aHa_inv <- as.numeric(t(a) %*% H_inv) # row vec length p, cache for speed
  A_term <- as.numeric(aHa_inv %*% Sigma_YY %*% H_inv %*% a)
  B_term <- as.numeric(aHa_inv %*% Sigma_ff %*% H_inv %*% a)
  C_term <- as.numeric(aHa_inv %*% Sigma_Yf %*% H_inv %*% a)

  list(
    model_type = model_type,
    H = H_hat,
    H_inv = H_inv,
    beta_star = beta_hat,
    Sigma_YY = Sigma_YY,
    Sigma_ff = Sigma_ff,
    Sigma_Yf = Sigma_Yf,
    A_term = A_term,   # a^T H^{-1} Σ_YY H^{-1} a
    B_term = B_term,   # a^T H^{-1} Σ_ff H^{-1} a
    C_term = C_term    # a^T H^{-1} Σ_Yf H^{-1} a
  )
}

lambda_star_ppi_pp <- function(A_term, B_term, C_term, n, N) {
  r <- n / N
  lambda_star <- C_term / ((1 + r) * B_term)
  lambda_star
}

var_contrast_ppi <- function(A_term, B_term, n, N) {
  (A_term / n) + (B_term / N)
}

var_contrast_ppi_pp <- function(A_term, B_term, C_term, n, N, lambda_val) {
  A_term / n +
    (lambda_val^2) * (B_term / N + B_term / n) -
    (2 * lambda_val / n) * C_term
}

theoretical_power_ols <- function(delta, var_hat, alpha = 0.05) {
  se_hat <- sqrt(var_hat)
  z_alpha <- qnorm(1 - alpha/2)
  z_stat  <- delta / se_hat

  power_val <-
    1 - pnorm(z_alpha - z_stat) +
        pnorm(-z_alpha - z_stat)

  list(
    se = se_hat,
    power = power_val
  )
}



# set.seed(123)
# X_sampler_L <- function(n) matrix(rnorm(n * 3), n, 3)
# f_true <- function(X) 2*X[,1] - X[,2]^2 + 0.5*X[,3]
# a <- c(1, 0, 0)

# ppi_ols_population_moments(
#   X_sampler_L = X_sampler_L,
#   f_generator = f_true,
#   model_type = "glm_correct",
#   a = a,
#   delta = 0.1
# )