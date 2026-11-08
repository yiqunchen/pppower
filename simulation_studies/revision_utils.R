# =============================================================================
# Utilities for the R1 revision simulations (Setting 23)
# =============================================================================
# Data generation, population moments, plug-in variance factor (kappa), and
# vectorized PPI++ estimators for several choices of the tuning parameter:
#   oracle     lambda* from population moments
#   plugin     lambda_hat from the labeled data (current pppower default)
#   plugin_adj plugin estimate, variance inflated by the second-order plug-in
#              term c * s_e2 * kappa_hat / (n (n - 3)), t_{n-2} reference
#   crossfit   2-fold cross-fitted lambda_hat
#   pilot      lambda fixed from an independent pilot sample of size m
# Every function works on R replicates at once (one column per replicate).

# ---------------------------------------------------------------------------
# Distributions
# ---------------------------------------------------------------------------
# A distribution spec is a list with `type` in
#   "gaussian", "lognormal", "t5", "binary"
# and the fields used by that type (rho, sigma_log, p, sens, spec).
# Continuous types are generated with mean zero shift and then shifted by
# `shift`; binary types use prevalence `p`.

lognormal_rho_log <- function(rho, sigma_log) {
  # log-scale correlation that gives original-scale correlation rho
  log(1 + rho * (exp(sigma_log^2) - 1)) / sigma_log^2
}

gen_pairs <- function(dist, n, R, shift = 0) {
  if (dist$type == "binary") {
    Y <- matrix(rbinom(n * R, 1, dist$p), n, R)
    f <- ifelse(Y == 1, rbinom(n * R, 1, dist$sens), rbinom(n * R, 1, 1 - dist$spec))
    return(list(Y = Y, f = matrix(f, n, R)))
  }
  z1 <- matrix(rnorm(n * R), n, R)
  z2 <- matrix(rnorm(n * R), n, R)
  if (dist$type == "gaussian") {
    Y <- z1
    f <- dist$rho * z1 + sqrt(1 - dist$rho^2) * z2
  } else if (dist$type == "lognormal") {
    rl <- lognormal_rho_log(dist$rho, dist$sigma_log)
    Y <- exp(dist$sigma_log * z1)
    f <- exp(dist$sigma_log * (rl * z1 + sqrt(1 - rl^2) * z2))
  } else if (dist$type == "t5") {
    rg <- dist$rho_gauss
    Y <- qt(pnorm(z1), df = 5)
    f <- qt(pnorm(rg * z1 + sqrt(1 - rg^2) * z2), df = 5)
  } else stop("unknown type")
  list(Y = Y + shift, f = f + shift)
}

# Mean of N unlabeled predictions, one value per replicate
gen_fU_bar <- function(dist, N, R, shift = 0, chunk = 4e6) {
  if (dist$type == "binary") {
    p_f <- dist$sens * dist$p + (1 - dist$spec) * (1 - dist$p)
    return(rbinom(R, N, p_f) / N)
  }
  if (dist$type == "gaussian") return(rnorm(R, shift, sqrt(1 / N)))
  out <- numeric(R)
  per <- max(1, floor(chunk / N))
  idx <- split(seq_len(R), ceiling(seq_len(R) / per))
  for (ii in idx) out[ii] <- colMeans(gen_pairs(dist, N, length(ii), shift)$f)
  out
}

# ---------------------------------------------------------------------------
# Population moments and the plug-in factor kappa
# ---------------------------------------------------------------------------
# kappa (linear conditional mean approximation):
#   kappa = 2 - E[(f-mu)^2 e^2] / (sf2 se2) + 2 Cov(f, e^2) E[(f-mu)^3] / (sf2^2 se2)
# where e is the residual of the linear projection of Y on f.
kappa_from_sample <- function(Y, f) {
  fc <- f - mean(f)
  b <- sum((Y - mean(Y)) * fc) / sum(fc^2)
  e <- Y - mean(Y) - b * fc
  sf2 <- mean(fc^2); se2 <- mean(e^2); m3 <- mean(fc^3)
  2 - mean(fc^2 * e^2) / (sf2 * se2) + 2 * mean(fc * (e^2 - se2)) * m3 / (sf2^2 * se2)
}

kappa_binary <- function(p, sens, spec) {
  p_f <- sens * p + (1 - spec) * (1 - p)
  mu1 <- sens * p / p_f                 # P(Y = 1 | f = 1)
  mu0 <- (1 - sens) * p / (1 - p_f)     # P(Y = 1 | f = 0)
  t1 <- mu1 * (1 - mu1); t0 <- mu0 * (1 - mu0)
  se2 <- (1 - p_f) * t0 + p_f * t1
  ((1 - p_f) * t1 + p_f * t0) / se2
}

pop_moments <- function(dist, shift = 0, big = 4e6) {
  if (dist$type == "binary") {
    m <- pppower::binary_moments_from_sens_spec(p = dist$p, sens = dist$sens, spec = dist$spec)
    sY2 <- m$sigma_y2; sf2 <- m$sigma_f2; cv <- m$cov_y_f
    kap <- kappa_binary(dist$p, dist$sens, dist$spec)
    theta <- dist$p
  } else if (dist$type == "gaussian") {
    sY2 <- 1; sf2 <- 1; cv <- dist$rho; kap <- 1; theta <- shift
  } else if (dist$type == "lognormal") {
    s2 <- dist$sigma_log^2
    sY2 <- exp(s2) * (exp(s2) - 1); sf2 <- sY2; cv <- dist$rho * sY2
    theta <- exp(s2 / 2) + shift
    set.seed(99); d <- gen_pairs(dist, big, 1)
    kap <- kappa_from_sample(d$Y, d$f)
  } else {
    set.seed(99); d <- gen_pairs(dist, big, 1)
    sY2 <- var(d$Y[, 1]); sf2 <- var(d$f[, 1]); cv <- cov(d$Y[, 1], d$f[, 1])
    kap <- kappa_from_sample(d$Y, d$f); theta <- shift
  }
  # kappa_plan: the formula assumes a linear E[Y | f]; when it does not hold
  # (e.g. strongly skewed log-normal outcomes) it can turn negative, so we
  # never plan with less than the Gaussian value of 1.
  list(theta = theta, sY2 = sY2, sf2 = sf2, cov = cv,
       rho = cv / sqrt(sY2 * sf2), se2 = sY2 - cv^2 / sf2, kappa = kap,
       kappa_plan = max(1, kap))
}

# t5 copula: Gaussian correlation that gives the target t5 correlation
t5_rho_gauss <- function(rho, big = 2e6) {
  set.seed(7); z1 <- rnorm(big); z2 <- rnorm(big)
  g <- function(rg) cor(qt(pnorm(z1), 5), qt(pnorm(rg * z1 + sqrt(1 - rg^2) * z2), 5)) - rho
  uniroot(g, c(0.01, 0.999))$root
}

# ---------------------------------------------------------------------------
# Analytical variances and power
# ---------------------------------------------------------------------------
var_opt <- function(mom, n, N) mom$sY2 / n - mom$cov^2 / mom$sf2 * N / (n * (n + N))

var_plugin <- function(mom, n, N) {
  cc <- N / (n + N)
  var_opt(mom, n, N) + cc * mom$se2 * mom$kappa_plan / (n * (n - 3))
}

wald_power <- function(delta, V, alpha = 0.05) {
  z <- qnorm(1 - alpha / 2)
  pnorm(-z + abs(delta) / sqrt(V)) + pnorm(-z - abs(delta) / sqrt(V))
}

# smallest n with power >= target, for a variance function V(n)
plan_n <- function(Vfun, delta, target = 0.8, alpha = 0.05, n_min = 5, n_max = 1e5) {
  pw <- function(n) wald_power(delta, Vfun(n), alpha) - target
  if (pw(n_max) < 0) return(NA_integer_)
  if (pw(n_min) >= 0) return(as.integer(n_min))
  ceiling(uniroot(pw, c(n_min, n_max))$root)
}

# ---------------------------------------------------------------------------
# Vectorized estimators
# ---------------------------------------------------------------------------
col_moments <- function(Y, f) {
  n <- nrow(Y)
  Yb <- colMeans(Y); fb <- colMeans(f)
  Yc <- Y - rep(Yb, each = n); fc <- f - rep(fb, each = n)
  list(n = n, Yb = Yb, fb = fb,
       sY2 = colSums(Yc^2) / (n - 1), sf2 = colSums(fc^2) / (n - 1),
       cYf = colSums(Yc * fc) / (n - 1), Yc = Yc, fc = fc)
}

ppi_var_hat <- function(m, lam, N) {
  pmax(m$sY2 / m$n + lam^2 * m$sf2 * (1 / N + 1 / m$n) - 2 * lam * m$cYf / m$n, 1e-12)
}

safe_lambda <- function(cYf, sf2, n, N) {
  lam <- cYf / ((1 + n / N) * sf2)
  lam[!is.finite(lam)] <- 0
  lam
}

# kappa_hat per replicate, from labeled data
kappa_hat_cols <- function(m) {
  b <- m$cYf / m$sf2; b[!is.finite(b)] <- 0
  e <- m$Yc - rep(b, each = m$n) * m$fc
  sf2b <- colMeans(m$fc^2); se2b <- colMeans(e^2); m3 <- colMeans(m$fc^3)
  k <- 2 - colMeans(m$fc^2 * e^2) / (sf2b * se2b) +
    2 * colMeans(m$fc * (e^2 - rep(se2b, each = m$n))) * m3 / (sf2b^2 * se2b)
  k[!is.finite(k)] <- 1
  pmin(pmax(k, 1), 50)
}

# Returns list of estimates and variances for each method (R-vectors).
ppi_methods <- function(Y, f, fU_bar, N, lam_oracle, lam_pilot) {
  m <- col_moments(Y, f); n <- m$n; cc <- N / (n + N)
  D <- fU_bar - m$fb
  out <- list()

  out$oracle <- list(th = m$Yb + lam_oracle * D, v = ppi_var_hat(m, lam_oracle, N), df = Inf)

  lam_p <- safe_lambda(m$cYf, m$sf2, n, N)
  th_p <- m$Yb + lam_p * D
  v_p <- ppi_var_hat(m, lam_p, N)
  out$plugin <- list(th = th_p, v = v_p, df = Inf)

  # plug-in variance correction: c * s_e2 * kappa_hat / (n (n - 3))
  se2_hat <- pmax(m$sY2 - m$cYf^2 / m$sf2, 0); se2_hat[!is.finite(se2_hat)] <- m$sY2[!is.finite(se2_hat)]
  se2_hat <- se2_hat * (n - 1) / (n - 2)
  k_hat <- kappa_hat_cols(m)
  out$plugin_adj <- list(th = th_p, v = v_p + cc * se2_hat * k_hat / (n * (n - 3)), df = n - 2)

  # 2-fold cross-fitting: lambda for each half is estimated on the other half
  h <- floor(n / 2); i1 <- seq_len(h); i2 <- (h + 1):n
  m1 <- col_moments(Y[i1, , drop = FALSE], f[i1, , drop = FALSE])
  m2 <- col_moments(Y[i2, , drop = FALSE], f[i2, , drop = FALSE])
  l1 <- safe_lambda(m2$cYf, m2$sf2, n, N)    # used on fold 1
  l2 <- safe_lambda(m1$cYf, m1$sf2, n, N)    # used on fold 2
  w1 <- h / n; w2 <- 1 - w1
  th_cf <- w1 * (m1$Yb + l1 * (fU_bar - m1$fb)) + w2 * (m2$Yb + l2 * (fU_bar - m2$fb))
  part <- function(mk, lk) pmax(mk$sY2 / mk$n + lk^2 * mk$sf2 / mk$n - 2 * lk * mk$cYf / mk$n, 0)
  v_cf <- w1^2 * part(m1, l1) + w2^2 * part(m2, l2) + (w1 * l1 + w2 * l2)^2 * m$sf2 / N
  out$crossfit <- list(th = th_cf, v = pmax(v_cf, 1e-12), df = Inf)

  out$pilot <- list(th = m$Yb + lam_pilot * D, v = ppi_var_hat(m, lam_pilot, N), df = Inf)

  out$classical <- list(th = m$Yb, v = m$sY2 / n, df = n - 1)
  out
}

pilot_lambda <- function(dist, m_pilot, R, n, N, shift = 0) {
  d <- gen_pairs(dist, m_pilot, R, shift)
  mm <- col_moments(d$Y, d$f)
  safe_lambda(mm$cYf, mm$sf2, n, N)
}

# Summaries for one method: rejection, bias, variance, coverage
summarize_method <- function(est, theta0, theta_true, alpha = 0.05) {
  df <- est$df; df[!is.finite(df) & !is.infinite(df)] <- Inf   # NaN df (zero variance) -> normal
  crit <- qt(1 - alpha / 2, df)              # df = Inf gives the normal quantile
  se <- sqrt(est$v)
  # a replicate with zero estimated variance (e.g. all labels 0) is counted
  # as a non-rejection
  z <- abs(est$th - theta0) / se
  z[!is.finite(z)] <- 0
  data.frame(
    reject = mean(z > crit),
    bias = mean(est$th) - theta_true,
    mc_var = var(est$th),
    mean_vhat = mean(est$v),
    coverage = mean(abs(est$th - theta_true) <= crit * se)
  )
}
