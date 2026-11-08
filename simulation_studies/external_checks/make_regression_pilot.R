# Pilot data and pppower regression sample sizes for the ppi_py comparison
# (check_regression_ppi_py.py). Writes pilot CSVs and pppower_regression_nstar.csv.
# Usage: Rscript simulation_studies/external_checks/make_regression_pilot.R
suppressPackageStartupMessages(library(pppower))
set.seed(2026)
out_dir <- file.path("simulation_studies", "external_checks")
M <- 100000
z <- qnorm(0.975); zb <- qnorm(0.8)
rows <- list()

make_ols <- function(m, rho) {
  X <- cbind(1, rnorm(m), rnorm(m)); eps <- rnorm(m)
  Y <- 0.3 * X[, 2] + eps
  f <- 0.3 * X[, 2] + rho * eps + sqrt(1 - rho^2) * rnorm(m)
  list(X = X, Y = Y, f = f)
}
make_logit <- function(m, acc) {
  X <- cbind(1, rnorm(m), rnorm(m))
  Y <- rbinom(m, 1, plogis(0.5 * X[, 2]))
  f <- ifelse(runif(m) < acc, Y, 1 - Y)
  list(X = X, Y = Y, f = f)
}

for (cfg in list(list("ols", 0.5), list("ols", 0.8), list("logistic", 0.8), list("logistic", 0.9))) {
  type <- cfg[[1]]; q <- cfg[[2]]
  gen <- if (type == "ols") make_ols else make_logit
  lab <- gen(M, q); unl <- gen(M, q)
  if (type == "ols") {
    beta <- qr.solve(lab$X, lab$Y)
    blk <- compute_ppi_blocks("ols", lab$X, lab$Y, lab$f, unl$X, unl$f, beta = beta)
  } else {
    beta <- coef(glm(lab$Y ~ lab$X - 1, family = binomial()))
    blk <- compute_ppi_blocks("glm", lab$X, lab$Y, lab$f, unl$X, unl$f, beta = beta, family = "binomial")
  }
  tag <- sprintf("%s_%s", type, q)
  write.csv(data.frame(x1 = lab$X[, 2], x2 = lab$X[, 3], Y = lab$Y, f = lab$f),
            file.path(out_dir, sprintf("pilot_%s.csv", tag)), row.names = FALSE)
  a <- c(0, 1, 0)                                # coefficient of x1
  Hinv <- solve(blk$H_L)
  var_cl <- drop(t(a) %*% Hinv %*% blk$Sigma_YY %*% Hinv %*% a)
  delta <- if (type == "ols") 0.3 else 0.5
  n_cl <- ceiling(var_cl / (delta / (z + zb))^2)
  for (N in c(500, 2000, 10000)) {
    n_star <- ceiling(suppressWarnings(power_ppi_regression(
      delta = delta, N = N, n = NULL, power = 0.8, lambda_mode = "oracle", c = a,
      H_L = blk$H_L, H_U = blk$H_U, Sigma_YY = blk$Sigma_YY,
      Sigma_ff_l = blk$Sigma_ff_l, Sigma_ff_u = blk$Sigma_ff_u, Sigma_Yf = blk$Sigma_Yf)))
    rows[[length(rows) + 1]] <- data.frame(tag = tag, type = type, quality = q, N = N,
                                           n_cl = n_cl, n_star = n_star)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "pppower_regression_nstar.csv"), row.names = FALSE)
