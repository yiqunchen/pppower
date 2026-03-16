# Construct Hessian/Fisher and Covariance Blocks for PPI/PPI++ Regression

Computes all required Hessian/Fisher information matrices and covariance
blocks used by Prediction-Powered Inference (PPI) and PPI++ for
regression-based estimands (OLS or GLM contrasts).

This helper wraps:

- `compute_hessian_fisher` — model Hessian or Fisher information

- `compute_sigma_blocks` — Sigma-block covariance components

and produces a unified object that can be passed directly to
[`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)
for labeled sample size calculations.

## Usage

``` r
compute_ppi_blocks(
  model_type = c("ols", "glm"),
  X_l,
  Y_l,
  f_l,
  X_u,
  f_u,
  beta = NULL,
  family = NULL
)
```

## Arguments

- model_type:

  Character string: `"ols"` or `"glm"`. Determines whether the estimator
  is linear regression (OLS) or a GLM.

- X_l:

  Matrix of labeled covariates (n x p).

- Y_l:

  Vector of labeled responses (required for OLS and GLM).

- f_l:

  Vector of model predictions on labeled data.

- X_u:

  Matrix of unlabeled covariates (N x p).

- f_u:

  Vector of model predictions on unlabeled data.

- beta:

  Numeric vector of regression coefficients used for Hessian/Fisher
  evaluation. Required for GLM; optional for OLS.

- family:

  GLM family (`"binomial"` or `"gaussian"`). Only used for
  `model_type = "glm"`.

## Value

A named list of matrices:

- H_L:

  Labeled-data Hessian / Fisher information (p x p)

- H_U:

  Unlabeled-data Hessian / Fisher information (p x p)

- Sigma_YY:

  Covariance of labeled score (`Y - X beta`) (p x p)

- Sigma_ff_l:

  Covariance of prediction score on labeled data (p x p)

- Sigma_ff_u:

  Covariance of prediction score on unlabeled data (p x p)

- Sigma_Yf:

  Cross-covariance between labeled scores and prediction scores (p x p)

These matrices are exactly the inputs needed for
[`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)
in `type = "regression"` mode.

## See also

[power_ppi_regression](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)

## Examples

``` r
set.seed(1)
p <- 3
n_l <- 400
n_u <- 2000

X_l <- matrix(rnorm(n_l * p), n_l, p)
X_u <- matrix(rnorm(n_u * p), n_u, p)
beta <- c(1, -0.5, 0.3)

# Labeled response and predictions
Y_l <- drop(X_l %*% beta + rnorm(n_l))
f_l <- drop(X_l %*% beta + rnorm(n_l, sd = 0.3))
f_u <- drop(X_u %*% beta + rnorm(n_u, sd = 0.3))

blocks <- compute_ppi_blocks(
  model_type = "ols",
  X_l = X_l, Y_l = Y_l, f_l = f_l,
  X_u = X_u, f_u = f_u,
  beta = beta
)

# Use in sample size solver:
c_vec <- c(1, 0, 0)
delta <- as.numeric(t(c_vec) %*% beta)

power_ppi_regression(
  delta = delta, N = n_u, power = 0.80,
  lambda_mode = "oracle",
  c = c_vec,
  H_L = blocks$H_L, H_U = blocks$H_U,
  Sigma_YY = blocks$Sigma_YY,
  Sigma_ff_l = blocks$Sigma_ff_l,
  Sigma_ff_u = blocks$Sigma_ff_u,
  Sigma_Yf = blocks$Sigma_Yf
)
#> [1] 8
```
