# Power or required sample size for PPI++ regression contrast

Computes the minimum labeled sample size \\n\\ for a PPI++ regression
contrast, or computes power given `n`. Follows the pwr convention:
exactly one of `n` or `power` must be `NULL`.

Supports both OLS and canonical GLM contrasts via sandwich/Fisher blocks
computed on labeled and unlabeled samples.

## Usage

``` r
power_ppi_regression(
  delta,
  N,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  c = NULL,
  H_L = NULL,
  H_U = NULL,
  Sigma_YY = NULL,
  Sigma_ff_l = NULL,
  Sigma_ff_u = NULL,
  Sigma_Yf = NULL,
  warn_smallN = TRUE,
  smallN_threshold = 500,
  mode = c("error", "cap")
)
```

## Arguments

- delta:

  Numeric scalar. Effect size \\a^\top \beta - \theta_0\\. Must be
  nonzero when solving for `n`.

- N:

  Positive scalar. Number of unlabeled observations.

- n:

  Labeled sample size. Set to `NULL` to solve for required `n`.

- power:

  Desired power. Set to `NULL` to compute power from `n`.

- alpha:

  Significance level (default `0.05`).

- lambda_mode:

  `"vanilla"`, `"oracle"`, or `"user"`.

- lambda_user:

  Scalar specifying \\\lambda\\ when `lambda_mode = "user"`.

- c:

  Contrast vector \\a\\.

- H_L:

  Labeled Hessian/Fisher block.

- H_U:

  Unlabeled Hessian/Fisher block.

- Sigma_YY:

  Labeled covariance block.

- Sigma_ff_l:

  Labeled prediction covariance block.

- Sigma_ff_u:

  Unlabeled prediction covariance block.

- Sigma_Yf:

  Cross-covariance block.

- warn_smallN:

  Logical. Warn when `N < smallN_threshold`.

- smallN_threshold:

  Numeric. Default `500`.

- mode:

  `"error"` or `"cap"`.

## Value

When computing power: scalar in \[0, 1\]. When computing sample size:
integer `n` required.

## Examples

``` r
p <- 3
cvec <- c(1, 0, -1)
Hpop <- diag(p)
SYY  <- diag(p)
Sff  <- 0.5 * diag(p)
SYf  <- 0.3 * diag(p)

# Compute required sample size
power_ppi_regression(
  delta = 0.15, N = 2000, power = 0.90,
  c = cvec, H_L = Hpop, H_U = Hpop,
  Sigma_YY = SYY, Sigma_ff_l = Sff, Sigma_ff_u = Sff,
  Sigma_Yf = SYf, lambda_mode = "oracle"
)
#> [1] 815
```
