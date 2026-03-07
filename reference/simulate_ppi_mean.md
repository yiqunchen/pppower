# Monte Carlo power for the PPI++ mean estimator

Monte Carlo power for the PPI++ mean estimator

## Usage

``` r
simulate_ppi_mean(
  R,
  n,
  N,
  alpha = 0.05,
  delta = NULL,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f = NULL,
  var_f = NULL,
  var_res = NULL,
  lambda = NULL,
  lambda_type = c("oracle", "plugin", "user"),
  moments = NULL,
  seed = 1,
  family = stats::gaussian(),
  theta0 = 0
)
```

## Arguments

- R:

  Integer, number of Monte Carlo draws.

- n:

  Labeled sample size.

- N:

  Unlabeled sample size.

- alpha:

  Two-sided test level.

- delta:

  Effect size \\\theta - \theta_0\\.

- sigma_y2:

  Total variance of \\Y\\.

- sigma_f2:

  Variance of \\f(X)\\.

- cov_y_f:

  Covariance between \\Y\\ and \\f(X)\\.

- var_f:

  Variance of \\f(X)\\ (optional, overrides `sigma_f2`).

- var_res:

  Residual variance \\Y - f(X)\\.

- lambda:

  Optional user-specified blend weight.

- lambda_type:

  One of `"oracle"`, `"plugin"`, or `"user"`.

- moments:

  Optional list with entries for population moments.

- seed:

  RNG seed.

- family:

  A [`stats::family()`](https://rdrr.io/r/stats/family.html) object
  ([`gaussian()`](https://rdrr.io/r/stats/family.html) or
  [`binomial()`](https://rdrr.io/r/stats/family.html)).

- theta0:

  Null hypothesis value \\\theta_0\\.

## Value

A list with empirical and analytical power, average standard error,
average lambda, and simulation details.
