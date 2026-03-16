# Monte Carlo vs. analytical power for PPI / PPI++ mean estimation

**\[soft-deprecated\]**

Prefer
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md)
for PPI++ power calculations. This function is retained for backward
compatibility with
[`power_curve_mean()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean.md)
and
[`type1_error_curve_mean()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean.md).

## Usage

``` r
simulate_power(
  delta,
  N,
  n,
  alpha = 0.05,
  R = 1e+05,
  var_f = NULL,
  var_res = NULL,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f = NULL,
  metrics = NULL,
  metric_type = NULL,
  m_labeled = n,
  correction = TRUE
)
```

## Arguments

- delta:

  Effect size \\\theta - \theta_0\\.

- N:

  Unlabeled sample size.

- n:

  Labeled sample size.

- alpha:

  Two-sided significance level.

- R:

  Number of Monte Carlo draws (default 100000).

- var_f:

  Variance of \\f(X)\\.

- var_res:

  Variance of residuals \\Y - f(X)\\.

- sigma_y2:

  Optional outcome variance (needed for PP++ if `cov_y_f` is supplied).

- sigma_f2:

  Optional prediction variance (needed for PP++ if `cov_y_f` is
  supplied).

- cov_y_f:

  Optional covariance \\\operatorname{Cov}(Y, f(X))\\; when present,
  PP++ power is returned.

- metrics:

  Optional predictive summaries to recover missing moments.

- metric_type:

  Character string describing `metrics`.

- m_labeled:

  Labeled sample size associated with `metrics` (defaults to `n`).

- correction:

  Apply finite-sample corrections when recovering moments.

## Value

Named numeric vector with `Exact_PP`, `Empirical_PP`, and when possible
`Exact_PPplus`, `Empirical_PPplus`.
