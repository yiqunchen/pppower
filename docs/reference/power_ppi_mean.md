# Power or required sample size for the PPI++ mean estimator

Unified interface for power analysis of the PPI++ mean estimator.
Follows the pwr convention: supply **either** `n` or `power` (but not
both). The function computes whichever is left `NULL`.

- `power = NULL` (default): compute power given `n`.

- `n = NULL`: solve for the minimum labeled sample size achieving
  `power`.

You can supply raw variance pieces (`sigma_y2`, `sigma_f2`, `cov_y_f`,
`var_f`, `var_res`) or reuse the metrics interface (`metrics`,
`metric_type`, `m_labeled`, `correction`) to back out the required
inputs.

## Usage

``` r
power_ppi_mean(
  delta,
  N,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f = NULL,
  var_f = NULL,
  var_res = NULL,
  metrics = NULL,
  metric_type = NULL,
  m_labeled = NULL,
  correction = TRUE,
  lambda = NULL,
  lambda_type = c("oracle", "user"),
  lambda_mode = NULL,
  lambda_user = NULL,
  mode = c("error", "cap"),
  warn_smallN = TRUE,
  smallN_threshold = 500
)
```

## Arguments

- delta:

  Effect size \\\theta - \theta_0\\.

- N:

  Unlabeled sample size.

- n:

  Labeled sample size. Set to `NULL` to solve for the required `n`.

- power:

  Target power. Set to `NULL` (default) to compute power from `n`.

- alpha:

  Two-sided significance level (default 0.05).

- sigma_y2:

  Optional outcome variance; overrides anything implied by `metrics`.

- sigma_f2:

  Optional prediction variance; overrides anything implied by `metrics`.

- cov_y_f:

  Optional covariance \\\mathrm{Cov}(Y, f(X))\\. When supplied (directly
  or via `metrics$cov_y_f`) the PPI++ power is returned in addition to
  the PP quantities.

- var_f:

  Optional variance of \\f(X)\\.

- var_res:

  Optional residual variance of residual \\Y - f(X)\\.

- metrics:

  Optional list of predictive-performance summaries.

- metric_type:

  Character string describing the metric bundle (e.g., `"continuous"`,
  `"hard"`, `"prob"`).

- m_labeled:

  Labeled sample size associated with the metrics (defaults to `n`).

- correction:

  Logical; apply finite-sample variance corrections (default `TRUE`).

- lambda:

  Optional user-specified blend weight.

- lambda_type:

  `"oracle"` (default) uses the closed-form blend; `"user"` evaluates
  power at the supplied `lambda`.

- lambda_mode:

  Alias for `lambda_type` used in sample-size mode; one of `"vanilla"`,
  `"oracle"`, or `"user"`.

- lambda_user:

  Scalar specifying \\\lambda\\ when `lambda_mode = "user"` (sample-size
  mode only).

- mode:

  `"error"` (default) or `"cap"`. `"cap"` returns `n = N` with an
  `"achieved_power"` attribute when the required sample size exceeds
  `N`. Only used when solving for `n`.

- warn_smallN:

  Logical. Warn when `N < smallN_threshold` (default `TRUE`).

- smallN_threshold:

  Numeric. Default `500`.

## Value

When computing power: scalar in \[0, 1\]. When computing sample size:
integer `n` required, possibly with `"achieved_power"` attribute if
capped.

## Examples

``` r
# Compute power given n:
power_ppi_mean(delta = 0.3, N = 5000, n = 200,
               sigma_y2 = 1, sigma_f2 = 0.5, cov_y_f = 0.4)
#> [1] 0.9991525

# Compute required n given target power:
power_ppi_mean(delta = 0.3, N = 5000, power = 0.8,
               sigma_y2 = 1, sigma_f2 = 0.5, cov_y_f = 0.4)
#> [1] 60
```
