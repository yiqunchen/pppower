# Power for PPI++ paired t-test (within-subject differences)

Computes power for a paired t-test using PPI++. This is equivalent to
mean estimation on the differences \\D_i = Y_i^A - Y_i^B\\ with
auxiliary predictions \\f^D_i = f^A(X_i) - f^B(X_i)\\.

## Usage

``` r
power_ppi_paired(
  delta,
  N,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  sigma_D2,
  rho_D,
  sigma_fD2 = NULL,
  pilot_n = NULL,
  confidence = 0.9
)
```

## Arguments

- delta:

  Effect size: the true mean difference \\E\[Y^A - Y^B\]\\.

- N:

  Number of unlabeled pairs with predictions.

- n:

  Number of labeled pairs.

- power:

  Target power. Set to `NULL` (default) to compute power from `n`. Set
  to a value in (0, 1) to solve for `n` instead.

- alpha:

  Two-sided significance level (default 0.05).

- sigma_D2:

  Variance of the differences \\D = Y^A - Y^B\\.

- rho_D:

  Correlation between \\D\\ and \\f^D\\.

- sigma_fD2:

  Optional variance of predicted differences. If NULL, assumed equal to
  `sigma_D2`.

- pilot_n:

  Optional number of labeled pairs in the pilot that estimated `rho_D`.
  If given, the calculation uses the lower confidence bound for `rho_D`
  from
  [`rho_lower_bound()`](https://yiqunchen.github.io/pppower/reference/rho_lower_bound.md)
  instead of the point estimate.

- confidence:

  One-sided confidence level for the bound (default 0.9).

## Value

When `power = NULL`: scalar power in \\\[0, 1\]\\. When `n = NULL`:
required number of labeled pairs (integer).

## Details

The paired t-test reduces to a one-sample problem on the differences.
The PPI++ variance for paired data is:
\$\$\mathrm{Var}(\hat\Delta\_{\lambda^\*}) \approx \frac{\sigma_D^2(1 -
\rho_D^2)}{n}\$\$ where \\\rho_D\\ is the correlation between the
observed difference and the predicted difference.

## Examples

``` r
# Compute power
power_ppi_paired(
  delta = 0.3,
  N = 1000,
  n = 50,
  sigma_D2 = 1,
  rho_D = 0.7
)
#> [1] 0.8276133

# Compute required sample size
power_ppi_paired(
  delta = 0.3,
  N = 1000,
  power = 0.80,
  sigma_D2 = 1,
  rho_D = 0.7
)
#> [1] 47
```
