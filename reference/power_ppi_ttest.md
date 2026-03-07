# Power for PPI++ two-sample t-test (difference in means)

Computes power for a two-sample t-test using PPI++ where we estimate the
difference in means between two independent groups.

## Usage

``` r
power_ppi_ttest(
  delta,
  n_A,
  n_B,
  N_A,
  N_B,
  alpha = 0.05,
  sigma_y2_A,
  sigma_f2_A,
  cov_yf_A,
  sigma_y2_B,
  sigma_f2_B,
  cov_yf_B
)
```

## Arguments

- delta:

  Effect size: the true difference in means \\\mu_A - \mu_B\\.

- n_A:

  Labeled sample size for group A.

- n_B:

  Labeled sample size for group B.

- N_A:

  Unlabeled sample size for group A.

- N_B:

  Unlabeled sample size for group B.

- alpha:

  Two-sided significance level (default 0.05).

- sigma_y2_A:

  Variance of outcomes in group A.

- sigma_f2_A:

  Variance of predictions in group A.

- cov_yf_A:

  Covariance between Y and f in group A.

- sigma_y2_B:

  Variance of outcomes in group B.

- sigma_f2_B:

  Variance of predictions in group B.

- cov_yf_B:

  Covariance between Y and f in group B.

## Value

Scalar power in \\\[0, 1\]\\.

## Details

The variance of the difference estimator is the sum of the PPI++
variances from each group (assuming independence):
\$\$\mathrm{Var}(\hat\mu_A - \hat\mu_B) =
\mathrm{Var}\_{PPI++}(\hat\mu_A) + \mathrm{Var}\_{PPI++}(\hat\mu_B)\$\$

## Examples

``` r
# Two-sample test with equal groups
power_ppi_ttest(
  delta = 0.3,
  n_A = 100, n_B = 100,
  N_A = 5000, N_B = 5000,
  sigma_y2_A = 1, sigma_f2_A = 1, cov_yf_A = 0.7,
  sigma_y2_B = 1, sigma_f2_B = 1, cov_yf_B = 0.7
)
#> [1] 0.8371692
```
