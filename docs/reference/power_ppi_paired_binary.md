# Power for PPI++ paired test with binary outcomes

Computes power for a paired test when outcomes are binary and
predictions are derived from a classifier with sensitivity and
specificity. The paired differences \\D = Y^A - Y^B\\ and prediction
differences \\f^D = f^A - f^B\\ are treated via a CLT approximation
using population moments implied by a joint Bernoulli model with
within-pair correlation.

## Usage

``` r
power_ppi_paired_binary(
  delta,
  N,
  n,
  alpha = 0.05,
  p_A,
  p_B,
  rho_within,
  sens,
  spec
)
```

## Arguments

- delta:

  Effect size: mean difference \\E\[Y^A - Y^B\]\\. If NULL, uses \\p_A -
  p_B\\.

- N:

  Number of unlabeled pairs with predictions.

- n:

  Number of labeled pairs.

- alpha:

  Two-sided significance level (default 0.05).

- p_A:

  Marginal probability \\P(Y^A=1)\\.

- p_B:

  Marginal probability \\P(Y^B=1)\\.

- rho_within:

  Pearson correlation between \\Y^A\\ and \\Y^B\\.

- sens:

  Classifier sensitivity.

- spec:

  Classifier specificity.

## Value

Scalar power in \\\[0, 1\]\\.

## Details

The joint probability \\P(Y^A=1, Y^B=1)\\ is set to \$\$p\_{11} = p_A
p_B + \rho\_{\text{within}} \sqrt{p_A(1-p_A)p_B(1-p_B)}\$\$ and clipped
to the feasible range. Predictions are generated conditionally on \\Y\\
with sensitivity/specificity, and population moments are computed by
marginalizing over the joint distribution of \\(Y^A, Y^B)\\.

## Examples

``` r
power_ppi_paired_binary(
  delta = 0.08,
  N = 5000,
  n = 500,
  p_A = 0.34,
  p_B = 0.26,
  rho_within = 0.3,
  sens = 0.85,
  spec = 0.85
)
#> [1] 0.9809592
```
