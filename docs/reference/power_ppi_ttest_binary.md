# Power for PPI++ two-sample test with binary outcomes (sens/spec input)

Computes PPI++ power for a two-sample difference in proportions when
predictions come from a binary classifier summarized by sensitivity and
specificity. This helper converts the classification metrics into the
required moment inputs and delegates to
[`power_ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest.md).

## Usage

``` r
power_ppi_ttest_binary(
  p_A,
  p_B,
  n_A,
  n_B,
  N_A,
  N_B,
  sens_A,
  spec_A,
  sens_B = sens_A,
  spec_B = spec_A,
  delta = NULL,
  alpha = 0.05
)
```

## Arguments

- p_A, p_B:

  Group-specific prevalences \\P(Y=1)\\.

- n_A, n_B:

  Labeled sample sizes per group.

- N_A, N_B:

  Unlabeled sample sizes per group.

- sens_A, spec_A:

  Sensitivity and specificity for group A predictions.

- sens_B, spec_B:

  Sensitivity and specificity for group B predictions (defaults match
  group A).

- delta:

  Optional effect size \\p_A - p_B\\. If `NULL`, it is computed from
  `p_A` and `p_B`.

- alpha:

  Two-sided significance level (default 0.05).

## Value

Scalar power in \[0, 1\] with a `"moments"` attribute containing the
per-group moment lists.
