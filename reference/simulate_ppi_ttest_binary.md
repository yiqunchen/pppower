# Monte Carlo power for binary two-sample PPI++ test

Monte Carlo power for binary two-sample PPI++ test

## Usage

``` r
simulate_ppi_ttest_binary(
  R,
  p_A,
  p_B,
  sens_A,
  spec_A,
  n_A,
  n_B,
  N_A,
  N_B,
  sens_B = sens_A,
  spec_B = spec_A,
  alpha = 0.05,
  delta = NULL,
  seed = 1
)
```

## Arguments

- R:

  Number of Monte Carlo replicates.

- p_A, p_B:

  Group-specific prevalences \\P(Y=1)\\.

- sens_A, spec_A:

  Sensitivity/specificity for group A predictions.

- n_A, n_B:

  Labeled sample sizes.

- N_A, N_B:

  Unlabeled sample sizes.

- sens_B, spec_B:

  Sensitivity/specificity for group B predictions (defaults match group
  A).

- alpha:

  Two-sided significance level.

- delta:

  Optional effect size (defaults to `p_A - p_B`).

- seed:

  RNG seed (set to `NULL` to avoid seeding).

## Value

List with empirical and analytical power plus simulation details.
