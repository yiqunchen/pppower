# Monte Carlo power for the EIF binary estimator

Monte Carlo power for the EIF binary estimator

## Usage

``` r
simulate_eif_binary(
  R,
  delta,
  N,
  m_cal,
  alpha = 0.05,
  mu0 = NULL,
  mu1 = NULL,
  p_hat = NULL,
  p = NULL,
  sens = NULL,
  spec = NULL,
  seed = 1
)
```

## Arguments

- R:

  Number of Monte Carlo replicates.

- delta:

  Effect size \\\theta - \theta_0\\.

- N:

  Unlabeled/test sample size with surrogate predictions.

- m_cal:

  Calibration sample size.

- alpha:

  Two-sided significance level.

- mu0, mu1:

  Conditional means \\E\[Y \mid \hat Y = 0\]\\ and \\E\[Y \mid \hat Y =
  1\]\\.

- p_hat:

  Prevalence of predicted positives \\P(\hat Y = 1)\\.

- p:

  True prevalence \\P(Y = 1)\\.

- sens, spec:

  Sensitivity and specificity of the surrogate classifier.

- seed:

  RNG seed for reproducibility.

## Value

List with empirical and analytical power plus simulation details.
