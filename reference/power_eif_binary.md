# Power or required calibration size for EIF estimator with binary surrogates

Unified interface for the efficient influence function (EIF) estimator
used with binary surrogate predictions (e.g., LLM-as-a-judge hard
votes). Follows the pwr convention: supply **either** `m_cal` or `power`
(but not both). The function computes whichever is left `NULL`.

You can provide either direct conditional means `(mu0, mu1, p_hat)` or
sensitivity/specificity plus prevalence `(p, sens, spec)`.

## Usage

``` r
power_eif_binary(
  delta,
  N,
  m_cal = NULL,
  power = NULL,
  alpha = 0.05,
  mu0 = NULL,
  mu1 = NULL,
  p_hat = NULL,
  p = NULL,
  sens = NULL,
  spec = NULL
)
```

## Arguments

- delta:

  Effect size \\\theta - \theta_0\\.

- N:

  Unlabeled/test sample size with surrogate predictions.

- m_cal:

  Calibration (labeled) sample size. Set to `NULL` to solve for the
  required `m_cal`.

- power:

  Target power. Set to `NULL` (default) to compute power from `m_cal`.

- alpha:

  Two-sided significance level (default 0.05).

- mu0, mu1:

  Conditional means \\E\[Y \mid \hat Y = 0\]\\ and \\E\[Y \mid \hat Y =
  1\]\\.

- p_hat:

  Prevalence of predicted positives \\P(\hat Y = 1)\\.

- p:

  True prevalence \\P(Y = 1)\\.

- sens, spec:

  Sensitivity and specificity of the surrogate classifier.

## Value

When computing power: scalar in \[0, 1\] with attributes `theta`,
`p_hat`, `var_unlabeled`, and `var_cal`. When computing sample size:
integer required calibration size.
