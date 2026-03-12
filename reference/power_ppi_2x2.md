# Power or required sample size for a 2x2 table with PPI++ surrogates

Computes power or solves for the required labeled sample size for
testing an odds ratio (OR) or relative risk (RR) from a 2x2 contingency
table when a binary PPI++ surrogate is available.

Follows the pwr convention: exactly one of `n` or `power` must be
`NULL`.

**Odds ratio** (default): Uses logistic regression with design matrix
\\\[1, X\]\\ where \\X \in \\0,1\\\\ is the group indicator. The effect
is \\\beta_1 = \log(\text{OR})\\. Analytical variance blocks are
computed from the cell probabilities and surrogate performance, then
passed to
[`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md).

**Relative risk**: Uses the delta method on per-group PPI++ mean
estimates. The test statistic is a Wald test for \\\log(\hat p_1 / \hat
p_0)\\.

## Usage

``` r
power_ppi_2x2(
  p_exp,
  p_ctrl,
  N,
  n = NULL,
  power = NULL,
  alpha = 0.05,
  prev_exp = 0.5,
  sens,
  spec,
  effect_measure = c("OR", "RR"),
  lambda_mode = "oracle",
  lambda_user = NULL,
  warn_smallN = TRUE,
  ...
)
```

## Arguments

- p_exp:

  Numeric scalar. Probability of outcome in the exposed group, \\P(Y = 1
  \mid X = 1)\\.

- p_ctrl:

  Numeric scalar. Probability of outcome in the control group, \\P(Y = 1
  \mid X = 0)\\.

- N:

  Positive integer. Total number of unlabeled observations.

- n:

  Total number of labeled observations. Set to `NULL` to solve for the
  required `n`.

- power:

  Target power. Set to `NULL` to compute power from `n`.

- alpha:

  Significance level (default `0.05`).

- prev_exp:

  Numeric scalar. Prevalence of exposure, \\P(X = 1)\\. Default `0.5`
  (balanced groups).

- sens:

  Sensitivity of the binary surrogate classifier.

- spec:

  Specificity of the binary surrogate classifier.

- effect_measure:

  Character. `"OR"` (odds ratio, default) or `"RR"` (relative risk).

- lambda_mode:

  `"oracle"` (default), `"vanilla"`, or `"user"`.

- lambda_user:

  Scalar lambda when `lambda_mode = "user"`.

- warn_smallN:

  Logical. Warn when `N` is small.

- ...:

  Additional arguments (currently unused).

## Value

When computing power: scalar in \[0, 1\]. When computing sample size:
integer `n` required.

## Examples

``` r
# OR test: exposed group has higher risk
power_ppi_2x2(
  p_exp = 0.40, p_ctrl = 0.20,
  N = 500, n = 80, alpha = 0.05,
  sens = 0.85, spec = 0.85,
  effect_measure = "OR"
)
#> [1] 0.6631492

# Solve for required n
power_ppi_2x2(
  p_exp = 0.40, p_ctrl = 0.20,
  N = 500, power = 0.80, alpha = 0.05,
  sens = 0.85, spec = 0.85,
  effect_measure = "OR"
)
#> [1] 115

# RR test
power_ppi_2x2(
  p_exp = 0.40, p_ctrl = 0.20,
  N = 500, n = 80, alpha = 0.05,
  sens = 0.85, spec = 0.85,
  effect_measure = "RR"
)
#> [1] 0.6410825
```
