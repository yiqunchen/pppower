# pppower

**pppower** is an R package for power analysis and sample-size
calculation under the Prediction-Powered Inference (PPI/`PPI++`)
framework.

## Prevalence Estimation with PPI

For binary outcomes (prevalence), you can estimate required **labeled**
sample size while leveraging model predictions on a large unlabeled set.
This is exactly what PPI is for:

- Inputs: target effect size (), unlabeled size (N), and model quality
  (e.g., sensitivity/specificity)
- Output: required labeled (n) for your target power
- Benefit: better predictions can reduce labeled annotation burden

Use `power_ppi_mean(..., n = NULL)` for one-sample prevalence and
[`power_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest_binary.md)
for two-group prevalence differences.

## Interactive Calculator

[Open Interactive PPI Sample Size
Calculator](https://yiqunchen.github.io/pppower/articles/sample-size-calculator.html)

``` R
Direct link to the live calculator page.
```

## What is Prediction-Powered Inference?

Many modern studies have access to:

- **Labeled data** $`(X_i, Y_i)`$ for $`i = 1, \ldots, n`$ — expensive
  to obtain (human annotations, lab assays, expert reviews)
- **Unlabeled data** $`\tilde{X}_j`$ for $`j = 1, \ldots, N`$ — cheap
  and abundant, with ML predictions $`f(\tilde{X}_j)`$ available

PPI combines both sources for valid, efficient inference. The `PPI++`
estimator for the population mean $`\theta^* = E[Y]`$ is:

``` math
\hat{\theta}_\lambda = \bar{Y}_L + \lambda(\bar{f}_U - \bar{f}_L)
```

where $`\lambda`$ controls the blend between labeled data and
predictions. The key insight: **good predictions reduce the amount of
labeled data needed** for a given level of statistical power.

## Key Formulas

The `PPI++` variance under the oracle $`\lambda^*`$ is:

``` math
\text{Var}(\hat{\theta}_{\lambda^*}) = \frac{\sigma_Y^2}{n} - \frac{\text{Cov}(Y,f)^2}{\sigma_f^2} \cdot \frac{N}{n(n+N)}
```

Solving for the required labeled sample size (quadratic in $`n`$):

``` math
n \geq \frac{\sigma_Y^2 - S^2 N + \sqrt{(\sigma_Y^2 - S^2 N)^2 + 4 S^2 N \sigma_Y^2 (1 - \rho_{Yf}^2)}}{2S^2}
```

where $`S^2 = (\Delta / (z_{1-\alpha/2} + z_{1-\beta}))^2`$.

**Rule of thumb:** The sample size ratio satisfies
$`n_{\texttt{PPI++}} / n_{\text{classical}} \approx 1 - R^2`$, where
$`R^2`$ is the squared correlation between $`Y`$ and $`f(X)`$. A
predictor explaining 80% of the variance cuts labeled data needs by
~80%.

## Supported Designs

| Design | Power / Sample Size | Simulation |
|----|----|----|
| One-sample mean (continuous) | [`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md) | [`simulate_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_mean.md) |
| One-sample mean (binary) | [`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md) | `simulate_ppi_vanilla_mean()` |
| Two-sample $`t`$-test | [`power_ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest.md) | [`simulate_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_ttest_binary.md) |
| Paired $`t`$-test | [`power_ppi_paired()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired.md) | — |
| 2x2 contingency table | [`power_ppi_2x2()`](https://yiqunchen.github.io/pppower/reference/power_ppi_2x2.md) | — |
| Regression contrast | [`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md) | — |

## Installation

You can install the development version of `pppower` from
[GitHub](https://github.com/) with:

``` r

# install.packages("devtools")
devtools::install_github("yiqunchen/pppower")
```

## Quick Example

A three-step workflow: compute power → solve for required $`n`$ → verify
with Monte Carlo.

``` r

library(pppower)

# Step 1: Compute power for a given sample size
power_ppi_mean(
  delta    = 0.2,
  N        = 5000,
  n        = 200,
  var_f    = 0.6,
  var_res  = 0.4,
  cov_y_f  = 0.6,
  alpha    = 0.05
)

# Step 2: Solve for required labeled n to achieve 80% power
power_ppi_mean(
  delta       = 0.2,
  N           = 5000,
  n           = NULL,
  power       = 0.80,
  sigma_y2    = 1.0,
  sigma_f2    = 0.49,
  cov_y_f     = 0.63,
  lambda_mode = "oracle"
)

# Step 3: Verify with Monte Carlo simulation
simulate_ppi_mean(
  R     = 2000,
  n     = 100,
  N     = 5000,
  delta = 0.2,
  var_f   = 0.49,
  var_res = 0.51,
  alpha = 0.05
)
```

## Input Modes

`pppower` supports three ways to specify variance components:

**1. Direct moments** — supply $`\sigma_Y^2`$, $`\sigma_f^2`$, and
$`\text{Cov}(Y, f)`$ directly:

``` r

power_ppi_mean(delta = 0.2, N = 5000, n = 200,
                  sigma_y2 = 1.0, sigma_f2 = 0.49, cov_y_f = 0.63)
```

**2. Prediction variance / residual variance** — decompose outcome
variance via `var_f` and `var_res`:

``` r

power_ppi_mean(delta = 0.2, N = 5000, n = 200,
                  var_f = 0.49, var_res = 0.51, cov_y_f = 0.49)
```

**3. Metrics interface** — supply sensitivity/specificity (binary) or
MSE/$`R^2`$ (continuous):

``` r

power_ppi_mean(delta = 0.05, N = 5000, n = 200,
                  metrics = list(sensitivity = 0.85, specificity = 0.90,
                                 p_y = 0.3, m_obs = 200),
                  metric_type = "classification")
```

## Vignettes

- **Quickstart: Prediction-Powered Power & Sample Size** — end-to-end
  walkthrough of power curves and Type I error checks
  ([`vignette("intro-ppi")`](https://yiqunchen.github.io/pppower/articles/intro-ppi.md))
- **Detailed Dive: Sample Size & Metrics Interface** — how
  `power_ppi_mean(..., n = NULL)` works across input modes, including
  regression contrasts
  ([`vignette("ppi-sample-size")`](https://yiqunchen.github.io/pppower/articles/ppi-sample-size.md))
- **Detailed Dive: Variance Formulas and lambda-star** — mathematical
  derivations behind the `PPI++` variance and oracle $`\lambda^*`$
  ([`vignette("deep-dive-math")`](https://yiqunchen.github.io/pppower/articles/deep-dive-math.md))
- **2x2 Contingency Table Calculator** — odds-ratio and relative-risk
  planning with binary surrogates
  ([`vignette("calculator-2x2")`](https://yiqunchen.github.io/pppower/articles/calculator-2x2.md))

## Citation

If you use `pppower` in your research, please cite:

    @software{pppower,
      title = {pppower: Prediction-Powered Inference Power Calculations},
      author = {Guo, Moran and Chen, Yiqun},
      url = {https://github.com/yiqunchen/pppower},
      year = {2025}
    }
