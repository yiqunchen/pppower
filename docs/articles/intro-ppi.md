# Quickstart: Prediction-Powered Power & Sample Size

## What this vignette covers

This is the **basic-use** guide for `pppower`. It shows the minimum
calls you need for power, sample size, paired designs, and 2x2
contingency-table planning. See the “Detailed Dive” and “2x2 Calculator”
tabs on the pkgdown site for derivations, regression contrasts, and
odds-ratio / relative-risk setups.

``` r

library(pppower)
```

### What it does (analogy to `pwr`)

- Drop-in power / sample-size calculators that leverage predictions on
  unlabeled data (`PPI++`).
- Same spirit as
  [`pwr::pwr.t.test`](https://rdrr.io/pkg/pwr/man/pwr.t.test.html), but
  variance comes from predictions + residuals.
- Works for means, paired differences, regression contrasts, and 2x2
  contingency tables.

### 1. Power for a prediction-powered mean test

``` r

power_ppi_mean(
  delta    = 0.2,   # effect size
  N        = 5000,  # unlabeled size
  n        = 200,   # labeled size
  var_f    = 0.6,   # Var(f(X))
  var_res  = 0.4,   # Var(Y - f(X))
  cov_y_f  = 0.6    # Cov(Y, f)
)
#> [1] 0.9915412
```

### 2. Required labeled size for target power

``` r

power_ppi_mean(
  delta    = 0.2,
  N        = 5000,
  n        = NULL,
  power    = 0.80,
  sigma_y2 = 1.0,
  sigma_f2 = 0.6,
  cov_y_f  = 0.6
)
#> [1] 81
```

You can also supply prediction metrics (MSE, $`R^2`$,
sensitivity/specificity) instead of raw moments; see the “Detailed Dive”
vignette.

### 3. Paired designs in one line

``` r

power_ppi_paired(
  delta    = 0.3,   # mean difference
  N        = 1000,  # unlabeled pairs
  n        = 50,    # labeled pairs
  sigma_D2 = 1.0,
  rho_D    = 0.7
)
#> [1] 0.8276133
```

### 4. 2x2 designs in one line

``` r

power_ppi_2x2(
  p_exp          = 0.40,
  p_ctrl         = 0.20,
  N              = 1000,
  power          = 0.80,
  sens           = 0.85,
  spec           = 0.85,
  effect_measure = "OR"
)
#> [1] 108
```

### Where next?

- **Detailed Dive:** formulas, metrics-to-moments mapping, regression
  contrasts.
- **2x2 Calculator:** odds-ratio and relative-risk planning from
  contingency tables.
