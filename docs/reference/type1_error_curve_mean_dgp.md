# Type I error and power curve for the PPI mean estimator

Computes empirical and analytical Type I error (or power) across a grid
of effect sizes for the prediction-powered mean estimator, **without
simulating a superpopulation**. Instead of generating data, you supply
the key population quantities directly: \\\sigma_f^2 = Var(f(X))\\ and
\\\sigma\_{\mathrm{res}}^2 = Var(Y - f(X))\\. Analytical power is
computed from the normal theory formula \\\sqrt{\sigma_f^2/N +
\sigma\_{\mathrm{res}}^2/n}\\, while empirical power is estimated via
repeated Monte Carlo simulations using `simulate_ppi_vanilla_mean()`.

## Usage

``` r
type1_error_curve_mean_dgp(
  effect_grid,
  N,
  n,
  var_f,
  var_res,
  theta = 0,
  alpha = 0.05,
  R = 2000,
  seed = 1
)
```

## Arguments

- effect_grid:

  Numeric vector of effect sizes \\\delta = \theta - \theta_0\\ over
  which to evaluate Type I error / power.

- N:

  Unlabeled sample size.

- n:

  Labeled sample size.

- var_f:

  Variance of predictions \\Var(f(X))\\.

- var_res:

  Residual variance \\Var(Y - f(X))\\.

- theta:

  True mean of \\Y\\ under the alternative. Defaults to 0. Only affects
  \\\theta_0 = \theta - \delta\\; it does not change power.

- alpha:

  Two-sided significance level (default 0.05).

- R:

  Number of Monte Carlo replicates (default 2000).

- seed:

  RNG seed for reproducibility (default 1).

## Value

A data frame with one row per value in `effect_grid`, containing:

- effect_size:

  The effect size \\\delta\\.

- type1_empirical:

  Empirical power (or Type I error when \\\delta = 0\\).

- type1_exact:

  Analytical power from the normal theory formula.

- theta:

  True mean of \\Y\\.

- theta0:

  Null value \\\theta_0\\.

- N:

  Unlabeled sample size.

- n:

  Labeled sample size.

- alpha:

  Significance level.

- var_f:

  Prediction variance.

- var_res:

  Residual variance.

- family:

  Distribution family used (`"gaussian"`).

## Examples

``` r
type1_error_curve_mean_dgp(
  effect_grid = seq(-0.4, 0.4, by = 0.05),
  N           = 4000,
  n           = 200,
  var_f       = 0.45,
  var_res     = 0.80,
  theta       = 0,
  R           = 2000
)
#>    effect_size type1_empirical type1_exact theta theta0    N   n alpha var_f
#> 1        -0.40          1.0000   0.9999905     0   0.40 4000 200  0.05  0.45
#> 2        -0.35          1.0000   0.9997654     0   0.35 4000 200  0.05  0.45
#> 3        -0.30          0.9955   0.9967173     0   0.30 4000 200  0.05  0.45
#> 4        -0.25          0.9695   0.9737153     0   0.25 4000 200  0.05  0.45
#> 5        -0.20          0.8810   0.8767233     0   0.20 4000 200  0.05  0.45
#> 6        -0.15          0.6420   0.6476942     0   0.15 4000 200  0.05  0.45
#> 7        -0.10          0.3370   0.3445730     0   0.10 4000 200  0.05  0.45
#> 8        -0.05          0.1270   0.1220192     0   0.05 4000 200  0.05  0.45
#> 9         0.00          0.0530   0.0500000     0   0.00 4000 200  0.05  0.45
#> 10        0.05          0.1280   0.1220192     0  -0.05 4000 200  0.05  0.45
#> 11        0.10          0.3595   0.3445730     0  -0.10 4000 200  0.05  0.45
#> 12        0.15          0.6465   0.6476942     0  -0.15 4000 200  0.05  0.45
#> 13        0.20          0.8715   0.8767233     0  -0.20 4000 200  0.05  0.45
#> 14        0.25          0.9775   0.9737153     0  -0.25 4000 200  0.05  0.45
#> 15        0.30          0.9985   0.9967173     0  -0.30 4000 200  0.05  0.45
#> 16        0.35          0.9995   0.9997654     0  -0.35 4000 200  0.05  0.45
#> 17        0.40          1.0000   0.9999905     0  -0.40 4000 200  0.05  0.45
#>    var_res   family
#> 1      0.8 gaussian
#> 2      0.8 gaussian
#> 3      0.8 gaussian
#> 4      0.8 gaussian
#> 5      0.8 gaussian
#> 6      0.8 gaussian
#> 7      0.8 gaussian
#> 8      0.8 gaussian
#> 9      0.8 gaussian
#> 10     0.8 gaussian
#> 11     0.8 gaussian
#> 12     0.8 gaussian
#> 13     0.8 gaussian
#> 14     0.8 gaussian
#> 15     0.8 gaussian
#> 16     0.8 gaussian
#> 17     0.8 gaussian
```
