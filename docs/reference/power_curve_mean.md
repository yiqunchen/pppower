# Power curve for the PP mean estimator

Computes Monte Carlo power values for a grid of labeled sample sizes
using
[`simulate_power()`](https://yiqunchen.github.io/pppower/reference/simulate_power.md)
given pre-computed variance components.

## Usage

``` r
power_curve_mean(
  n_grid,
  delta,
  N,
  var_f,
  var_res,
  alpha = 0.05,
  R = 10000,
  seed = NULL
)
```

## Arguments

- n_grid:

  Integer vector of candidate labeled sample sizes.

- delta:

  Effect size \\\theta - \theta_0\\.

- N:

  Unlabeled sample size.

- var_f:

  Variance of \\f(X)\\.

- var_res:

  Variance of residuals \\Y - f(X)\\.

- alpha:

  Two-sided significance level.

- R:

  Number of Monte Carlo replicates passed to
  [`simulate_power()`](https://yiqunchen.github.io/pppower/reference/simulate_power.md).

- seed:

  Optional RNG seed for reproducibility.

## Value

Data frame with columns `n`, `power_empirical`, `power_exact`, `delta`,
`N`, `alpha`, `var_f`, and `var_res` ready for plotting power curves.

## Examples

``` r
power_curve_mean(
  n_grid = seq(50, 200, by = 25),
  delta = 0.2,
  N = 2000,
  var_f = 0.4,
  var_res = 1.0,
  R = 5000,
  seed = 123
)
#>     n power_empirical power_exact delta    N alpha var_f var_res
#> 1  50          0.2906   0.2905906   0.2 2000  0.05   0.4       1
#> 2  75          0.4050   0.4049879   0.2 2000  0.05   0.4       1
#> 3 100          0.5082   0.5081511   0.2 2000  0.05   0.4       1
#> 4 125          0.5982   0.5982060   0.2 2000  0.05   0.4       1
#> 5 150          0.6750   0.6749441   0.2 2000  0.05   0.4       1
#> 6 175          0.7392   0.7391332   0.2 2000  0.05   0.4       1
#> 7 200          0.7920   0.7920460   0.2 2000  0.05   0.4       1
```
