# Type I error curve for the PP mean estimator

Computes empirical and analytical Type I error estimates across a grid
of effect sizes using Monte Carlo via
[`simulate_power()`](https://yiqunchen.github.io/pppower/reference/simulate_power.md).
When the null is true (`effect_size = 0`), the curve reports the Type I
error; for other effect sizes the values coincide with the rejection
probability (i.e., power).

## Usage

``` r
type1_error_curve_mean(
  effect_grid,
  N,
  n,
  var_f,
  var_res,
  alpha = 0.05,
  R = 10000,
  seed = NULL
)
```

## Arguments

- effect_grid:

  Numeric vector of effect sizes \\\theta - \theta_0\\ to evaluate.

- N:

  Unlabeled sample size.

- n:

  Labeled sample size.

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

Data frame with columns `effect_size`, `type1_empirical`, `type1_exact`,
`N`, `n`, `alpha`, `var_f`, and `var_res`.

## Examples

``` r
type1_error_curve_mean(
  effect_grid = seq(-0.4, 0.4, by = 0.05),
  N = 4000,
  n = 200,
  var_f = 0.4,
  var_res = 1.1,
  R = 2000
)
#>    effect_size type1_empirical type1_exact    N   n alpha var_f var_res
#> 1        -0.40          0.9995   0.9996444 4000 200  0.05   0.4     1.1
#> 2        -0.35          0.9965   0.9967072 4000 200  0.05   0.4     1.1
#> 3        -0.30          0.9800   0.9797667 4000 200  0.05   0.4     1.1
#> 4        -0.25          0.9165   0.9163301 4000 200  0.05   0.4     1.1
#> 5        -0.20          0.7620   0.7619701 4000 200  0.05   0.4     1.1
#> 6        -0.15          0.5175   0.5177820 4000 200  0.05   0.4     1.1
#> 7        -0.10          0.2670   0.2669161 4000 200  0.05   0.4     1.1
#> 8        -0.05          0.1025   0.1025043 4000 200  0.05   0.4     1.1
#> 9         0.00          0.0500   0.0500000 4000 200  0.05   0.4     1.1
#> 10        0.05          0.1025   0.1025043 4000 200  0.05   0.4     1.1
#> 11        0.10          0.2670   0.2669161 4000 200  0.05   0.4     1.1
#> 12        0.15          0.5175   0.5177820 4000 200  0.05   0.4     1.1
#> 13        0.20          0.7620   0.7619701 4000 200  0.05   0.4     1.1
#> 14        0.25          0.9165   0.9163301 4000 200  0.05   0.4     1.1
#> 15        0.30          0.9800   0.9797667 4000 200  0.05   0.4     1.1
#> 16        0.35          0.9965   0.9967072 4000 200  0.05   0.4     1.1
#> 17        0.40          0.9995   0.9996444 4000 200  0.05   0.4     1.1
```
