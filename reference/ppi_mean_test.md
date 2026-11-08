# PPI++ One-Sample Mean Test

Performs a one-sample hypothesis test for the population mean using the
PPI++ (Prediction-Powered Inference Plus) estimator with EIF-optimal
lambda.

## Usage

``` r
ppi_mean_test(
  Y_L,
  f_L,
  f_U,
  theta0 = 0,
  alpha = 0.05,
  lambda = NULL,
  se_type = c("standard", "corrected")
)
```

## Arguments

- Y_L:

  Numeric vector of labeled outcomes.

- f_L:

  Numeric vector of predictions for labeled samples.

- f_U:

  Numeric vector of predictions for unlabeled samples.

- theta0:

  Null hypothesis value (default 0).

- alpha:

  Significance level (default 0.05).

- lambda:

  Optional user-specified lambda. If NULL, uses EIF-optimal lambda.

- se_type:

  `"standard"` (default) uses the usual plug-in standard error and a
  normal reference. `"corrected"` adds the second-order term
  \\c\\\hat\sigma\_\varepsilon^2\hat\kappa/\\n(n-3)\\\\ for the
  estimated \\\lambda\\, with \\\hat\kappa\\ bounded below by 1, and
  uses a \\t\_{n-2}\\ reference. It is meant for small `n`; it is
  ignored when `lambda` is supplied, because a fixed \\\lambda\\ needs
  no correction.

## Value

A list with components:

- estimate:

  The PPI++ point estimate of the mean.

- se:

  Standard error of the estimate.

- z_stat:

  Test statistic (compared with a \\t\_{n-2}\\ distribution when
  `se_type = "corrected"`).

- df:

  Degrees of freedom of the reference distribution (`Inf` for normal).

- p_value:

  Two-sided p-value.

- reject:

  Logical indicating whether null is rejected at level alpha.

- lambda:

  The lambda value used.

- ci:

  Confidence interval at level 1-alpha.

## Examples

``` r
set.seed(123)
n <- 100; N <- 1000
Y_L <- rnorm(n, mean = 0.5)
f_L <- Y_L + rnorm(n, sd = 0.5)
f_U <- rnorm(N, mean = 0.5) + rnorm(N, sd = 0.5)
ppi_mean_test(Y_L, f_L, f_U, theta0 = 0)
#> $estimate
#> [1] 0.5867406
#> 
#> $se
#> [1] 0.04983671
#> 
#> $z_stat
#> [1] 11.77326
#> 
#> $df
#> [1] Inf
#> 
#> $p_value
#> [1] 5.361031e-32
#> 
#> $reject
#> [1] TRUE
#> 
#> $lambda
#> [1] 0.720832
#> 
#> $ci
#> [1] 0.4890625 0.6844188
#> 
```
