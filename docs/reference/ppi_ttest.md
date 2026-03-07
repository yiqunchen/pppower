# PPI++ Two-Sample t-Test

Performs a two-sample hypothesis test for the difference in means
between two independent groups using PPI++ estimators.

## Usage

``` r
ppi_ttest(
  Y_A,
  f_A_L,
  f_A_U,
  Y_B,
  f_B_L,
  f_B_U,
  delta0 = 0,
  alpha = 0.05,
  lambda_A = NULL,
  lambda_B = NULL
)
```

## Arguments

- Y_A:

  Numeric vector of labeled outcomes for group A.

- f_A_L:

  Numeric vector of predictions for labeled samples in group A.

- f_A_U:

  Numeric vector of predictions for unlabeled samples in group A.

- Y_B:

  Numeric vector of labeled outcomes for group B.

- f_B_L:

  Numeric vector of predictions for labeled samples in group B.

- f_B_U:

  Numeric vector of predictions for unlabeled samples in group B.

- delta0:

  Null hypothesis value for the difference (default 0).

- alpha:

  Significance level (default 0.05).

- lambda_A:

  Optional lambda for group A. If NULL, uses EIF-optimal.

- lambda_B:

  Optional lambda for group B. If NULL, uses EIF-optimal.

## Value

A list with components:

- estimate:

  The PPI++ estimate of the difference (mean_A - mean_B).

- se:

  Standard error of the difference estimate.

- z_stat:

  Z test statistic.

- p_value:

  Two-sided p-value.

- reject:

  Logical indicating whether null is rejected at level alpha.

- lambda_A:

  Lambda used for group A.

- lambda_B:

  Lambda used for group B.

- ci:

  Confidence interval at level 1-alpha.

## Examples

``` r
set.seed(123)
n <- 50; N <- 500
# Group A (treatment)
Y_A <- rnorm(n, mean = 0.3)
f_A_L <- Y_A + rnorm(n, sd = 0.5)
f_A_U <- rnorm(N, mean = 0.3) + rnorm(N, sd = 0.5)
# Group B (control)
Y_B <- rnorm(n, mean = 0)
f_B_L <- Y_B + rnorm(n, sd = 0.5)
f_B_U <- rnorm(N, mean = 0) + rnorm(N, sd = 0.5)

ppi_ttest(Y_A, f_A_L, f_A_U, Y_B, f_B_L, f_B_U)
#> $estimate
#> [1] 0.2489567
#> 
#> $se
#> [1] 0.1016156
#> 
#> $z_stat
#> [1] 2.449985
#> 
#> $p_value
#> [1] 0.01428622
#> 
#> $reject
#> [1] TRUE
#> 
#> $lambda_A
#> [1] 0.741807
#> 
#> $lambda_B
#> [1] 0.766848
#> 
#> $ci
#> [1] 0.04979377 0.44811963
#> 
```
