# PPI++ Paired t-Test

Performs a paired hypothesis test for the mean difference using PPI++.
This is equivalent to a one-sample test on the differences D = Y_A -
Y_B.

## Usage

``` r
ppi_paired_test(D_L, fD_L, fD_U, delta0 = 0, alpha = 0.05, lambda = NULL)
```

## Arguments

- D_L:

  Numeric vector of labeled differences (Y_A - Y_B).

- fD_L:

  Numeric vector of predicted differences for labeled samples.

- fD_U:

  Numeric vector of predicted differences for unlabeled samples.

- delta0:

  Null hypothesis value for the mean difference (default 0).

- alpha:

  Significance level (default 0.05).

- lambda:

  Optional lambda. If NULL, uses EIF-optimal.

## Value

A list with components (same as `ppi_mean_test`).

## Examples

``` r
set.seed(123)
n <- 50; N <- 500
D_L <- rnorm(n, mean = 0.3)  # True difference
fD_L <- D_L + rnorm(n, sd = 0.5)
fD_U <- rnorm(N, mean = 0.3) + rnorm(N, sd = 0.5)
ppi_paired_test(D_L, fD_L, fD_U)
#> $estimate
#> [1] 0.2717196
#> 
#> $se
#> [1] 0.06818877
#> 
#> $z_stat
#> [1] 3.984814
#> 
#> $p_value
#> [1] 6.753292e-05
#> 
#> $reject
#> [1] TRUE
#> 
#> $lambda
#> [1] 0.741807
#> 
#> $ci
#> [1] 0.1380721 0.4053671
#> 
```
