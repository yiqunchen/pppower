# Power curve for PPI mean with supplied population moments

Computes **analytical** and **Monte Carlo** power for the
prediction-powered mean estimator across a grid of labeled sample sizes
`n_grid`, given population moments and an effect size. This function
does **not** simulate a superpopulation; instead, you provide the
required moments: \\\sigma_f^2 = Var(f(X))\\, \\\sigma\_{\mathrm{res}}^2
= Var(Y - f(X))\\, and the effect size \\\delta = \theta - \theta_0\\.
For each `n` in `n_grid`, analytical power is computed from the normal
theory formula using \\\sqrt{\sigma_f^2/N +
\sigma\_{\mathrm{res}}^2/n}\\, and empirical power is estimated via
Monte Carlo using `simulate_ppi_vanilla_mean()`.

## Usage

``` r
power_curve_mean_dgp(
  n_grid = seq(50, 200, by = 25),
  N,
  theta0 = 0,
  var_f = 0.45,
  var_res = 0.8,
  delta = 0.2,
  family = stats::gaussian(),
  R = 1000,
  alpha = 0.05,
  seed = 1
)
```

## Arguments

- n_grid:

  Numeric vector of labeled sample sizes to evaluate.

- N:

  Unlabeled sample size.

- theta0:

  Null value for the mean estimand \\\theta_0\\.

- var_f:

  Population variance of predictions \\Var(f(X))\\.

- var_res:

  Population residual variance \\Var(Y - f(X))\\.

- delta:

  Effect size \\\theta - \theta_0\\.

- family:

  GLM family (default
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html)). Included
  for API consistency; the analytical calculation depends only on
  `var_f`, `var_res`, `N`, `n`, and `delta`.

- R:

  Number of Monte Carlo replicates for empirical power (default 1000).

- alpha:

  Two-sided test level (default 0.05).

- seed:

  RNG seed used inside the per-`n` simulations.

## Value

A data frame with one row per `n` in `n_grid`, containing:

- n:

  Labeled sample size.

- analytical:

  Analytical power from the normal theory formula.

- empirical:

  Monte Carlo power (rejection rate).

- abs_diff:

  Absolute difference \|empirical - analytical\|.

## Examples

``` r
power_curve_mean_dgp(
  n_grid = seq(50, 200, by = 25),
  N = 3000,
  theta0 = 0,
  var_f = 0.45,
  var_res = 0.80,
  delta = 0.2,
  R = 1000
)
#>     n analytical empirical     abs_diff
#> 1  50  0.3498847     0.373 0.0231153218
#> 2  75  0.4853155     0.510 0.0246845069
#> 3 100  0.6008198     0.582 0.0188198063
#> 4 125  0.6954147     0.696 0.0005853261
#> 5 150  0.7706339     0.754 0.0166338854
#> 6 175  0.8291243     0.836 0.0068757262
#> 7 200  0.8738207     0.871 0.0028206571
```
