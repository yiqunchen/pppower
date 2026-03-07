# Detailed Dive: Variance Formulas and lambda-star

## Mean estimator and PPI++

For labeled size $n$ and unlabeled size $N$:

*Estimator*
$${\widehat{\theta}}_{\lambda} = {\bar{Y}}_{L} + \lambda\left( {\bar{f}}_{U} - {\bar{f}}_{L} \right).$$

*Variance*
$$\operatorname{Var}\left( {\widehat{\theta}}_{\lambda} \right) = \frac{\sigma_{Y}^{2}}{n} + \lambda^{2}\sigma_{f}^{2}\left( \frac{1}{N} + \frac{1}{n} \right) - \frac{2\lambda\,\operatorname{Cov}(Y,f)}{n}.$$

Minimizing yields \$\$ \lambda^\\ = \frac{\operatorname{Cov}(Y,f)}{(1 +
n/N)\sigma_f^2}, \qquad \operatorname{Var}(\hat\theta\_{\lambda^\\}) =
\frac{\sigma_Y^2}{n} -
\frac{\operatorname{Cov}(Y,f)^2}{\sigma_f^2}\frac{N}{n(n+N)}. \$\$

``` r
lambda_opt <- function(n, N, cov_y_f, var_f) cov_y_f / ((1 + n / N) * var_f)
lambda_opt(n = 200, N = 5000, cov_y_f = 0.6, var_f = 0.6)
#> [1] 0.9615385
```

## Power and required $n$

Given effect size $\Delta$, target power $1 - \beta$, and level
$\alpha$, \$\$
\frac{\|\Delta\|}{\sqrt{\operatorname{Var}(\hat\theta\_{\lambda^\\})}}
\ge z\_{1-\alpha/2} + z\_{1-\beta}. \$\$
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md)
solves this directly; use it instead of hand algebra:

``` r
power_ppi_mean(
  delta      = 0.2,
  N          = 5000,
  n          = NULL,
  power      = 0.9,
  sigma_y2   = 1.0,
  sigma_f2   = 0.6,
  cov_y_f    = 0.6,
  lambda_mode = "oracle"
)
#> [1] 109
```

## Metrics-to-moments mapping (binary)

For sensitivity/specificity $\left( \text{sens},\text{spec} \right)$ and
prevalence $p$,
$$p_{\widehat{Y}} = \text{sens}\, p + \left( 1 - \text{spec} \right)(1 - p),\qquad\operatorname{Cov}(Y,f) = p(1 - p)\left( \text{sens} + \text{spec} - 1 \right).$$
The helper
[`binary_moments_from_sens_spec()`](https://yiqunchen.github.io/pppower/reference/binary_moments_from_sens_spec.md)
returns $\operatorname{Var}(Y)$, $\operatorname{Var}(f)$, and
$\operatorname{Cov}(Y,f)$:

``` r
binary_moments_from_sens_spec(p = 0.3, sens = 0.85, spec = 0.85)
#> $sigma_y2
#> [1] 0.21
#> 
#> $sigma_f2
#> [1] 0.2304
#> 
#> $cov_y_f
#> [1] 0.147
#> 
#> $p_hat
#> [1] 0.36
```

Use these outputs directly in
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md).
