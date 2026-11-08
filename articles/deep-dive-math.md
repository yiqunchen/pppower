# Detailed Dive: Variance Formulas and lambda-star

## Mean estimator and `PPI++`

For labeled size $`n`$ and unlabeled size $`N`$:

*Estimator*
``` math
\hat\theta_\lambda = \bar{Y}_L + \lambda(\bar{f}_U - \bar{f}_L).
```

*Variance*
``` math
\operatorname{Var}(\hat\theta_\lambda) =
\frac{\sigma_Y^2}{n}
 + \lambda^2 \sigma_f^2 \left(\frac{1}{N} + \frac{1}{n}\right)
 - \frac{2\lambda\,\operatorname{Cov}(Y,f)}{n}.
```

Minimizing yields
``` math
\lambda^* = \frac{\operatorname{Cov}(Y,f)}{(1 + n/N)\sigma_f^2}, \qquad
\operatorname{Var}(\hat\theta_{\lambda^*}) =
\frac{\sigma_Y^2}{n} - \frac{\operatorname{Cov}(Y,f)^2}{\sigma_f^2}\frac{N}{n(n+N)}.
```

``` r

lambda_opt <- function(n, N, cov_y_f, var_f) cov_y_f / ((1 + n / N) * var_f)
lambda_opt(n = 200, N = 5000, cov_y_f = 0.6, var_f = 0.6)
#> [1] 0.9615385
```

## Power and required $`n`$

Given effect size $`\Delta`$, target power $`1-\beta`$, and level
$`\alpha`$,
``` math
\frac{\lvert\Delta\rvert}{\sqrt{\operatorname{Var}(\hat\theta_{\lambda^*})}}
\ge z_{1-\alpha/2} + z_{1-\beta}.
```
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

For sensitivity/specificity $`(\text{sens}, \text{spec})`$ and
prevalence $`p`$,
``` math
p_{\hat{Y}} = \text{sens}\,p + (1-\text{spec})(1-p), \qquad
\operatorname{Cov}(Y,f) = p(1-p)(\text{sens} + \text{spec} - 1).
```
The helper
[`binary_moments_from_sens_spec()`](https://yiqunchen.github.io/pppower/reference/binary_moments_from_sens_spec.md)
returns $`\operatorname{Var}(Y)`$, $`\operatorname{Var}(f)`$, and
$`\operatorname{Cov}(Y,f)`$:

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
