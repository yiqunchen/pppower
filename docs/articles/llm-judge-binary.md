# Real Data: LLM-as-a-Judge (Binary Surrogates)

## Setup

This vignette mirrors the GPT-4o-mini judge example from
`../llm-as-a-judge-debias`. We treat the LLM decision
$`\hat{Y}\in\{0,1\}`$ as a binary surrogate for the human label $`Y`$,
and we plan calibration size and power for the EIF estimator.

Assume we have estimates:

- Prevalence $`p = P(Y=1)`$ from pilot human labels
- Sensitivity $`q_1 = P(\hat{Y}=1 \mid Y=1)`$
- Specificity $`q_0 = P(\hat{Y}=0 \mid Y=0)`$

``` r

p    <- 0.236
sens <- 0.25
spec <- 0.875
N    <- 5000   # unlabeled/test comparisons
```

## Required calibration size for a win-rate shift

Suppose we want $`80\%`$ power to detect a $`\Delta = 0.05`$ shift in
win rate.

``` r

n_req <- power_eif_binary(
  delta = 0.05,
  N     = N,
  m_cal = NULL,
  power = 0.80,
  p     = p,
  sens  = sens,
  spec  = spec
)
n_req
#> [1] 556
```

With the above numbers, $`\approx 550`$ labeled comparisons are needed
for $`80\%`$ power; using only 100 or 250 labeled comparisons delivers
much lower power:

``` r

data.frame(
  m_cal = c(100, 250, n_req),
  power = sapply(c(100, 250, n_req), function(m) {
    power_eif_binary(
      delta = 0.05,
      N     = N,
      m_cal = m,
      p     = p,
      sens  = sens,
      spec  = spec
    )
  })
)
#>   m_cal     power
#> 1   100 0.2215246
#> 2   250 0.4686672
#> 3   556 0.8005491
```

## Monte Carlo check

``` r

sim <- simulate_eif_binary(
  R      = 500,
  delta  = 0.05,
  N      = N,
  m_cal  = n_req,
  p      = p,
  sens   = sens,
  spec   = spec,
  seed   = 2025
)
c(empirical = sim$empirical_power, analytical = sim$analytical_power)
#>  empirical analytical 
#>  0.8300000  0.8005491
```

## Using sensitivity/specificity vs. conditional means

If you already have $`\hat\mu_0 = E[Y \mid \hat{Y}=0]`$,
$`\hat\mu_1 = E[Y \mid \hat{Y}=1]`$, and $`p_{\hat{Y}}`$, you can pass
them directly:

``` r

power_eif_binary(
  delta = 0.05,
  N     = N,
  m_cal = n_req,
  mu0   = 0.05,
  mu1   = 0.55,
  p_hat = 0.16
)
#> [1] 0.9832705
#> attr(,"theta")
#> [1] 0.13
#> attr(,"p_hat")
#> [1] 0.16
#> attr(,"var_unlabeled")
#> [1] 6.72e-06
#> attr(,"var_cal")
#> [1] 0.0001429856
```

Either path yields the same variance pieces; choose whichever matches
your data summary.
