# Moments for binary outcome/prediction pairs from sensitivity/specificity

Moments for binary outcome/prediction pairs from sensitivity/specificity

## Usage

``` r
binary_moments_from_sens_spec(p, sens, spec)
```

## Arguments

- p:

  Prevalence \\P(Y = 1)\\.

- sens:

  Classifier sensitivity \\P(\hat Y = 1 \mid Y = 1)\\.

- spec:

  Classifier specificity \\P(\hat Y = 0 \mid Y = 0)\\.

## Value

List with `sigma_y2`, `sigma_f2`, `cov_y_f`, and `p_hat`.
