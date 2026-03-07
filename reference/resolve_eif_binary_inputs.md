# Resolve EIF inputs for binary surrogates

Accepts either the trio `(mu0, mu1, p_hat)` or a measurement-error
parameterization via `(p, sens, spec)` and returns the quantities needed
by EIF power/sample-size routines.

## Usage

``` r
resolve_eif_binary_inputs(
  mu0 = NULL,
  mu1 = NULL,
  p_hat = NULL,
  p = NULL,
  sens = NULL,
  spec = NULL
)
```

## Arguments

- mu0, mu1:

  Conditional means \\E\[Y \mid \hat Y = 0\]\\ and \\E\[Y \mid \hat Y =
  1\]\\.

- p_hat:

  Prevalence of predicted positives \\P(\hat Y = 1)\\.

- p:

  True prevalence \\P(Y = 1)\\.

- sens, spec:

  Sensitivity and specificity of the surrogate classifier.

## Value

List with `mu0`, `mu1`, `p_hat`, and `theta` (prevalence).
