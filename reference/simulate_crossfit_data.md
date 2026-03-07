# Simulate cross-fitted predictive data

Generates a synthetic regression/classification dataset with covariates,
outcomes, and out-of-fold predictions from a K-fold cross-fitted GLM.

## Usage

``` r
simulate_crossfit_data(
  n = 20000,
  p = 5,
  family = stats::binomial(),
  K = 5,
  seed = 1
)
```

## Arguments

- n:

  Number of observations to simulate.

- p:

  Number of covariates (columns in the design matrix).

- family:

  GLM family (defaults to
  [`stats::binomial()`](https://rdrr.io/r/stats/family.html);
  [`stats::gaussian()`](https://rdrr.io/r/stats/family.html) gives a
  Gaussian response).

- K:

  Number of folds used for cross-fitting.

- seed:

  RNG seed for reproducibility.

## Value

A data frame containing the outcome `y`, covariate columns
`x1, ..., xp`, and `fhat_cf`, the cross-fitted predictions. Attributes
`family` and `K` record the generating family and fold count.
