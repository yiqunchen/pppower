# Generate One Labeled + Unlabeled Sample

Generate One Labeled + Unlabeled Sample

## Usage

``` r
simulate_one_draw(
  n_labeled,
  n_unlabeled,
  X_sampler_L,
  X_sampler_U = NULL,
  f_generator,
  eps_sampler = function(n) rnorm(n, 0, 1),
  delta = 0
)
```

## Arguments

- n_labeled:

  Number of labeled samples

- n_unlabeled:

  Number of unlabeled samples

- X_sampler_L:

  Function sampling labeled covariates

- X_sampler_U:

  Function sampling unlabeled covariates (defaults to X_sampler_L)

- f_generator:

  Function generating f(X)

- eps_sampler:

  Error sampler, default N(0,1)

- delta:

  Mean shift applied to f

## Value

A list with labeled and unlabeled data and true mean.
