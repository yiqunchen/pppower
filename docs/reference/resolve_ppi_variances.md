# Resolve PP variance components from metrics

Resolve PP variance components from metrics

## Usage

``` r
resolve_ppi_variances(
  var_f = NULL,
  var_res = NULL,
  metrics = NULL,
  metric_type = NULL,
  m_labeled = NULL,
  correction = TRUE
)
```

## Arguments

- var_f:

  Optional numeric scalar supplying \\\operatorname{Var}(f)\\ directly.

- var_res:

  Optional numeric scalar supplying \\\operatorname{Var}(Y-f)\\
  directly.

- metrics:

  Optional named list of predictive-performance summaries. The required
  fields depend on `metric_type`.

- metric_type:

  Character string identifying the metric bundle. Supported values are:
  `"continuous"` (regression-style metrics), `"prob"` (binary
  probabilistic metrics such as the Brier score), and `"classification"`
  (binary classification metrics such as a confusion matrix or
  precision/recall).

- m_labeled:

  Labeled sample size associated with `metrics`; defaults to
  `metrics$m_obs` when present.

- correction:

  Logical; apply the finite-sample adjustment when deriving moments from
  metrics (default `TRUE`).

## Value

List with numeric elements `var_f` and `var_res`.
