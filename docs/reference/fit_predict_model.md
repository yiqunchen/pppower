# Fit predictive model and return predictions for Mean Estimation for PPI / PPI++

Internal unified interface for fitting predictive models used in the PPI
and PPI++ estimators. Supports correctly specified, misspecified, and
incorrectly specified linear models, as well as random forests. Returns
fitted model object and predictions on the unlabeled covariates.

## Usage

``` r
fit_predict_model(
  model_type,
  X_L,
  y_L,
  X_U,
  mtry = NULL,
  rf_engine = c("ranger", "randomForest"),
  rf_trees = 200,
  rf_min_node_size = 5,
  rf_num_threads = NULL,
  rf_seed = NULL
)
```

## Arguments

- model_type:

  Type of model ("glm_correct", "glm_mis", "glm_wrong", "rf")

- X_L:

  Labeled covariates

- y_L:

  Labeled outcomes

- X_U:

  Unlabeled covariates

- mtry:

  RF mtry value

- rf_engine:

  Random forest engine ("ranger" or "randomForest")

- rf_trees:

  Number of trees

- rf_min_node_size:

  Minimum node size for ranger

- rf_num_threads:

  Threads (ranger)

- rf_seed:

  Random seed

## Value

List containing model fit object and predictions on X_U.
