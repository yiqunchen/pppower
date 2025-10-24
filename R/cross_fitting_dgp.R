crossfit_glm <- function(X, y, K = 2, family = stats::binomial(), seed = 1) {
  n <- nrow(X)
  folds <- kfold_split(n, K, seed)
  fhat <- numeric(n)

  fam_name <- family$family  # "binomial" or "gaussian" currently

  for (k in seq_along(folds)) {
    te <- folds[[k]]
    tr <- setdiff(seq_len(n), te)

    df_tr <- data.frame(y = y, X)[tr, , drop = FALSE]
    df_te <- data.frame(X)[te, , drop = FALSE]

    # Fit on training folds only, predict on the holdout fold
    fit <- stats::glm(y ~ .,
                      data   = df_tr,
                      family = family)
    pred <- stats::predict(fit, newdata = df_te, type = "response")

    # clamp only for binomial (probabilities); not for gaussian
    if (identical(fam_name, "binomial")) {
      pred <- pmin(pmax(pred, 0), 1)
    }
    fhat[te] <- pred
  }

  fhat
}
