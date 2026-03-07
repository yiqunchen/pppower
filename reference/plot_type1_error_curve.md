# Plot Type I error curves

Creates a simple line plot of empirical and/or analytical Type I error
(or rejection probabilities) against effect size.

## Usage

``` r
plot_type1_error_curve(
  curve_df,
  empirical = TRUE,
  exact = TRUE,
  add_reference = TRUE,
  empirical_col = "type1_empirical",
  exact_col = "type1_exact",
  legend_pos = "topright",
  ...
)
```

## Arguments

- curve_df:

  Data frame returned by
  [`type1_error_curve_mean()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean.md)
  or
  [`type1_error_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean_dgp.md).

- empirical:

  Logical; include the Monte Carlo estimate (`type1_empirical`).

- exact:

  Logical; include the analytical estimate (`type1_exact`).

- add_reference:

  Logical; add a horizontal line at the nominal level `alpha` when
  available.

- empirical_col, exact_col:

  Column names to use for empirical and exact curves. Override only if
  you have renamed the defaults.

- legend_pos:

  Character or numeric legend position passed to
  [`graphics::legend()`](https://rdrr.io/r/graphics/legend.html),
  ignored when fewer than two curves are drawn.

- ...:

  Additional arguments forwarded to the initial
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html)
  call (e.g., `main`, `xlab`, `ylab`).

## Value

The input data frame, invisibly.
