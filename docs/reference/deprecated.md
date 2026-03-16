# Deprecated functions in pppower

These functions have been renamed or merged into unified interfaces
following the pwr-style convention of "leave one NULL" (compute
whichever of power or sample size is left `NULL`).

- `power_ppi_pp_mean`:

  Use
  [`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md)
  instead.

- `power_ppi_pp_paired`:

  Use
  [`power_ppi_paired()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired.md)
  instead.

- `power_ppi_pp_paired_binary`:

  Use
  [`power_ppi_paired_binary()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired_binary.md)
  instead.

- `power_ppi_pp_ttest`:

  Use
  [`power_ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest.md)
  instead.

- `power_ppi_pp_ttest_binary`:

  Use
  [`power_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest_binary.md)
  instead.

- `n_required_ppi_pp`:

  Use
  [`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md)
  or
  [`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)
  with `n = NULL`.

- `n_required_ppi_pp_paired`:

  Use
  [`power_ppi_paired()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired.md)
  with `n = NULL`.

- `simulate_power_ppiplus_mean`:

  Use
  [`simulate_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_mean.md)
  instead.

- `simulate_power_ppi_pp_ttest_binary`:

  Use
  [`simulate_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_ttest_binary.md)
  instead.

- `simulate_power_ppi_mean`:

  Use `simulate_ppi_vanilla_mean()` or
  [`simulate_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_mean.md).

## Usage

``` r
simulate_power_ppi_mean(...)

power_ppi_pp_mean(delta, N, n, ...)

power_ppi_pp_paired(...)

power_ppi_pp_paired_binary(...)

n_required_ppi_pp_paired(
  delta,
  N,
  power = 0.8,
  alpha = 0.05,
  sigma_D2,
  rho_D,
  sigma_fD2 = NULL
)

power_ppi_pp_ttest(...)

power_ppi_pp_ttest_binary(...)

simulate_power_ppiplus_mean(...)

simulate_power_ppi_pp_ttest_binary(...)

n_required_ppi_pp(
  delta,
  N,
  alpha = 0.05,
  power = 0.8,
  type = c("mean", "regression"),
  lambda_mode = c("vanilla", "oracle", "user"),
  lambda_user = NULL,
  sigma_y2 = NULL,
  sigma_f2 = NULL,
  cov_y_f = NULL,
  var_f = NULL,
  var_res = NULL,
  metrics = NULL,
  metric_type = NULL,
  m_labeled = NULL,
  correction = TRUE,
  c = NULL,
  H_L = NULL,
  H_U = NULL,
  Sigma_YY = NULL,
  Sigma_ff_l = NULL,
  Sigma_ff_u = NULL,
  Sigma_Yf = NULL,
  warn_smallN = TRUE,
  smallN_threshold = 500,
  mode = c("error", "cap")
)
```
