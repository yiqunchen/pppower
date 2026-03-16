# Package index

## Core power & sample size

- [`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md)
  : Power or required sample size for the PPI++ mean estimator
- [`power_ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest.md)
  : Power for PPI++ two-sample t-test (difference in means)
- [`power_ppi_paired()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired.md)
  : Power for PPI++ paired t-test (within-subject differences)
- [`power_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/power_ppi_ttest_binary.md)
  : Power for PPI++ two-sample test with binary outcomes (sens/spec
  input)
- [`power_ppi_paired_binary()`](https://yiqunchen.github.io/pppower/reference/power_ppi_paired_binary.md)
  : Power for PPI++ paired test with binary outcomes
- [`power_ppi_regression()`](https://yiqunchen.github.io/pppower/reference/power_ppi_regression.md)
  : Power or required sample size for PPI++ regression contrast
- [`power_ppi_2x2()`](https://yiqunchen.github.io/pppower/reference/power_ppi_2x2.md)
  : Power or required sample size for a 2x2 table with PPI++ surrogates

## Simulation helpers

- [`simulate_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_mean.md)
  : Monte Carlo power for the PPI++ mean estimator
- [`simulate_ppi_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/simulate_ppi_ttest_binary.md)
  : Monte Carlo power for binary two-sample PPI++ test
- [`simulate_crossfit_data()`](https://yiqunchen.github.io/pppower/reference/simulate_crossfit_data.md)
  : Simulate cross-fitted predictive data

## Curves & plotting

- [`power_curve_mean()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean.md)
  : Power curve for the PP mean estimator
- [`power_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean_dgp.md)
  : Power curve for PPI mean with supplied population moments
- [`type1_error_curve_mean()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean.md)
  : Type I error curve for the PP mean estimator
- [`type1_error_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean_dgp.md)
  : Type I error and power curve for the PPI mean estimator
- [`plot_type1_error_curve()`](https://yiqunchen.github.io/pppower/reference/plot_type1_error_curve.md)
  : Plot Type I error curves

## Utilities

- [`binary_moments_from_sens_spec()`](https://yiqunchen.github.io/pppower/reference/binary_moments_from_sens_spec.md)
  : Moments for binary outcome/prediction pairs from
  sensitivity/specificity
- [`resolve_ppi_variances()`](https://yiqunchen.github.io/pppower/reference/resolve_ppi_variances.md)
  : Resolve PP variance components from metrics
- [`compute_ppi_blocks()`](https://yiqunchen.github.io/pppower/reference/compute_ppi_blocks.md)
  : Construct Hessian/Fisher and Covariance Blocks for PPI/PPI++
  Regression
- [`theme_pppower()`](https://yiqunchen.github.io/pppower/reference/theme_pppower.md)
  : Publication-ready ggplot2 theme for pppower
- [`save_pppower_figure()`](https://yiqunchen.github.io/pppower/reference/save_pppower_figure.md)
  : Save a ggplot figure in PDF and/or PNG formats
- [`ppi_mean_test()`](https://yiqunchen.github.io/pppower/reference/ppi_mean_test.md)
  : PPI++ One-Sample Mean Test
- [`ppi_paired_test()`](https://yiqunchen.github.io/pppower/reference/ppi_paired_test.md)
  : PPI++ Paired t-Test
- [`ppi_ttest()`](https://yiqunchen.github.io/pppower/reference/ppi_ttest.md)
  : PPI++ Two-Sample t-Test

## Deprecated

Soft-deprecated shims; use the renamed functions above.

- [`simulate_power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`power_ppi_pp_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`power_ppi_pp_paired()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`power_ppi_pp_paired_binary()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`n_required_ppi_pp_paired()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`power_ppi_pp_ttest()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`power_ppi_pp_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`simulate_power_ppiplus_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`simulate_power_ppi_pp_ttest_binary()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  [`n_required_ppi_pp()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  : Deprecated functions in pppower
- [`simulate_power()`](https://yiqunchen.github.io/pppower/reference/simulate_power.md)
  **\[soft-deprecated\]** : Monte Carlo vs. analytical power for PPI /
  PPI++ mean estimation
