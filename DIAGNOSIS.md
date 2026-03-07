# Package Diagnosis: pppower

## Current State Analysis

### Exported Functions (27 total - TOO MANY!)

**Power Calculations:** -
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md) -
Analytical power for PPI mean ✓ (KEEP) - `power_ppi_ols()` - Analytical
power for PPI OLS ✓ (KEEP) -
[`power_ppi_pp_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
Power for PPI++ mean ✓ (KEEP) - `power_ppi_pp_ols()` - Power for PPI++
OLS ✓ (KEEP)

**Sample Size:** - `n_required_pp()` - Required n for PPI
(mean/ols/custom) ✓ (KEEP) -
[`n_required_ppi_pp()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
Required n for PPI++ ✓ (KEEP)

**Simulation Functions:** -
[`simulate_power()`](https://yiqunchen.github.io/pppower/reference/simulate_power.md) -
Monte Carlo power (generic) ⚠️ (CONSOLIDATE) -
[`simulate_power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
MC power for mean ⚠️ (CONSOLIDATE) -
[`simulate_power_ppiplus_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
MC power for PPI++ mean ⚠️ (CONSOLIDATE) -
[`simulate_crossfit_data()`](https://yiqunchen.github.io/pppower/reference/simulate_crossfit_data.md) -
Generate synthetic data ⚠️ (MOVE TO SIMULATIONS)

**Plotting Functions:** -
[`power_curve_mean()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean.md) -
Power curve (uses simulate_power) ⚠️ (CONSOLIDATE) -
[`power_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean_dgp.md) -
Power curve (uses simulate_power_ppi_mean) ⚠️ (CONSOLIDATE) -
`power_curve_gaussian()` - Wrapper for gaussian ⚠️ (REMOVE - thin
wrapper) - `power_curve_binomial()` - Wrapper for binomial ⚠️ (REMOVE -
thin wrapper) -
[`type1_error_curve_mean()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean.md) -
Type I error curve ⚠️ (CONSOLIDATE) -
[`type1_error_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean_dgp.md) -
Type I error curve (DGP version) ⚠️ (CONSOLIDATE) -
[`plot_type1_error_curve()`](https://yiqunchen.github.io/pppower/reference/plot_type1_error_curve.md) -
Plot helper ✓ (KEEP)

**Fitting Functions:** - `ppi_ols_fit()` - Fit PPI OLS ✓ (KEEP - useful
for inference) - `ppi_ols_wald()` - Wald test ⚠️ (MAYBE KEEP - useful
for inference) - `ppi_plus_ols_fit()` - Fit PPI++ OLS ✓ (KEEP)

## Issues Identified

### 1. Too Many Exported Functions

- **27 exported functions** is excessive for a power analysis package
- The `pwr` package has ~10 main functions
- Many functions are internal helpers that shouldn’t be exported

### 2. Simulation Code Mixed with Core Code

- Simulation functions (`simulate_*`) are in `R/` alongside analytical
  functions
- Should be separated into `inst/simulations/` or `inst/scripts/`
- Makes the package API confusing

### 3. Redundant Wrapper Functions

- `power_curve_gaussian()` and `power_curve_binomial()` are thin
  wrappers
- [`power_curve_mean()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean.md)
  vs
  [`power_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/power_curve_mean_dgp.md) -
  two ways to do similar things
- [`type1_error_curve_mean()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean.md)
  vs
  [`type1_error_curve_mean_dgp()`](https://yiqunchen.github.io/pppower/reference/type1_error_curve_mean_dgp.md) -
  same issue

### 4. Inconsistent Naming

- Some functions use `_ppi_`, some use `_pp_`, some use neither
- [`simulate_power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md)
  vs
  [`simulate_power_ppiplus_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
  inconsistent

### 5. Premature Optimization

- Multiple ways to compute the same thing
- Complex internal helpers exposed to users
- Overly flexible APIs that confuse rather than help

## Proposed Refactoring

### Core API (Following `pwr` Package Pattern)

**Essential Functions to Keep:** 1.
[`power_ppi_mean()`](https://yiqunchen.github.io/pppower/reference/power_ppi_mean.md) -
Analytical power for PPI mean 2. `power_ppi_ols()` - Analytical power
for PPI OLS  
3.
[`power_ppi_pp_mean()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
Analytical power for PPI++ mean 4. `power_ppi_pp_ols()` - Analytical
power for PPI++ OLS 5. `n_required_pp()` - Sample size for PPI (unified)
6.
[`n_required_ppi_pp()`](https://yiqunchen.github.io/pppower/reference/deprecated.md) -
Sample size for PPI++ 7. `ppi_ols_fit()` - Fit PPI OLS (for inference)
8. `ppi_plus_ols_fit()` - Fit PPI++ OLS (for inference)

**Optional but Useful:** 9.
[`plot_type1_error_curve()`](https://yiqunchen.github.io/pppower/reference/plot_type1_error_curve.md) -
Plotting helper

**Functions to Consolidate/Remove:** - Consolidate `power_curve_*`
functions into one flexible function - Consolidate `type1_error_curve_*`
functions into one flexible function - Consolidate `simulate_power_*`
functions (or move to simulations folder) - Remove thin wrappers
(`power_curve_gaussian`, `power_curve_binomial`)

### Proposed Structure

    R/
      ├── power-mean.R          # power_ppi_mean(), power_ppi_pp_mean()
      ├── power-ols.R           # power_ppi_ols(), power_ppi_pp_ols()
      ├── sample-size.R         # n_required_pp(), n_required_ppi_pp()
      ├── fit-ols.R             # ppi_ols_fit(), ppi_plus_ols_fit()
      ├── plotting.R            # plot_type1_error_curve(), power_curve(), type1_curve()
      ├── utils.R               # Internal helpers (resolve_ppi_variances, etc.)
      └── pppower-package.R     # Package documentation

    inst/
      └── simulations/          # Simulation scripts (not exported)
          ├── simulate-power-mean.R
          ├── simulate-power-ols.R
          ├── simulate-crossfit-data.R
          └── README.md          # How to use simulation scripts

### Key Principles

1.  **Simplicity**: One clear way to do each task
2.  **Separation**: Core power calculations separate from simulation
    code
3.  **Consistency**: Follow `pwr` package naming and structure patterns
4.  **Minimal API**: Only export what users need
5.  **Internal Helpers**: Keep utility functions internal unless truly
    needed

## Next Steps

1.  Create `inst/simulations/` folder
2.  Move simulation code to simulations folder
3.  Consolidate similar functions
4.  Update NAMESPACE to only export essential functions
5.  Update documentation
6.  Update tests
