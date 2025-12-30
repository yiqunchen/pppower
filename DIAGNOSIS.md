# Package Diagnosis: pppower

## Current State Analysis

### Exported Functions (27 total - TOO MANY!)

**Power Calculations:**
- `power_ppi_mean()` - Analytical power for PPI mean ✓ (KEEP)
- `power_ppi_ols()` - Analytical power for PPI OLS ✓ (KEEP)
- `power_ppi_pp_mean()` - Power for PPI++ mean ✓ (KEEP)
- `power_ppi_pp_ols()` - Power for PPI++ OLS ✓ (KEEP)

**Sample Size:**
- `n_required_pp()` - Required n for PPI (mean/ols/custom) ✓ (KEEP)
- `n_required_ppi_pp()` - Required n for PPI++ ✓ (KEEP)

**Simulation Functions:**
- `simulate_power()` - Monte Carlo power (generic) ⚠️ (CONSOLIDATE)
- `simulate_power_ppi_mean()` - MC power for mean ⚠️ (CONSOLIDATE)
- `simulate_power_ppiplus_mean()` - MC power for PPI++ mean ⚠️ (CONSOLIDATE)
- `simulate_crossfit_data()` - Generate synthetic data ⚠️ (MOVE TO SIMULATIONS)

**Plotting Functions:**
- `power_curve_mean()` - Power curve (uses simulate_power) ⚠️ (CONSOLIDATE)
- `power_curve_mean_dgp()` - Power curve (uses simulate_power_ppi_mean) ⚠️ (CONSOLIDATE)
- `power_curve_gaussian()` - Wrapper for gaussian ⚠️ (REMOVE - thin wrapper)
- `power_curve_binomial()` - Wrapper for binomial ⚠️ (REMOVE - thin wrapper)
- `type1_error_curve_mean()` - Type I error curve ⚠️ (CONSOLIDATE)
- `type1_error_curve_mean_dgp()` - Type I error curve (DGP version) ⚠️ (CONSOLIDATE)
- `plot_type1_error_curve()` - Plot helper ✓ (KEEP)

**Fitting Functions:**
- `ppi_ols_fit()` - Fit PPI OLS ✓ (KEEP - useful for inference)
- `ppi_ols_wald()` - Wald test ⚠️ (MAYBE KEEP - useful for inference)
- `ppi_plus_ols_fit()` - Fit PPI++ OLS ✓ (KEEP)

## Issues Identified

### 1. Too Many Exported Functions
- **27 exported functions** is excessive for a power analysis package
- The `pwr` package has ~10 main functions
- Many functions are internal helpers that shouldn't be exported

### 2. Simulation Code Mixed with Core Code
- Simulation functions (`simulate_*`) are in `R/` alongside analytical functions
- Should be separated into `inst/simulations/` or `inst/scripts/`
- Makes the package API confusing

### 3. Redundant Wrapper Functions
- `power_curve_gaussian()` and `power_curve_binomial()` are thin wrappers
- `power_curve_mean()` vs `power_curve_mean_dgp()` - two ways to do similar things
- `type1_error_curve_mean()` vs `type1_error_curve_mean_dgp()` - same issue

### 4. Inconsistent Naming
- Some functions use `_ppi_`, some use `_pp_`, some use neither
- `simulate_power_ppi_mean()` vs `simulate_power_ppiplus_mean()` - inconsistent

### 5. Premature Optimization
- Multiple ways to compute the same thing
- Complex internal helpers exposed to users
- Overly flexible APIs that confuse rather than help

## Proposed Refactoring

### Core API (Following `pwr` Package Pattern)

**Essential Functions to Keep:**
1. `power_ppi_mean()` - Analytical power for PPI mean
2. `power_ppi_ols()` - Analytical power for PPI OLS  
3. `power_ppi_pp_mean()` - Analytical power for PPI++ mean
4. `power_ppi_pp_ols()` - Analytical power for PPI++ OLS
5. `n_required_pp()` - Sample size for PPI (unified)
6. `n_required_ppi_pp()` - Sample size for PPI++
7. `ppi_ols_fit()` - Fit PPI OLS (for inference)
8. `ppi_plus_ols_fit()` - Fit PPI++ OLS (for inference)

**Optional but Useful:**
9. `plot_type1_error_curve()` - Plotting helper

**Functions to Consolidate/Remove:**
- Consolidate `power_curve_*` functions into one flexible function
- Consolidate `type1_error_curve_*` functions into one flexible function
- Consolidate `simulate_power_*` functions (or move to simulations folder)
- Remove thin wrappers (`power_curve_gaussian`, `power_curve_binomial`)

### Proposed Structure

```
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
```

### Key Principles

1. **Simplicity**: One clear way to do each task
2. **Separation**: Core power calculations separate from simulation code
3. **Consistency**: Follow `pwr` package naming and structure patterns
4. **Minimal API**: Only export what users need
5. **Internal Helpers**: Keep utility functions internal unless truly needed

## Next Steps

1. Create `inst/simulations/` folder
2. Move simulation code to simulations folder
3. Consolidate similar functions
4. Update NAMESPACE to only export essential functions
5. Update documentation
6. Update tests

