# Refactoring Plan: Simplifying pppower Package

## Goal
Simplify the package to follow the `pwr` package model: clean, focused API with minimal exported functions.

## Current Problems
1. **27 exported functions** - way too many
2. **Simulation code mixed with core code** - should be separated
3. **Redundant wrapper functions** - multiple ways to do the same thing
4. **Inconsistent naming** - `_ppi_`, `_pp_`, etc.
5. **Premature optimization** - overly complex APIs

## Proposed Core API (8-10 functions)

### Essential Functions (Keep & Export)
1. `power_ppi_mean()` - Analytical power for PPI mean ✓
2. `power_ppi_ols()` - Analytical power for PPI OLS ✓
3. `power_ppi_pp_mean()` - Analytical power for PPI++ mean ✓
4. `power_ppi_pp_ols()` - Analytical power for PPI++ OLS ✓
5. `n_required_pp()` - Sample size for PPI (unified) ✓
6. `n_required_ppi_pp()` - Sample size for PPI++ ✓
7. `ppi_ols_fit()` - Fit PPI OLS (for inference) ✓
8. `ppi_plus_ols_fit()` - Fit PPI++ OLS (for inference) ✓

### Optional Helper Functions
9. `power_curve()` - Unified power curve function (consolidates all `power_curve_*` functions)
10. `type1_curve()` - Unified type I error curve (consolidates all `type1_error_curve_*` functions)

### Functions to Make Internal (Not Export)
- `simulate_power()` - Keep internal, used by plotting functions
- `simulate_power_ppi_mean()` - Keep internal
- `simulate_power_ppiplus_mean()` - Keep internal
- `simulate_power_ppi_ols()` - Keep internal
- `simulate_power_ppi_pp_ols()` - Keep internal
- `simulate_crossfit_data()` - Keep internal or move to inst/simulations

### Functions to Remove/Consolidate
- `power_curve_mean()` → consolidate into `power_curve()`
- `power_curve_mean_dgp()` → consolidate into `power_curve()`
- `power_curve_gaussian()` → REMOVE (thin wrapper)
- `power_curve_binomial()` → REMOVE (thin wrapper)
- `type1_error_curve_mean()` → consolidate into `type1_curve()`
- `type1_error_curve_mean_dgp()` → consolidate into `type1_curve()`
- `ppi_ols_wald()` → Make internal helper

### Research Code to Move
- `R/simulation_mean.R` → `inst/simulations/simulation_mean.R`
- `R/simulation_ols.R` → `inst/simulations/simulation_ols.R`

## New File Structure

```
R/
  ├── power-mean.R          # power_ppi_mean(), power_ppi_pp_mean()
  ├── power-ols.R           # power_ppi_ols(), power_ppi_pp_ols()
  ├── sample-size.R         # n_required_pp(), n_required_ppi_pp()
  ├── fit-ols.R             # ppi_ols_fit(), ppi_plus_ols_fit()
  ├── plotting.R             # power_curve(), type1_curve(), plot_type1_error_curve()
  ├── simulation-internal.R  # Internal simulation helpers (not exported)
  ├── utils.R                # Internal helpers (resolve_ppi_variances, etc.)
  └── pppower-package.R      # Package documentation

inst/
  └── simulations/           # Research/analysis scripts (not part of package)
      ├── simulation_mean.R
      ├── simulation_ols.R
      └── README.md
```

## Implementation Steps

1. ✅ Create `inst/simulations/` folder
2. Move research code to `inst/simulations/`
3. Consolidate `power_curve_*` functions into single `power_curve()`
4. Consolidate `type1_error_curve_*` functions into single `type1_curve()`
5. Update NAMESPACE to only export essential functions
6. Update documentation
7. Update tests

## Key Principles

1. **Simplicity**: One clear way to do each task
2. **Separation**: Core power calculations separate from simulation code
3. **Consistency**: Follow `pwr` package naming and structure patterns
4. **Minimal API**: Only export what users need
5. **Internal Helpers**: Keep utility functions internal unless truly needed

