# Refactoring Summary

## What I've Done So Far

### ✅ Completed
1. **Created diagnosis document** (`DIAGNOSIS.md`) - Analyzed current package structure
2. **Created refactoring plan** (`REFACTORING_PLAN.md`) - Proposed simplification strategy
3. **Created `inst/simulations/` folder** - Separated research code from package code
4. **Moved research scripts** - `simulation_mean.R` and `simulation_ols.R` moved to `inst/simulations/`
5. **Created README** - Documented what the simulations folder contains

## Current Package State

### Exported Functions: 27 (Too Many!)

**Core Power Functions (8 - Keep):**
- `power_ppi_mean()` ✓
- `power_ppi_ols()` ✓
- `power_ppi_pp_mean()` ✓
- `power_ppi_pp_ols()` ✓
- `n_required_pp()` ✓
- `n_required_ppi_pp()` ✓
- `ppi_ols_fit()` ✓
- `ppi_plus_ols_fit()` ✓

**Plotting Functions (7 - Consolidate):**
- `power_curve_mean()` → Consolidate
- `power_curve_mean_dgp()` → Consolidate
- `power_curve_gaussian()` → **REMOVE** (thin wrapper)
- `power_curve_binomial()` → **REMOVE** (thin wrapper)
- `type1_error_curve_mean()` → Consolidate
- `type1_error_curve_mean_dgp()` → Consolidate
- `plot_type1_error_curve()` ✓ (Keep)

**Simulation Functions (5 - Make Internal):**
- `simulate_power()` → Make internal
- `simulate_power_ppi_mean()` → Make internal
- `simulate_power_ppiplus_mean()` → Make internal
- `simulate_power_ppi_ols()` → Make internal
- `simulate_crossfit_data()` → Make internal

**Other (7):**
- `ppi_ols_wald()` → Make internal
- Various other functions

## Recommended Next Steps

### Option 1: Conservative (Recommended)
1. Remove thin wrappers (`power_curve_gaussian`, `power_curve_binomial`)
2. Make simulation functions internal (not exported)
3. Keep plotting functions as-is for now (can consolidate later)
4. **Result: ~15 exported functions** (down from 27)

### Option 2: Moderate
1. Everything in Option 1, plus:
2. Consolidate `power_curve_mean()` and `power_curve_mean_dgp()` into one function
3. Consolidate `type1_error_curve_mean()` and `type1_error_curve_mean_dgp()` into one function
4. **Result: ~10 exported functions**

### Option 3: Aggressive
1. Everything in Option 2, plus:
2. Consolidate all power functions into unified interface
3. **Result: ~6-8 exported functions** (like `pwr` package)

## My Recommendation

Start with **Option 1** - it's the safest and still reduces complexity significantly. We can always consolidate more later if needed.

The key improvements:
- ✅ Research code separated from package code
- ✅ Thin wrappers removed
- ✅ Simulation functions made internal (users can still access via `:::`)
- ✅ Clear separation of concerns

Would you like me to proceed with Option 1, or do you prefer a different approach?

