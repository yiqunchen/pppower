# Research Simulation Scripts (Legacy)

This directory contains comprehensive research simulation scripts used for the
paper's extensive experiments. These scripts include many configuration options
for comparing different estimator variants, cross-fitting schemes, and model types.

**For simple power analysis validation, see `simulation_studies/` instead.**

## Files

- `simulation_mean.R` - Comprehensive MC experiments for mean estimation
  (includes cross-fitting, external lambda, multiple model types, grid searches)
- `simulation_ols.R` - MC experiments for OLS regression

## When to Use This vs simulation_studies/

| Use Case | Folder |
|----------|--------|
| Validate analytical power formula | `simulation_studies/` |
| Simple power curve generation | `simulation_studies/` |
| Reproduce paper's comprehensive experiments | `inst/simulations/` |
| Compare cross-fitting vs non-crossfitting | `inst/simulations/` |
| Test multiple model types (GLM, RF) | `inst/simulations/` |

## Note

These scripts are **not part of the package API**. For power analysis, use:

**Analytical power:**
- `power_ppi_mean()`, `power_ppi_pp_mean()`
- `power_ppi_ols()`, `power_ppi_pp_ols()`

**Monte Carlo verification:**
- `simulate_power_ppi_mean()`, `simulate_power_ppiplus_mean()`
