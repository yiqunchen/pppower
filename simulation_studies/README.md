# Simulation Studies for PPPower

This folder contains simulation scripts for validating the analytical power
formulas in the pppower package. These are **research scripts**, not part
of the package API.

## Structure

```
simulation_studies/
├── README.md                      # This file
├── sim_utils.R                    # Simple DGP and MC utility functions
├── power_validation_mean.R        # Validate mean estimation power formulas
└── power_validation_ttest.R       # Validate two-sample t-test power formulas
```

## What Power Analysis Validation Needs

A power analysis validation is straightforward:

1. **Generate synthetic data** with known population parameters
2. **Compute estimator** (PPI or PPI++ mean)
3. **Perform test** at level alpha
4. **Repeat R times** and compute rejection rate (empirical power)
5. **Compare** to analytical power formula

That's it. No need for complex cross-fitting schemes, multiple model types,
or grid searches over dozens of parameters.

## Usage

```r
source("simulation_studies/sim_utils.R")
source("simulation_studies/power_validation_mean.R")

# Run a simple validation
result <- validate_ppi_mean_power(
  n = 100,           # labeled sample size
  N = 1000,          # unlabeled sample size
  delta = 0.3,       # effect size
  sigma_y = 1,       # outcome SD
  rho = 0.7,         # correlation between Y and f
  R = 500            # MC replicates
)
```

## Why This is Simpler

The original `inst/simulations/` code was designed for comprehensive
research experiments comparing many configurations. This folder provides
simplified scripts focused on **validating the analytical power formulas**
in the paper.

Key simplifications:
- No cross-fitting (assumes externally trained predictor)
- No model fitting within simulation (predictor is oracle)
- No tidyverse dependencies
- Clear separation of DGP, estimation, and reporting
