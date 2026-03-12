#!/usr/bin/env bash
# Run all power validation simulations (Settings 1-13)
set -euo pipefail

cd /Users/yiqun/Desktop/chen-lab/pppower

echo "=========================================="
echo "pppower: Power Validation Simulations"
echo "Started: $(date)"
echo "=========================================="

# Install dependencies + pppower from local source
echo "Installing dependencies..."
Rscript -e 'if (!requireNamespace("pwr", quietly=TRUE)) install.packages("pwr", repos="https://cloud.r-project.org")'
Rscript -e "devtools::install(quick = TRUE, quiet = TRUE)"

Rscript simulation_studies/power_validation_comprehensive.R 2>&1 | tee simulation_studies/sim_output_$(date +%Y%m%d_%H%M%S).log

echo ""
echo "=========================================="
echo "Finished: $(date)"
echo "=========================================="
