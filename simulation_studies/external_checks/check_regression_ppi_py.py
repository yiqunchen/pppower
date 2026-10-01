"""Compare pppower::power_ppi_regression with ppi_py's ppi_ols_power and
ppi_logistic_power on the same pilot data, for one regression coefficient.

Requires the ppi_py package (pip install ppi-python).
Usage (from the repo root):
  Rscript simulation_studies/external_checks/make_regression_pilot.R
  python3 simulation_studies/external_checks/check_regression_ppi_py.py
"""
import csv, pathlib
import numpy as np
from ppi_py.ppi_power_analysis import ppi_ols_power, ppi_logistic_power

here = pathlib.Path(__file__).parent
rows = list(csv.DictReader(open(here / "pppower_regression_nstar.csv")))
print(f"{'model':>12} {'quality':>7} {'N':>6} {'n_cl':>6} {'pppower':>8} {'ppi_py':>7} {'ppi_corr':>8}")
maxdiff = 0
cache = {}
for r in rows:
    tag = r["tag"]
    if tag not in cache:
        d = np.genfromtxt(here / f"pilot_{tag}.csv", delimiter=",", names=True)
        X = np.column_stack([np.ones(len(d)), d["x1"], d["x2"]])
        cache[tag] = (X, d["Y"], d["f"])
    X, Y, f = cache[tag]
    fun = ppi_ols_power if r["type"] == "ols" else ppi_logistic_power
    n_star, N, n_cl = int(r["n_star"]), int(r["N"]), int(r["n_cl"])
    res = fun(X, Y, f, cost_X=0.0, cost_Y=1.0, cost_Yhat=1e-9, coord=1,
              effective_n=n_cl, n_max=n_star + N)
    maxdiff = max(maxdiff, abs(res["n"] - n_star))
    print(f"{r['type']:>12} {float(r['quality']):7.1f} {N:6d} {n_cl:6d} {n_star:8d} {res['n']:7d} {res['ppi_corr']:8.3f}")
print("max |pppower - ppi_py| =", maxdiff, "label(s)")
