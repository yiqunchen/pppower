"""Cross-check pppower's labeled sample size against ppi_py's power analysis.

ppi_py (Angelopoulos et al.) implements the budget-based power analysis of
Broska, Howes and van Loon (2025) in ppi_py/ppi_power_analysis.py. With
predictions nearly free and the total pool capped at n + N, its "cheapest
pair" solves effective_n(n, N) = n_classical, which is the same equation that
pppower solves for fixed N.

Usage (from the repo root):
  Rscript simulation_studies/external_checks/make_pppower_nstar.R
  python3 simulation_studies/external_checks/check_ppi_py.py
"""
import csv, re, urllib.request, pathlib

here = pathlib.Path(__file__).parent
src_path = here / "ppi_power_analysis.py"
if not src_path.exists():
    url = "https://raw.githubusercontent.com/aangelopoulos/ppi_py/main/ppi_py/ppi_power_analysis.py"
    src_path.write_text(urllib.request.urlopen(url).read().decode())

# load only the pure power-analysis functions (drop package-relative and sklearn imports)
src = src_path.read_text()
src = re.sub(r"^from \..*$", "", src, flags=re.M)
src = re.sub(r"^from sklearn.*$", "", src, flags=re.M)
ns = {}
exec(src, ns)

rows = list(csv.DictReader(open(here / "pppower_nstar.csv")))
maxdiff = 0
print(f"{'rho':>4} {'N':>6} {'delta':>5} {'n_cl':>5} {'pppower':>8} {'ppi_py':>7}")
for r in rows:
    rho, N, d = float(r["rho"]), int(r["N"]), float(r["delta"])
    n_cl, n_star = int(r["n_cl"]), int(r["n_star"])
    res = ns["ppi_power"](cost_X=0.0, cost_Y=1.0, cost_Yhat=1e-9, ppi_corr=rho,
                          effective_n=n_cl, n_max=n_star + N)
    maxdiff = max(maxdiff, abs(res["n"] - n_star))
    print(f"{rho:4.1f} {N:6d} {d:5.2f} {n_cl:5d} {n_star:8d} {res['n']:7d}")
print("max |pppower - ppi_py| =", maxdiff, "label(s)")
