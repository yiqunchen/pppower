# 2x2 Contingency Table Calculator

This calculator uses **prediction-powered inference (PPI/`PPI++`)**
formulas from this package for planning studies that compare two groups
via a 2x2 contingency table. It supports both **odds ratio** (OR) and
**relative risk** (RR) effect measures.

Interactive PPI Tool

## 2×2 contingency table calculator

Compare two groups using odds ratio or relative risk with a binary PPI++
surrogate.

**Setup:** Two groups (exposed/control) with a binary outcome Y and a
binary surrogate prediction f with known sensitivity and specificity.

**Odds ratio:** Tested via logistic regression with PPI++ sandwich
variance. Equivalent to a Wald test for log(OR).

**Relative risk:** Tested via delta method on per-group PPI++ mean
estimates. Wald test for log(RR).

**Methods:** details and derivations are on [Variance &
Methods](https://yiqunchen.github.io/pppower/articles/deep-dive-math.md).

Odds Ratio (OR)

Relative Risk (RR)

### Inputs

Group probabilities

Control P(Y=1) (%)

Exposed P(Y=1) (%)

Prevalence of exposure (%)

Unlabeled N (total)

Alpha (two-sided)

Target power (%)

Expected non-response (%)

Surrogate quality

Sensitivity

Specificity

### Results

PPI++ n (total labeled) **-**

PPI++ recruitment n (adj.) **-**

Classical n (no surrogate) **-**

Labeled saving vs classical **-**

True OR **-**

Oracle lambda **-**

Classic

-

PPI++

-

PPI adj.

-

Formula note will appear here.

### Inputs

Group probabilities

Control P(Y=1) (%)

Exposed P(Y=1) (%)

Prevalence of exposure (%)

Unlabeled N (total)

Alpha (two-sided)

Target power (%)

Expected non-response (%)

Surrogate quality

Sensitivity

Specificity

### Results

PPI++ n (total labeled) **-**

PPI++ recruitment n (adj.) **-**

Classical n (no surrogate) **-**

Labeled saving vs classical **-**

True RR **-**

log(RR) **-**

Classic

-

PPI++

-

PPI adj.

-

Formula note will appear here.

Uses the same formulas as
[`power_ppi_2x2()`](https://yiqunchen.github.io/pppower/reference/power_ppi_2x2.md)
in the pppower R package. See [Variance &
Methods](https://yiqunchen.github.io/pppower/articles/deep-dive-math.md)
for derivations.
